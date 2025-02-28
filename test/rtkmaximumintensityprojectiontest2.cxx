#include "rtkTest.h"
#include "rtkThreeDCircularProjectionGeometryXMLFile.h"
#include "rtkRayBoxIntersectionImageFilter.h"
#include "rtkConstantImageSource.h"
#include "rtkJosephForwardProjectionImageFilter2.h"
#include <itkImageRegionSplitterDirection.h>
#include <itkImageRegionIterator.h>
#include <itkImageFileWriter.h>
#include <itkShiftScaleImageFilter.h>
#include <cmath>


/**
 * \file rtkmaximumintensityprojectiontest2.cxx
 *
 * \brief Functional test for MIP forward filter
 *
 * The test projects a volume filled with ones.
 *
 * \author Mikhail Polkovnikov
 */

int
main(int, char **)
{
  constexpr unsigned int Dimension = 3;
  using OutputPixelType = float;
  using OutputImageType = itk::Image<OutputPixelType, Dimension>;
  constexpr unsigned int NumberOfProjectionImages = 1;

  // Constant image sources
  using ConstantImageSourceType = rtk::ConstantImageSource<OutputImageType>;
  ConstantImageSourceType::PointType   origin;
  ConstantImageSourceType::SizeType    size;
  ConstantImageSourceType::SizeType    sizeOutput;
  ConstantImageSourceType::SpacingType spacing;
  ConstantImageSourceType::SpacingType spacingOutput;

  // Create MIP Forward Projector volume input.
  const ConstantImageSourceType::Pointer volInput = ConstantImageSourceType::New();
  origin[0] = 0;
  origin[1] = 0;
  origin[2] = 0;
  size[0] = 64;
  size[1] = 64;
  size[2] = 64;
  spacing[0] = 4.;
  spacing[1] = 4.;
  spacing[2] = 4.;

  sizeOutput[0] = 200;
  sizeOutput[1] = 200;
  sizeOutput[2] = 200;
  spacingOutput[0] = 1.;
  spacingOutput[1] = 1.;
  spacingOutput[2] = 1.;

  volInput->SetOrigin(origin);
  volInput->SetSpacing(spacing);
  volInput->SetSize(size);
  volInput->SetConstant(1.);
  volInput->UpdateOutputInformation();

  // Initialization Imager Volume
  const ConstantImageSourceType::Pointer projInput = ConstantImageSourceType::New();
  sizeOutput[2] = NumberOfProjectionImages;
  projInput->SetOrigin(origin);
  projInput->SetSpacing(spacingOutput);
  projInput->SetSize(sizeOutput);
  projInput->SetConstant(0.);
  projInput->Update();

  // Joseph Forward Projection filter for MIP projection
  using MIPType = rtk::JosephForwardProjectionImageFilter2<OutputImageType, OutputImageType>;
  MIPType::Pointer mipfp = MIPType::New();
  mipfp->SetInput(projInput->GetOutput());
  mipfp->SetInput(1, volInput->GetOutput());

  // Function to compute the maximum intensity (MIP) value along the ray projection.
  MIPType::SumAlongRayFunc sumAlongFunc = [](const itk::ThreadIdType,
                                             MIPType::OutputPixelType &    mipValue,
                                             const MIPType::InputPixelType volumeValue,
                                             const itk::Vector<double, 3> &) -> void {
    MIPType::OutputPixelType tmp = static_cast<MIPType::OutputPixelType>(volumeValue);
    if (tmp > mipValue)
    {
      mipValue = tmp;
    }
  };
  mipfp->SetSumAlongRay(sumAlongFunc);

  // Performs a MIP forward projection, i.e. calculation of a maximum intensity
  // step along the x-ray line.
  MIPType::ProjectedValueAccumulationFunc projAccumFunc = [](const itk::ThreadIdType          itkNotUsed(threadId),
                                                             const MIPType::InputPixelType &  input,
                                                             MIPType::OutputPixelType &       output,
                                                             const MIPType::OutputPixelType & rayCastValue,
                                                             const itk::Vector<double, 3> &   stepInMM,
                                                             const itk::Vector<double, 3> &   itkNotUsed(source),
                                                             const itk::Vector<double, 3> &   itkNotUsed(sourceToPixel),
                                                             const itk::Vector<double, 3> &   itkNotUsed(nearestPoint),
                                                             const itk::Vector<double, 3> & itkNotUsed(farthestPoint)) {
    MIPType::OutputPixelType tmp = static_cast<MIPType::OutputPixelType>(input);
    if (tmp < rayCastValue)
    {
      tmp = rayCastValue;
    }
    output = tmp * stepInMM.GetNorm();
  };
  mipfp->SetProjectedValueAccumulation(projAccumFunc);

  using GeometryType = rtk::ThreeDCircularProjectionGeometry;
  GeometryType::Pointer geometry = GeometryType::New();
  geometry->AddProjection(700, 800, 0);

  mipfp->SetGeometry(geometry);
  mipfp->Update();

  using ConstIteratorType = itk::ImageRegionConstIterator<OutputImageType>;
  ConstIteratorType inputIt(mipfp->GetOutput(), mipfp->GetOutput()->GetRequestedRegion());

  inputIt.GoToBegin();

  bool res = false;
  while (!inputIt.IsAtEnd())
  {
    OutputPixelType pixel = inputIt.Get();
    if (pixel < 4. || pixel > 4.25)
    {
      res = true;
    }
    ++inputIt;
  }
  if (res)
  {
    return EXIT_FAILURE;
  }

  std::cout << "\n\nTest PASSED! " << std::endl;
  return EXIT_SUCCESS;
}
