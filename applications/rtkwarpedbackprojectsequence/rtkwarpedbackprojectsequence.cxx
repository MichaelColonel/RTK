/*=========================================================================
 *
 *  Copyright RTK Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#include "rtkwarpedbackprojectsequence_ggo.h"
#include "rtkGgoFunctions.h"

#include "rtkWarpProjectionStackToFourDImageFilter.h"
#include "rtkThreeDCircularProjectionGeometryXMLFile.h"
#include "rtkPhasesToInterpolationWeights.h"

#ifdef RTK_USE_CUDA
  #include "itkCudaImage.h"
#endif
#include <itkImageFileWriter.h>

int main(int argc, char * argv[])
{
  GGO(rtkwarpedbackprojectsequence, args_info);

  using OutputPixelType = float;
  using DVFVectorType = itk::CovariantVector< OutputPixelType, 3 >;

#ifdef RTK_USE_CUDA
  using VolumeSeriesType = itk::CudaImage< OutputPixelType, 4 >;
  using ProjectionStackType = itk::CudaImage< OutputPixelType, 3 >;
  using DVFSequenceImageType = itk::CudaImage<DVFVectorType, VolumeSeriesType::ImageDimension>;
  using DVFImageType = itk::CudaImage<DVFVectorType, VolumeSeriesType::ImageDimension - 1>;
#else
  using VolumeSeriesType = itk::Image< OutputPixelType, 4 >;
  using ProjectionStackType = itk::Image< OutputPixelType, 3 >;
  using DVFSequenceImageType = itk::Image<DVFVectorType, VolumeSeriesType::ImageDimension>;
#endif
  using DVFReaderType = itk::ImageFileReader<  DVFSequenceImageType >;

  // Projections reader
  using ReaderType = rtk::ProjectionsReader< ProjectionStackType >;
  ReaderType::Pointer reader = ReaderType::New();
  rtk::SetProjectionsReaderFromGgo<ReaderType, args_info_rtkwarpedbackprojectsequence>(reader, args_info);

  // Geometry
  if(args_info.verbose_flag)
    std::cout << "Reading geometry information from "
              << args_info.geometry_arg
              << "..."
              << std::endl;
  rtk::ThreeDCircularProjectionGeometryXMLFileReader::Pointer geometryReader;
  geometryReader = rtk::ThreeDCircularProjectionGeometryXMLFileReader::New();
  geometryReader->SetFilename(args_info.geometry_arg);
  TRY_AND_EXIT_ON_ITK_EXCEPTION( geometryReader->GenerateOutputInformation() )

  // Create input: either an existing volume read from a file or a blank image
  itk::ImageSource< VolumeSeriesType >::Pointer inputFilter;
  if(args_info.input_given)
    {
    // Read an existing image to initialize the volume
    using InputReaderType = itk::ImageFileReader<  VolumeSeriesType >;
    InputReaderType::Pointer inputReader = InputReaderType::New();
    inputReader->SetFileName( args_info.input_arg );
    inputFilter = inputReader;
    }
  else
    {
    // Create new empty volume sequence
    using ConstantImageSourceType = rtk::ConstantImageSource< VolumeSeriesType >;
    ConstantImageSourceType::Pointer constantImageSource = ConstantImageSourceType::New();
    rtk::SetConstantImageSourceFromGgo<ConstantImageSourceType, args_info_rtkwarpedbackprojectsequence>(constantImageSource, args_info);

    // GenGetOpt can't handle default arguments for multiple arguments like dimension or spacing.
    // The only default it accepts is to set all components of a multiple argument to the same value.
    // Default dimension is 256^4, ie the number of reconstructed instants is 256. It has to be set to a more reasonable value
    // which is why a "frames" argument is introduced
    ConstantImageSourceType::SizeType inputSize = constantImageSource->GetSize();
    inputSize[3] = args_info.frames_arg;
    constantImageSource->SetSize(inputSize);

    inputFilter = constantImageSource;
    }
  TRY_AND_EXIT_ON_ITK_EXCEPTION( inputFilter->Update() )
  inputFilter->ReleaseDataFlagOn();

  // Read the phases file
  rtk::PhasesToInterpolationWeights::Pointer phaseReader = rtk::PhasesToInterpolationWeights::New();
  phaseReader->SetFileName(args_info.signal_arg);
  phaseReader->SetNumberOfReconstructedFrames(inputFilter->GetOutput()->GetLargestPossibleRegion().GetSize(3));
  phaseReader->Update();

  // Create the main filter, connect the basic inputs, and set the basic parameters
  using WarpForwardProjectSequenceFilterType = rtk::WarpProjectionStackToFourDImageFilter<VolumeSeriesType,
                                                     ProjectionStackType>;
  WarpForwardProjectSequenceFilterType::Pointer warpbackprojectsequence = WarpForwardProjectSequenceFilterType::New();
  warpbackprojectsequence->SetInputVolumeSeries(inputFilter->GetOutput() );
  warpbackprojectsequence->SetInputProjectionStack(reader->GetOutput());
  warpbackprojectsequence->SetGeometry( geometryReader->GetOutputObject() );
  warpbackprojectsequence->SetWeights(phaseReader->GetOutput());
  warpbackprojectsequence->SetSignal(rtk::ReadSignalFile(args_info.signal_arg));

  // Read DVF
  DVFReaderType::Pointer dvfReader = DVFReaderType::New();
  dvfReader->SetFileName( args_info.dvf_arg );
  TRY_AND_EXIT_ON_ITK_EXCEPTION( dvfReader->Update() )
  warpbackprojectsequence->SetDisplacementField(dvfReader->GetOutput());

  TRY_AND_EXIT_ON_ITK_EXCEPTION( warpbackprojectsequence->Update() )

  // Write
  using WriterType = itk::ImageFileWriter< VolumeSeriesType >;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( args_info.output_arg );
  writer->SetInput( warpbackprojectsequence->GetOutput() );
  TRY_AND_EXIT_ON_ITK_EXCEPTION( writer->Update() )

  return EXIT_SUCCESS;
}
