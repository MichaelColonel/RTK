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

#include "rtkwaveletsdenoising_ggo.h"
#include "rtkGgoFunctions.h"

#include "rtkDeconstructSoftThresholdReconstructImageFilter.h"
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <string>
#include <sstream>

int main(int argc, char * argv[])
{
  GGO(rtkwaveletsdenoising, args_info);

  using OutputPixelType = float;
  constexpr unsigned int Dimension = 3;

  using OutputImageType = itk::Image< OutputPixelType, Dimension >;

  // Read the input image
  using ReaderType = itk::ImageFileReader<OutputImageType>;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(args_info.input_arg);

  // Create the denoising filter
  using WaveletsSoftThresholdFilterType = rtk::DeconstructSoftThresholdReconstructImageFilter<OutputImageType>;
  WaveletsSoftThresholdFilterType::Pointer wst = WaveletsSoftThresholdFilterType::New();
  wst->SetInput(reader->GetOutput());
  wst->SetOrder(args_info.order_arg);
  wst->SetThreshold(args_info.threshold_arg);
  wst->SetNumberOfLevels(args_info.level_arg);

  // Write reconstruction
  using WriterType = itk::ImageFileWriter<OutputImageType>;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( wst->GetOutput() );
  writer->SetFileName(args_info.output_arg);
  TRY_AND_EXIT_ON_ITK_EXCEPTION( writer->Update() )

  return EXIT_SUCCESS;
}
