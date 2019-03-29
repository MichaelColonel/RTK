/*=========================================================================
*
*  Copyright Insight Software Consortium
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

/**
 * Test program for itkCudaImage class.
 * This program shows how to use Cuda image and Cuda program.
 */
#include "itkCudaImage.h"
#include "itkCudaKernelManager.h"
#include "itkCudaContextManager.h"
#include "itkCudaImageOps.h"

using ItkImage1f = itk::CudaImage<float, 2>;


int itkCudaImageTest(int argc, char *argv[])
{
  if (argc > 1)
  {
    std::cout << "received " << argc << " arguments, but didn't expect any."
              << "first ignored argument: " << argv[1] << std::endl;
  }
  unsigned int width, height;

  width  = 256;
  height = 256;

  //
  // create CudaImage
  //

  // set size & region
  ItkImage1f::Pointer srcA, srcB, dest;

  ItkImage1f::IndexType start;
  start[0] = 0;
  start[1] = 0;
  ItkImage1f::SizeType size;
  size[0] = width;
  size[1] = height;
  ItkImage1f::RegionType region;
  region.SetSize( size );
  region.SetIndex( start );

  // create
  srcA = ItkImage1f::New();
  srcA->SetRegions( region );
  srcA->Allocate();
  srcA->FillBuffer( 1.0f );

  srcB = ItkImage1f::New();
  srcB->SetRegions( region );
  srcB->Allocate();
  srcB->FillBuffer( 3.0f );

  dest = ItkImage1f::New();
  dest->SetRegions( region );
  dest->Allocate();
  dest->FillBuffer( 0.0f );

  // check pixel value
  ItkImage1f::IndexType idx;
  idx[0] = 0;
  idx[1] = 0;

  unsigned int nElem = width*height;

  //
  // create Cuda program object
  //
  itk::CudaKernelManager::Pointer kernelManager = itk::CudaKernelManager::New();

  // load program and compile
  std::string CudaSource = itk::CudaImageOps::GetCudaPTXSource();
  kernelManager->LoadProgramFromString(CudaSource.c_str());

  int2 imageSize = make_int2(width, height);

  //
  // create addition kernel
  //
  int kernel_add = kernelManager->CreateKernel("ImageAdd", typeid(float));

  std::cout << "======================" << std::endl;
  std::cout << "Kernel : Addition" << std::endl;
  std::cout << "------------------" << std::endl;
  std::cout << "Before Cuda kernel execution" << std::endl;
  std::cout << "SrcA : " << srcA->GetPixel( idx ) << std::endl;
  std::cout << "SrcB : " << srcB->GetPixel( idx ) << std::endl;
  std::cout << "Dest : " << dest->GetPixel( idx ) << std::endl;

  kernelManager->SetKernelArg(kernel_add, 0, 1, &imageSize);
  kernelManager->SetKernelArgWithImage(kernel_add, 1, srcA);
  kernelManager->SetKernelArgWithImage(kernel_add, 2, srcB);
  kernelManager->SetKernelArgWithImage(kernel_add, 3, dest);
  kernelManager->SetKernelArg(kernel_add, 4, sizeof(unsigned int), &nElem);

  kernelManager->LaunchKernel2D(kernel_add, width, height, 16, 16);

  std::cout << "------------------" << std::endl;
  std::cout << "After Cuda kernel execution" << std::endl;
  std::cout << "SrcA : " << srcA->GetPixel( idx ) << std::endl;
  std::cout << "SrcB : " << srcB->GetPixel( idx ) << std::endl;
  std::cout << "Des  : " << dest->GetPixel( idx ) << std::endl;
  std::cout << "======================" << std::endl;

  //
  // create multiplication kernel
  //
  int kernel_mult = kernelManager->CreateKernel("ImageMult", typeid(float));

  std::cout << "======================" << std::endl;
  std::cout << "Kernel : Multiplication" << std::endl;
  std::cout << "------------------" << std::endl;
  std::cout << "Before Cuda kernel execution" << std::endl;
  std::cout << "SrcA : " << srcA->GetPixel( idx ) << std::endl;
  std::cout << "SrcB : " << srcB->GetPixel( idx ) << std::endl;
  std::cout << "Dest : " << dest->GetPixel( idx ) << std::endl;

  kernelManager->SetKernelArg(kernel_mult, 0, 1, &imageSize);
  kernelManager->SetKernelArgWithImage(kernel_mult, 1, srcA);
  kernelManager->SetKernelArgWithImage(kernel_mult, 2, srcB);
  kernelManager->SetKernelArgWithImage(kernel_mult, 3, dest);
  kernelManager->SetKernelArg(kernel_mult, 4, sizeof(unsigned int), &nElem);
  kernelManager->LaunchKernel2D(kernel_mult, width, height, 16, 16);

  std::cout << "------------------" << std::endl;
  std::cout << "After Cuda kernel execution" << std::endl;
  std::cout << "SrcA : " << srcA->GetPixel( idx ) << std::endl;
  std::cout << "SrcB : " << srcB->GetPixel( idx ) << std::endl;
  std::cout << "Des  : " << dest->GetPixel( idx ) << std::endl;
  std::cout << "======================" << std::endl;

  itk::CudaContextManager *contextManager = itk::CudaContextManager::GetInstance();

  //
  // create subtraction kernel
  //
  int kernel_sub = kernelManager->CreateKernel("ImageSub", typeid(float));

  srcA->FillBuffer( 2.0f );
  srcB->FillBuffer( 4.0f );
  dest->FillBuffer( 1.0f );

  std::cout << "======================" << std::endl;
  std::cout << "Kernel : Subtraction" << std::endl;
  std::cout << "------------------" << std::endl;
  std::cout << "Before Cuda kernel execution" << std::endl;
  std::cout << "SrcA : " << srcA->GetPixel( idx ) << std::endl;
  std::cout << "SrcB : " << srcB->GetPixel( idx ) << std::endl;
  std::cout << "Dest : " << dest->GetPixel( idx ) << std::endl;

  kernelManager->SetKernelArg(kernel_sub, 0, 1, &imageSize);
  kernelManager->SetKernelArgWithImage(kernel_sub, 1, srcA);
  kernelManager->SetKernelArgWithImage(kernel_sub, 2, srcB);
  kernelManager->SetKernelArgWithImage(kernel_sub, 3, dest);
  kernelManager->SetKernelArg(kernel_sub, 4, sizeof(unsigned int), &nElem);
  kernelManager->LaunchKernel2D(kernel_sub, width, height, 16, 16);

  std::cout << "------------------" << std::endl;
  std::cout << "After Cuda kernel execution" << std::endl;
  std::cout << "SrcA : " << srcA->GetPixel( idx ) << std::endl;
  std::cout << "SrcB : " << srcB->GetPixel( idx ) << std::endl;
  std::cout << "Des  : " << dest->GetPixel( idx ) << std::endl;
  std::cout << "======================" << std::endl;

  return EXIT_SUCCESS;
}
