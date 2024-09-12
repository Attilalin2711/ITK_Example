/*=========================================================================
 *
 *  Copyright NumFOCUS
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

// Software Guide : BeginLatex
//
// This example illustrates the use of the
// \code{RayCastInterpolateImageFunction} class to generate digitally
// reconstructed radiographs (DRRs) from a 3D image volume such as CT
// or MR.
//
// Software Guide : EndLatex
#include <iostream>
#include <experimental/filesystem>


#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkPNGImageIO.h"
#include "itkPNGImageIOFactory.h"
#include "itkResampleImageFilter.h"
#include "itkCenteredEuler3DTransform.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkRescaleIntensityImageFilter.h"
#include <itkImageSeriesReader.h>
#include <itkGDCMImageIO.h>
#include <itkGDCMSeriesFileNames.h>
#include "itkRegionOfInterestImageFilter.h"
#include <itkVTKImageToImageFilter.h>
#include <itkImageToVTKImageFilter.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkImageResliceMapper.h>
#include <vtkImageSliceMapper.h>
#include <vtkImageSlice.h>
#include <vtkImageProperty.h>
#include <vtkInteractorStyleImage.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkImageActor.h>
#include <vtkImageStack.h>

// Software Guide : BeginLatex
//
// The \code{RayCastInterpolateImageFunction} class definition for
// this example is contained in the following header file.
//
// Software Guide : EndLatex

// Software Guide : BeginCodeSnippet
#include "itkRayCastInterpolateImageFunction.h"
// Software Guide : EndCodeSnippet

//#define WRITE_CUBE_IMAGE_TO_FILE


  //參數設置
typedef signed short InputPixelType;
typedef unsigned char OutputPixelType;
const unsigned int  Dim = 3;  //數據dimension

typedef itk::Image<InputPixelType, Dim> InputImageType;
typedef itk::Image<OutputPixelType, Dim> OutputImageType;
typedef itk::ImageSeriesReader<InputImageType> ReaderType;

//void usage()
//{
//  std::cerr << "\n";
//  std::cerr << "Usage: DRR <options> [input]\n";
//  std::cerr << "  calculates the Digitally Reconstructed Radiograph from a "
//               "volume. \n\n";
//  std::cerr << " where <options> is one or more of the following:\n\n";
//  std::cerr << "  <-h>                    Display (this) usage information\n";
//  std::cerr << "  <-v>                    Verbose output [default: no]\n";
//  std::cerr << "  <-res float float>      Pixel spacing of the output image "
//               "[default: "
//               "1x1mm]  \n";
//  std::cerr << "  <-size int int>         Dimension of the output image "
//               "[default: 501x501]  \n";
//  std::cerr
//    << "  <-sid float>            Distance of ray source (focal point) "
//       "[default: 400mm]\n";
//  std::cerr
//    << "  <-t float float float>  Translation parameter of the camera \n";
//  std::cerr
//    << "  <-rx float>             Rotation around x,y,z axis in degrees \n";
//  std::cerr << "  <-ry float>\n";
//  std::cerr << "  <-rz float>\n";
//  std::cerr << "  <-normal float float>   The 2D projection normal position "
//               "[default: 0x0mm]\n";
//  std::cerr
//    << "  <-cor float float float> The centre of rotation relative to centre "
//       "of volume\n";
//  std::cerr << "  <-threshold float>      Threshold [default: 0]\n";
//  std::cerr << "  <-o file>               Output image filename\n\n";
//  std::cerr << "                          by  thomas@hartkens.de\n";
//  std::cerr << "                          and john.hipwell@kcl.ac.uk (CISG "
//               "London)\n\n";
//  exit(1);
//}

//bool regionOfInterestImageFilter(InputImageType* image, InputImageType* outImage)
//{
//    InputImageType::IndexType desiredStart;
//    desiredStart.SetElement(0, 0);
//    desiredStart.SetElement(1, 0);
//    desiredStart.SetElement(2, 256);
//
//    InputImageType::SizeType desiredSize;
//    desiredSize.SetElement(0, 512);
//    desiredSize.SetElement(1, 512);
//    desiredSize.SetElement(2, 256);
//
//    InputImageType::RegionType desiredRegion;
//    desiredRegion.SetSize(desiredSize);
//    desiredRegion.SetIndex(desiredStart);
//    // Software Guide : EndCodeSnippet
//
//    typedef itk::RegionOfInterestImageFilter<InputImageType, InputImageType> RegionOfInterestFilterType;
//    typename RegionOfInterestFilterType::Pointer regionInterestFilter = RegionOfInterestFilterType::New();
//    regionInterestFilter->SetInput(image);
//    regionInterestFilter->SetRegionOfInterest(desiredRegion);
//    try
//    {
//        regionInterestFilter->Update();
//    }
//    catch (itk::ExceptionObject& ex)
//    {
//        //读取过程发生错误
//        std::cerr << "Error: " << ex << std::endl;
//        return false;
//    }
//
//    outImage = regionInterestFilter->GetOutput();
//    return true;
//}

int main(int argc, char * argv[])
{
  std::string seriesDicomDir = "D://Image//DRR3d//C396";
  std::string outputName = "D://InsightToolkit-5.2.1//example//DRR//bin//output//drrROI2.png";
  //確認不為空
  if (std::experimental::filesystem::exists(seriesDicomDir) == false ||
       std::experimental::filesystem::is_empty(seriesDicomDir))
  {
      std::cerr << "Error: Dicom Directory ERROR!Please Chick！" << std::endl;
  }


  bool ok;
  bool verbose = true;

  float rx = 90.;
  float ry = 0.;
  float rz = 180.;

  float tx = 0.;
  float ty = 0.;
  float tz = 0.;

  float cx = 0.;
  float cy = 0.;
  float cz = 0.;

  float sid = 400.;

  float sx = 1.0;
  float sy = 1.0;

  int dx = 512;
  int dy = 512;

  float o2Dx = 0;
  float o2Dy = 0;

  double threshold = -100;

  //dicom讀取
  itk::GDCMImageIO::Pointer gdcmIO = itk::GDCMImageIO::New();
  itk::GDCMSeriesFileNames::Pointer seriesFileNames = itk::GDCMSeriesFileNames::New();

  //獲取DICOM系列文件的文件名
  seriesFileNames->SetDirectory(seriesDicomDir);

  //如果文件夹中包含多个序列，使用GetSeriesUIDs获取所有的UID值，再在GetFileNames中指定需要返回的序列的UID，获取该序列的所有文件名
  //FileNamesContainer和SeriesUIDContainerType的定义为：std::vector<std::string> 
  //std::vector<std::string> seriesUIDs = seriesFileNames->GetSeriesUIDs();
  //std::vector<std::string> filenames = seriesFileNames->GetFileNames(seriesUIDs[0])
  const itk::GDCMSeriesFileNames::SeriesUIDContainerType& seriesUIDs = seriesFileNames->GetSeriesUIDs();
  const ReaderType::FileNamesContainer& filenames = seriesFileNames->GetFileNames(seriesUIDs[0]);

  typename ReaderType::Pointer reader = ReaderType::New();
  InputImageType::Pointer image;
  try
  {
      reader->SetImageIO(gdcmIO);
      reader->SetFileNames(filenames);
      reader->Update();
  }
  catch (itk::ExceptionObject& ex)
  {
      std::cerr << "Error: " << ex << std::endl;
      return false;
  }
  image = reader->GetOutput();

  //GET ROI VOLUME
  typedef itk::Image< InputPixelType, 3 > SCALAR_3D_IMAGE;
  InputImageType::IndexType desiredStart;
  desiredStart.SetElement(0, 0);
  desiredStart.SetElement(1, 0);
  desiredStart.SetElement(2, 0);

  InputImageType::SizeType desiredSize;
  desiredSize.SetElement(0, 512);
  desiredSize.SetElement(1, 512);
  desiredSize.SetElement(2, 512);

  InputImageType::RegionType desiredRegion(desiredStart, desiredSize);
  typedef itk::RegionOfInterestImageFilter<InputImageType, InputImageType> RegionOfInterestFilterType;
  typename RegionOfInterestFilterType::Pointer regionInterestFilter = RegionOfInterestFilterType::New();
  regionInterestFilter->SetInput(image);
  regionInterestFilter->SetRegionOfInterest(desiredRegion);
  regionInterestFilter->Update();
  image = regionInterestFilter->GetOutput();

  // Print out the details of the input volume
  if (verbose)
  {
    unsigned int i;
    const InputImageType::SpacingType spacing = image->GetSpacing();
    std::cout << std::endl << "Input ";

    InputImageType::RegionType region = image->GetBufferedRegion();
    region.Print(std::cout);

    std::cout << "  Resolution: [";

    for (i = 0; i < Dim; i++)
    {
      std::cout << spacing[i];
      if (i < Dim - 1)
        std::cout << ", ";
    }
    std::cout << "]" << std::endl;

    const InputImageType::PointType origin = image->GetOrigin();
    std::cout << "  Origin: [";
    for (i = 0; i < Dim; i++)
    {
      std::cout << origin[i];
      if (i < Dim - 1)
        std::cout << ", ";
    }
    std::cout << "]" << std::endl << std::endl;
  }

  // Software Guide : BeginLatex
  //
  // Creation of a \code{ResampleImageFilter} enables coordinates for
  // each of the pixels in the DRR image to be generated. These
  // coordinates are used by the \code{RayCastInterpolateImageFunction}
  // to determine the equation of each corresponding ray which is cast
  // through the input volume.
  //
  // Software Guide : EndLatex

  // Software Guide : BeginCodeSnippet
  using FilterType = itk::ResampleImageFilter<InputImageType, InputImageType>;

  FilterType::Pointer filter = FilterType::New();

  filter->SetInput(image);
  filter->SetDefaultPixelValue(0);
  // Software Guide : EndCodeSnippet

  // Software Guide : BeginLatex
  //
  // An Euler transformation is defined to position the input volume.
  // The \code{ResampleImageFilter} uses this transform to position the
  // output DRR image for the desired view.
  //
  // Software Guide : EndLatex

  // Software Guide : BeginCodeSnippet
  using TransformType = itk::CenteredEuler3DTransform<double>;

  TransformType::Pointer transform = TransformType::New();

  transform->SetComputeZYX(true);

  TransformType::OutputVectorType translation;

  translation[0] = tx;
  translation[1] = ty;
  translation[2] = tz;

  // constant for converting degrees into radians
  const double dtr = (std::atan(1.0) * 4.0) / 180.0;

  transform->SetTranslation(translation);
  transform->SetRotation(dtr * rx, dtr * ry, dtr * rz);

  InputImageType::PointType   imOrigin = image->GetOrigin();
  InputImageType::SpacingType imRes = image->GetSpacing();

  using InputImageRegionType = InputImageType::RegionType;
  using InputImageSizeType = InputImageRegionType::SizeType;

  InputImageRegionType imRegion = image->GetBufferedRegion();
  InputImageSizeType   imSize = imRegion.GetSize();

  imOrigin[0] += imRes[0] * static_cast<double>(imSize[0]) / 2.0;
  imOrigin[1] += imRes[1] * static_cast<double>(imSize[1]) / 2.0;
  imOrigin[2] += imRes[2] * static_cast<double>(imSize[2]) / 2.0;

  TransformType::InputPointType center;
  center[0] = cx + imOrigin[0];
  center[1] = cy + imOrigin[1];
  center[2] = cz + imOrigin[2];

  transform->SetCenter(center);

  if (verbose)
  {
    std::cout << "Image size: " << imSize[0] << ", " << imSize[1] << ", "
              << imSize[2] << std::endl
              << "   resolution: " << imRes[0] << ", " << imRes[1] << ", "
              << imRes[2] << std::endl
              << "   origin: " << imOrigin[0] << ", " << imOrigin[1] << ", "
              << imOrigin[2] << std::endl
              << "   center: " << center[0] << ", " << center[1] << ", "
              << center[2] << std::endl
              << "Transform: " << transform << std::endl;
  }
  // Software Guide : EndCodeSnippet

  // Software Guide : BeginLatex
  //
  // The \code{RayCastInterpolateImageFunction} is instantiated and passed the
  // transform object. The \code{RayCastInterpolateImageFunction} uses this
  // transform to reposition the x-ray source such that the DRR image
  // and x-ray source move as one around the input volume. This coupling
  // mimics the rigid geometry of the x-ray gantry.
  //
  // Software Guide : EndLatex

  // Software Guide : BeginCodeSnippet
  using InterpolatorType = itk::RayCastInterpolateImageFunction<InputImageType, double>;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetTransform(transform);
  // Software Guide : EndCodeSnippet

  // Software Guide : BeginLatex
  //
  // We can then specify a threshold above which the volume's
  // intensities will be integrated.
  //
  // Software Guide : EndLatex

  // Software Guide : BeginCodeSnippet
  interpolator->SetThreshold(threshold);
  // Software Guide : EndCodeSnippet

  // Software Guide : BeginLatex
  //
  // The ray-cast interpolator needs to know the initial position of the
  // ray source or focal point. In this example we place the input
  // volume at the origin and halfway between the ray source and the
  // screen. The distance between the ray source and the screen
  // is the "source to image distance" \code{sid} and is specified by
  // the user.
  //
  // Software Guide : EndLatex

  // Software Guide : BeginCodeSnippet
  InterpolatorType::InputPointType focalpoint;

  focalpoint[0] = imOrigin[0];
  focalpoint[1] = imOrigin[1];
  focalpoint[2] = imOrigin[2] - sid / 2.;

  interpolator->SetFocalPoint(focalpoint);
  // Software Guide : EndCodeSnippet

  if (verbose)
  {
    std::cout << "Focal Point: " << focalpoint[0] << ", " << focalpoint[1]
              << ", " << focalpoint[2] << std::endl;
  }

  // Software Guide : BeginLatex
  //
  // Having initialised the interpolator we pass the object to the
  // resample filter.
  //
  // Software Guide : EndLatex

  // Software Guide : BeginCodeSnippet
  interpolator->Print(std::cout);

  filter->SetInterpolator(interpolator);
  filter->SetTransform(transform);
  // Software Guide : EndCodeSnippet

  // Software Guide : BeginLatex
  //
  // The size and resolution of the output DRR image is specified via the
  // resample filter.
  //
  // Software Guide : EndLatex

  // Software Guide : BeginCodeSnippet

  // setup the scene
  InputImageType::SizeType size;

  size[0] = dx; // number of pixels along X of the 2D DRR image
  size[1] = dy; // number of pixels along Y of the 2D DRR image
  size[2] = 1;  // only one slice

  filter->SetSize(size);

  InputImageType::SpacingType spacing;

  spacing[0] = sx;  // pixel spacing along X of the 2D DRR image [mm]
  spacing[1] = sy;  // pixel spacing along Y of the 2D DRR image [mm]
  spacing[2] = 1.0; // slice thickness of the 2D DRR image [mm]

  filter->SetOutputSpacing(spacing);

  // Software Guide : EndCodeSnippet

  if (verbose)
  {
    std::cout << "Output image size: " << size[0] << ", " << size[1] << ", "
              << size[2] << std::endl;

    std::cout << "Output image spacing: " << spacing[0] << ", " << spacing[1]
              << ", " << spacing[2] << std::endl;
  }

  // Software Guide : BeginLatex
  //
  // In addition the position of the DRR is specified. The default
  // position of the input volume, prior to its transformation is
  // half-way between the ray source and screen and unless specified
  // otherwise the normal from the "screen" to the ray source passes
  // directly through the centre of the DRR.
  //
  // Software Guide : EndLatex

  // Software Guide : BeginCodeSnippet

  double origin[Dim];

  origin[0] = imOrigin[0] + o2Dx - sx * ((double)dx - 1.) / 2.;
  origin[1] = imOrigin[1] + o2Dy - sy * ((double)dy - 1.) / 2.;
  origin[2] = imOrigin[2] + sid / 2.;

  filter->SetOutputOrigin(origin);
  // Software Guide : EndCodeSnippet

  if (verbose)
  {
    std::cout << "Output image origin: " << origin[0] << ", " << origin[1]
              << ", " << origin[2] << std::endl;
  }

  // create writer

    // Software Guide : BeginLatex
    //
    // The output of the resample filter can then be passed to a writer to
    // save the DRR image to a file.
    //
    // Software Guide : EndLatex

    // Software Guide : BeginCodeSnippet
  using RescaleFilterType = itk::RescaleIntensityImageFilter<InputImageType, OutputImageType>;
  RescaleFilterType::Pointer rescaler = RescaleFilterType::New();
  rescaler->SetOutputMinimum(0);
  rescaler->SetOutputMaximum(255);
  rescaler->SetInput(filter->GetOutput());

  /*
  using WriterType = itk::ImageFileWriter<OutputImageType>;
  WriterType::Pointer writer = WriterType::New();

  using pngType = itk::PNGImageIO;
  pngType::Pointer pngio1 = pngType::New();
  itk::PNGImageIOFactory::RegisterOneFactory();
  writer->SetFileName(outputName);
  writer->SetImageIO(pngio1);
  writer->SetInput(rescaler->GetOutput());

  try
  {
    std::cout << "Writing image: " << outputName << std::endl;
    writer->Update();
  }
  catch (const itk::ExceptionObject & err)
  {
    std::cerr << "ERROR: ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
   
    // Software Guide : EndCodeSnippet
  }
  */

  typedef itk::ImageToVTKImageFilter<OutputImageType> ConnectorType;
  ConnectorType::Pointer connector = ConnectorType::New();
  connector->SetInput(rescaler->GetOutput());
  try
  {
      connector->Update();
  }
  catch (itk::ExceptionObject& e)
  {
      std::cerr << "exception in file reader" << std::endl;
      std::cerr << e << std::endl;
  }

  vtkImageData* img = vtkImageData::New();
  img = connector->GetOutput();
  
 /* vtkImageActor* actor = vtkImageActor::New();
  actor->SetInputData(img);*/
  
  vtkImageResliceMapper* imgresliceMapper = vtkImageResliceMapper::New();
  imgresliceMapper->SetInputData(img);
  //imgresliceMapper->SliceAtFocalPointOn();
  //imgresliceMapper->SliceFacesCameraOn();
  imgresliceMapper->SetSlabThickness(150);
  //imgresliceMapper->SetSlabTypeToMin();

  //vtkImageSliceMapper* sliceMapper = vtkImageSliceMapper::New();
  //sliceMapper->SetInputData(img);

  vtkImageSlice* imageSlice = vtkImageSlice::New();
  imageSlice->SetMapper(imgresliceMapper);
 
  vtkImageStack* imagestack = vtkImageStack::New();
  imagestack->AddImage(imageSlice);
  //imageSlice->GetProperty()->SetInterpolationTypeToNearest();


  // Setup renderers.
  vtkRenderer* renderer = vtkRenderer::New();
  renderer->AddViewProp(imagestack);
  renderer->ResetCamera();


  // Setup render window.
  vtkRenderWindow* renderWindow = vtkRenderWindow::New();
  renderWindow->SetSize(300, 300);
  renderWindow->AddRenderer(renderer);
  renderWindow->SetWindowName("Interpolation");

  // Setup render window interactor.
  vtkRenderWindowInteractor* renderWindowInteractor = vtkRenderWindowInteractor::New();

  vtkInteractorStyleImage* style = vtkInteractorStyleImage::New();

  renderWindowInteractor->SetInteractorStyle(style);

  // Render and start interaction.
  renderWindowInteractor->SetRenderWindow(renderWindow);
  renderWindow->Render();
  renderWindowInteractor->Initialize();

  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}
