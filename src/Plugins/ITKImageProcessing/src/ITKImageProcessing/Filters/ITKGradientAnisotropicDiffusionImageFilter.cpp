#include "ITKGradientAnisotropicDiffusionImageFilter.hpp"

#include "ITKImageProcessing/Common/ITKArrayHelper.hpp"
#include "ITKImageProcessing/Common/sitkCommon.hpp"

#include "simplnx/Parameters/ArraySelectionParameter.hpp"
#include "simplnx/Parameters/DataObjectNameParameter.hpp"
#include "simplnx/Parameters/GeometrySelectionParameter.hpp"
#include "simplnx/Parameters/NumberParameter.hpp"

#include <itkGradientAnisotropicDiffusionImageFilter.h>

using namespace nx::core;

namespace cxITKGradientAnisotropicDiffusionImageFilter
{
using ArrayOptionsType = ITK::FloatingScalarPixelIdTypeList;

struct ITKGradientAnisotropicDiffusionImageFilterFunctor
{
  using IntermediateType = float64;

  float64 timeStep = 0.125;
  float64 conductanceParameter = 3;
  uint32 conductanceScalingUpdateInterval = 1u;
  uint32 numberOfIterations = 5u;

  template <class InputImageT, class OutputImageT, uint32 Dimension>
  auto createFilter() const
  {
    using FilterType = itk::GradientAnisotropicDiffusionImageFilter<InputImageT, OutputImageT>;
    auto filter = FilterType::New();
    filter->SetTimeStep(timeStep);
    filter->SetConductanceParameter(conductanceParameter);
    filter->SetConductanceScalingUpdateInterval(conductanceScalingUpdateInterval);
    filter->SetNumberOfIterations(numberOfIterations);
    return filter;
  }
};
} // namespace cxITKGradientAnisotropicDiffusionImageFilter

namespace nx::core
{
//------------------------------------------------------------------------------
std::string ITKGradientAnisotropicDiffusionImageFilter::name() const
{
  return FilterTraits<ITKGradientAnisotropicDiffusionImageFilter>::name;
}

//------------------------------------------------------------------------------
std::string ITKGradientAnisotropicDiffusionImageFilter::className() const
{
  return FilterTraits<ITKGradientAnisotropicDiffusionImageFilter>::className;
}

//------------------------------------------------------------------------------
Uuid ITKGradientAnisotropicDiffusionImageFilter::uuid() const
{
  return FilterTraits<ITKGradientAnisotropicDiffusionImageFilter>::uuid;
}

//------------------------------------------------------------------------------
std::string ITKGradientAnisotropicDiffusionImageFilter::humanName() const
{
  return "ITK Gradient Anisotropic Diffusion Image Filter";
}

//------------------------------------------------------------------------------
std::vector<std::string> ITKGradientAnisotropicDiffusionImageFilter::defaultTags() const
{
  return {className(), "ITKImageProcessing", "ITKGradientAnisotropicDiffusionImageFilter", "ITKAnisotropicSmoothing", "AnisotropicSmoothing"};
}

//------------------------------------------------------------------------------
Parameters ITKGradientAnisotropicDiffusionImageFilter::parameters() const
{
  Parameters params;
  params.insertSeparator(Parameters::Separator{"Input Parameter(s)"});
  params.insert(std::make_unique<Float64Parameter>(k_TimeStep_Key, "Time Step", "The time step to be used for each iteration.", 0.125));
  params.insert(std::make_unique<Float64Parameter>(k_ConductanceParameter_Key, "Conductance Parameter",
                                                   "The conductance parameter controls the sensitivity of the conductance term in the basic anisotropic diffusion equation", 3.0));
  params.insert(std::make_unique<UInt32Parameter>(k_ConductanceScalingUpdateInterval_Key, "Conductance Scaling Update Interval", "The interval between conductance updates.", 1u));
  params.insert(std::make_unique<UInt32Parameter>(k_NumberOfIterations_Key, "Number Of Iterations",
                                                  "Specifies the number of iterations (time-step updates) that the solver will perform to produce a solution image", 5u));

  params.insertSeparator(Parameters::Separator{"Input Cell Data"});
  params.insert(std::make_unique<GeometrySelectionParameter>(k_InputImageGeomPath_Key, "Image Geometry", "Select the Image Geometry Group from the DataStructure.", DataPath({"Image Geometry"}),
                                                             GeometrySelectionParameter::AllowedTypes{IGeometry::Type::Image}));
  params.insert(std::make_unique<ArraySelectionParameter>(k_InputImageDataPath_Key, "Input Cell Data", "The image data that will be processed by this filter.", DataPath{},
                                                          nx::core::ITK::GetFloatingScalarPixelAllowedTypes()));

  params.insertSeparator(Parameters::Separator{"Output Cell Data"});
  params.insert(
      std::make_unique<DataObjectNameParameter>(k_OutputImageArrayName_Key, "Output Image Data Array", "The result of the processing will be stored in this Data Array.", "Output Image Data"));

  return params;
}

//------------------------------------------------------------------------------
IFilter::VersionType ITKGradientAnisotropicDiffusionImageFilter::parametersVersion() const
{
  return 1;
}

//------------------------------------------------------------------------------
IFilter::UniquePointer ITKGradientAnisotropicDiffusionImageFilter::clone() const
{
  return std::make_unique<ITKGradientAnisotropicDiffusionImageFilter>();
}

//------------------------------------------------------------------------------
IFilter::PreflightResult ITKGradientAnisotropicDiffusionImageFilter::preflightImpl(const DataStructure& dataStructure, const Arguments& filterArgs, const MessageHandler& messageHandler,
                                                                                   const std::atomic_bool& shouldCancel) const
{
  auto imageGeomPath = filterArgs.value<DataPath>(k_InputImageGeomPath_Key);
  auto selectedInputArray = filterArgs.value<DataPath>(k_InputImageDataPath_Key);
  auto outputArrayName = filterArgs.value<DataObjectNameParameter::ValueType>(k_OutputImageArrayName_Key);
  auto timeStep = filterArgs.value<float64>(k_TimeStep_Key);
  auto conductanceParameter = filterArgs.value<float64>(k_ConductanceParameter_Key);
  auto conductanceScalingUpdateInterval = filterArgs.value<uint32>(k_ConductanceScalingUpdateInterval_Key);
  auto numberOfIterations = filterArgs.value<uint32>(k_NumberOfIterations_Key);
  const DataPath outputArrayPath = selectedInputArray.replaceName(outputArrayName);

  Result<OutputActions> resultOutputActions = ITK::DataCheck<cxITKGradientAnisotropicDiffusionImageFilter::ArrayOptionsType>(dataStructure, selectedInputArray, imageGeomPath, outputArrayPath);

  return {std::move(resultOutputActions)};
}

//------------------------------------------------------------------------------
Result<> ITKGradientAnisotropicDiffusionImageFilter::executeImpl(DataStructure& dataStructure, const Arguments& filterArgs, const PipelineFilter* pipelineNode, const MessageHandler& messageHandler,
                                                                 const std::atomic_bool& shouldCancel) const
{
  auto imageGeomPath = filterArgs.value<DataPath>(k_InputImageGeomPath_Key);
  auto selectedInputArray = filterArgs.value<DataPath>(k_InputImageDataPath_Key);
  auto outputArrayName = filterArgs.value<DataObjectNameParameter::ValueType>(k_OutputImageArrayName_Key);
  const DataPath outputArrayPath = selectedInputArray.replaceName(outputArrayName);

  auto timeStep = filterArgs.value<float64>(k_TimeStep_Key);
  auto conductanceParameter = filterArgs.value<float64>(k_ConductanceParameter_Key);
  auto conductanceScalingUpdateInterval = filterArgs.value<uint32>(k_ConductanceScalingUpdateInterval_Key);
  auto numberOfIterations = filterArgs.value<uint32>(k_NumberOfIterations_Key);

  const cxITKGradientAnisotropicDiffusionImageFilter::ITKGradientAnisotropicDiffusionImageFilterFunctor itkFunctor = {timeStep, conductanceParameter, conductanceScalingUpdateInterval,
                                                                                                                      numberOfIterations};

  auto& imageGeom = dataStructure.getDataRefAs<ImageGeom>(imageGeomPath);

  return ITK::Execute<cxITKGradientAnisotropicDiffusionImageFilter::ArrayOptionsType>(dataStructure, selectedInputArray, imageGeomPath, outputArrayPath, itkFunctor, shouldCancel);
}
} // namespace nx::core
