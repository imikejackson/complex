#include "ITKMedianProjectionImageFilter.hpp"

#include "simplnx/Common/TypesUtility.hpp"
#include "simplnx/DataStructure/IDataArray.hpp"
#include "simplnx/Parameters/ArraySelectionParameter.hpp"
#include "simplnx/Parameters/BoolParameter.hpp"
#include "simplnx/Parameters/DataObjectNameParameter.hpp"
#include "simplnx/Parameters/GeometrySelectionParameter.hpp"
#include "simplnx/Parameters/NumberParameter.hpp"
#include "simplnx/Parameters/StringParameter.hpp"

#include "ITKImageProcessing/Common/ProjectionUtils.hpp"

#include <itkMedianProjectionImageFilter.h>

using namespace nx::core;

namespace cxITKMedianProjectionImageFilter
{
using ArrayOptionsType = ITK::ScalarPixelIdTypeList;

// Uncommenting below line enables RGB/ARGB images that are currently unsupported
// using ArrayOptionsType = ITK::ArrayOptions<ITK::ArrayComponentOptions<true, false, true>, ITK::ArrayUseAllTypes>;
// VectorPixelIDTypeList;

struct ITKMedianProjectionImageFilterFunctor
{
  uint32 projectionDimension = 0u;

  template <class InputImageT, class OutputImageT, uint32 Dimension>
  auto createFilter() const
  {
    using FilterType = itk::MedianProjectionImageFilter<InputImageT, OutputImageT>;
    auto filter = FilterType::New();
    filter->SetProjectionDimension(projectionDimension);
    return filter;
  }
};
} // namespace cxITKMedianProjectionImageFilter

namespace nx::core
{
//------------------------------------------------------------------------------
std::string ITKMedianProjectionImageFilter::name() const
{
  return FilterTraits<ITKMedianProjectionImageFilter>::name;
}

//------------------------------------------------------------------------------
std::string ITKMedianProjectionImageFilter::className() const
{
  return FilterTraits<ITKMedianProjectionImageFilter>::className;
}

//------------------------------------------------------------------------------
Uuid ITKMedianProjectionImageFilter::uuid() const
{
  return FilterTraits<ITKMedianProjectionImageFilter>::uuid;
}

//------------------------------------------------------------------------------
std::string ITKMedianProjectionImageFilter::humanName() const
{
  return "ITK Median Projection Image Filter";
}

//------------------------------------------------------------------------------
std::vector<std::string> ITKMedianProjectionImageFilter::defaultTags() const
{
  return {className(), "ITKImageProcessing", "ITKMedianProjectionImageFilter", "ITKImageStatistics", "ImageStatistics"};
}

//------------------------------------------------------------------------------
Parameters ITKMedianProjectionImageFilter::parameters() const
{
  Parameters params;
  params.insertSeparator(Parameters::Separator{"Input Parameter(s)"});
  params.insert(std::make_unique<UInt32Parameter>(k_ProjectionDimension_Key, "Projection Dimension", "The dimension index to project. 0=Slowest moving dimension.", 0u));
  params.insertLinkableParameter(std::make_unique<BoolParameter>(k_RemoveOriginalGeometry_Key, "Perform In-Place", "Performs the projection in-place for the given Image Geometry", true));

  params.insertSeparator(Parameters::Separator{"Input Cell Data"});
  params.insert(std::make_unique<GeometrySelectionParameter>(k_InputImageGeomPath_Key, "Image Geometry", "Select the Image Geometry Group from the DataStructure.", DataPath({"Image Geometry"}),
                                                             GeometrySelectionParameter::AllowedTypes{IGeometry::Type::Image}));
  params.insert(std::make_unique<ArraySelectionParameter>(k_InputImageDataPath_Key, "Input Cell Data", "The image data that will be processed by this filter.", DataPath{},
                                                          nx::core::ITK::GetScalarPixelAllowedTypes()));

  params.insertSeparator(Parameters::Separator{"Output Data"});
  params.insert(std::make_unique<StringParameter>(k_OutputImageGeomName_Key, "Created Image Geometry", "The name of the projected geometry", "Projected Image"));
  params.insert(
      std::make_unique<DataObjectNameParameter>(k_OutputImageArrayName_Key, "Output Image Data Array", "The result of the processing will be stored in this Data Array.", "Output Image Data"));

  params.linkParameters(k_RemoveOriginalGeometry_Key, k_OutputImageGeomName_Key, false);

  return params;
}

//------------------------------------------------------------------------------
IFilter::VersionType ITKMedianProjectionImageFilter::parametersVersion() const
{
  return 1;
}

//------------------------------------------------------------------------------
IFilter::UniquePointer ITKMedianProjectionImageFilter::clone() const
{
  return std::make_unique<ITKMedianProjectionImageFilter>();
}

//------------------------------------------------------------------------------
IFilter::PreflightResult ITKMedianProjectionImageFilter::preflightImpl(const DataStructure& dataStructure, const Arguments& filterArgs, const MessageHandler& messageHandler,
                                                                       const std::atomic_bool& shouldCancel) const
{
  auto imageGeomPath = filterArgs.value<DataPath>(k_InputImageGeomPath_Key);
  auto selectedInputArray = filterArgs.value<DataPath>(k_InputImageDataPath_Key);
  auto outputArrayName = filterArgs.value<DataObjectNameParameter::ValueType>(k_OutputImageArrayName_Key);
  auto projectionDimension = filterArgs.value<uint32>(k_ProjectionDimension_Key);
  auto performInPlace = filterArgs.value<bool>(k_RemoveOriginalGeometry_Key);
  auto outputGeomName = filterArgs.value<std::string>(k_OutputImageGeomName_Key);

  return ProjectionUtilities::RunITKProjectionDataCheck<cxITKMedianProjectionImageFilter::ArrayOptionsType>(dataStructure, selectedInputArray, imageGeomPath, outputGeomName, performInPlace,
                                                                                                            outputArrayName);
}

//------------------------------------------------------------------------------
Result<> ITKMedianProjectionImageFilter::executeImpl(DataStructure& dataStructure, const Arguments& filterArgs, const PipelineFilter* pipelineNode, const MessageHandler& messageHandler,
                                                     const std::atomic_bool& shouldCancel) const
{
  auto imageGeomPath = filterArgs.value<DataPath>(k_InputImageGeomPath_Key);
  auto selectedInputArray = filterArgs.value<DataPath>(k_InputImageDataPath_Key);
  auto outputArrayName = filterArgs.value<DataObjectNameParameter::ValueType>(k_OutputImageArrayName_Key);
  auto outputImageGeomName = filterArgs.value<std::string>(k_OutputImageGeomName_Key);
  auto performInPlace = filterArgs.value<bool>(k_RemoveOriginalGeometry_Key);
  auto projectionDimension = filterArgs.value<uint32>(k_ProjectionDimension_Key);

  const cxITKMedianProjectionImageFilter::ITKMedianProjectionImageFilterFunctor itkFunctor = {projectionDimension};

  return ProjectionUtilities::RunITKProjectionExecute<cxITKMedianProjectionImageFilter::ArrayOptionsType>(dataStructure, selectedInputArray, imageGeomPath, shouldCancel, outputArrayName,
                                                                                                          performInPlace, itkFunctor, outputImageGeomName);
}
} // namespace nx::core
