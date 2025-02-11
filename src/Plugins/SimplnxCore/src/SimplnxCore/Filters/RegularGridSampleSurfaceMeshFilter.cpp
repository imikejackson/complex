#include "RegularGridSampleSurfaceMeshFilter.hpp"

#include "SimplnxCore/Filters/Algorithms/RegularGridSampleSurfaceMesh.hpp"

#include "simplnx/DataStructure/DataPath.hpp"
#include "simplnx/Filter/Actions/CreateArrayAction.hpp"
#include "simplnx/Filter/Actions/CreateAttributeMatrixAction.hpp"
#include "simplnx/Filter/Actions/CreateGeometry1DAction.hpp"
#include "simplnx/Filter/Actions/CreateImageGeometryAction.hpp"
#include "simplnx/Filter/Actions/DeleteDataAction.hpp"
#include "simplnx/Parameters/ArraySelectionParameter.hpp"
#include "simplnx/Parameters/BoolParameter.hpp"
#include "simplnx/Parameters/ChoicesParameter.hpp"
#include "simplnx/Parameters/DataGroupCreationParameter.hpp"
#include "simplnx/Parameters/DataGroupSelectionParameter.hpp"
#include "simplnx/Parameters/DataObjectNameParameter.hpp"
#include "simplnx/Parameters/GeometrySelectionParameter.hpp"
#include "simplnx/Parameters/NumberParameter.hpp"
#include "simplnx/Parameters/VectorParameter.hpp"
#include "simplnx/Utilities/SIMPLConversion.hpp"

#include <sstream>

using namespace nx::core;

namespace nx::core
{
//------------------------------------------------------------------------------
std::string RegularGridSampleSurfaceMeshFilter::name() const
{
  return FilterTraits<RegularGridSampleSurfaceMeshFilter>::name.str();
}

//------------------------------------------------------------------------------
std::string RegularGridSampleSurfaceMeshFilter::className() const
{
  return FilterTraits<RegularGridSampleSurfaceMeshFilter>::className;
}

//------------------------------------------------------------------------------
Uuid RegularGridSampleSurfaceMeshFilter::uuid() const
{
  return FilterTraits<RegularGridSampleSurfaceMeshFilter>::uuid;
}

//------------------------------------------------------------------------------
std::string RegularGridSampleSurfaceMeshFilter::humanName() const
{
  return "Sample Triangle Geometry on Regular Grid";
}

//------------------------------------------------------------------------------
std::vector<std::string> RegularGridSampleSurfaceMeshFilter::defaultTags() const
{
  return {className(), "Sampling", "Spacing"};
}

//------------------------------------------------------------------------------
Parameters RegularGridSampleSurfaceMeshFilter::parameters() const
{
  Parameters params;

  // Create the parameter descriptors that are needed for this filter
  params.insertSeparator(Parameters::Separator{"Input Parameter(s)"});
  params.insert(std::make_unique<VectorUInt64Parameter>(k_Dimensions_Key, "Dimensions (Voxels)", "The dimensions of the created Image geometry", std::vector<uint64>{128, 128, 128},
                                                        std::vector<std::string>{"x", "y", "z"}));
  params.insert(
      std::make_unique<VectorFloat32Parameter>(k_Origin_Key, "Origin", "The origin of the created Image geometry", std::vector<float32>{0.0F, 0.0F, 0.0F}, std::vector<std::string>{"x", "y", "z"}));
  params.insert(
      std::make_unique<VectorFloat32Parameter>(k_Spacing_Key, "Spacing", "The spacing of the created Image geometry", std::vector<float32>{1.0F, 1.0F, 1.0F}, std::vector<std::string>{"x", "y", "z"}));
  params.insert(std::make_unique<ChoicesParameter>(k_LengthUnit_Key, "Length Units (For Description Only)", "The units to be displayed below", to_underlying(IGeometry::LengthUnit::Micrometer),
                                                   IGeometry::GetAllLengthUnitStrings()));

  params.insertSeparator(Parameters::Separator{"Input Data Objects"});
  params.insert(std::make_unique<GeometrySelectionParameter>(k_TriangleGeometryPath_Key, "Triangle Geometry", "The geometry to be sampled onto grid", DataPath{},
                                                             GeometrySelectionParameter::AllowedTypes{GeometrySelectionParameter::AllowedType ::Triangle}));
  params.insert(std::make_unique<ArraySelectionParameter>(k_SurfaceMeshFaceLabelsArrayPath_Key, "Face Labels", "Array specifying which Features are on either side of each Face", DataPath{},
                                                          ArraySelectionParameter::AllowedTypes{nx::core::DataType::int32}, ArraySelectionParameter::AllowedComponentShapes{{2}}));

  params.insertSeparator(Parameters::Separator{"Output Image Geometry"});
  params.insert(std::make_unique<DataGroupCreationParameter>(k_ImageGeomPath_Key, "Image Geometry", "The name and path for the image geometry to be created", DataPath{}));
  params.insertSeparator(Parameters::Separator{"Output Cell Attribute Matrix"});
  params.insert(std::make_unique<DataObjectNameParameter>(k_CellAMName_Key, "Cell Attribute Matrix", "The name for the cell data Attribute Matrix within the Image geometry", "Cell Data"));
  params.insertSeparator(Parameters::Separator{"Output Cell Data"});
  params.insert(std::make_unique<DataObjectNameParameter>(k_FeatureIdsArrayName_Key, "Feature Ids", "The name for the feature ids array in cell data Attribute Matrix", "Feature Ids"));

  return params;
}

//------------------------------------------------------------------------------
IFilter::VersionType RegularGridSampleSurfaceMeshFilter::parametersVersion() const
{
  return 1;
}

//------------------------------------------------------------------------------
IFilter::UniquePointer RegularGridSampleSurfaceMeshFilter::clone() const
{
  return std::make_unique<RegularGridSampleSurfaceMeshFilter>();
}

//------------------------------------------------------------------------------
IFilter::PreflightResult RegularGridSampleSurfaceMeshFilter::preflightImpl(const DataStructure& dataStructure, const Arguments& filterArgs, const MessageHandler& messageHandler,
                                                                           const std::atomic_bool& shouldCancel) const
{
  auto pSurfaceMeshFaceLabelsArrayPathValue = filterArgs.value<DataPath>(k_SurfaceMeshFaceLabelsArrayPath_Key);
  auto pDimensionsValue = filterArgs.value<VectorUInt64Parameter::ValueType>(k_Dimensions_Key);
  auto pSpacingValue = filterArgs.value<VectorFloat32Parameter::ValueType>(k_Spacing_Key);
  auto pOriginValue = filterArgs.value<VectorFloat32Parameter::ValueType>(k_Origin_Key);
  auto pLengthUnitValue = filterArgs.value<ChoicesParameter::ValueType>(k_LengthUnit_Key);
  auto pImageGeomPathValue = filterArgs.value<DataPath>(k_ImageGeomPath_Key);
  auto pCellAMNameValue = filterArgs.value<std::string>(k_CellAMName_Key);
  auto pFeatureIdsArrayNameValue = filterArgs.value<std::string>(k_FeatureIdsArrayName_Key);
  auto triangleGeometryPath = filterArgs.value<DataPath>(k_TriangleGeometryPath_Key);

  nx::core::Result<OutputActions> resultOutputActions;

  std::vector<usize> tupleDims = {static_cast<usize>(pDimensionsValue[0]), static_cast<usize>(pDimensionsValue[1]), static_cast<usize>(pDimensionsValue[2])};

  {
    auto createDataGroupAction = std::make_unique<CreateImageGeometryAction>(pImageGeomPathValue, tupleDims, std::vector<float32>(pOriginValue), std::vector<float32>(pSpacingValue), pCellAMNameValue);
    resultOutputActions.value().appendAction(std::move(createDataGroupAction));
  }
  DataPath featIdsPath = pImageGeomPathValue.createChildPath(pCellAMNameValue).createChildPath(pFeatureIdsArrayNameValue);
  {
    auto createDataGroupAction = std::make_unique<CreateArrayAction>(DataType::int32, tupleDims, std::vector<usize>{1}, featIdsPath);
    resultOutputActions.value().appendAction(std::move(createDataGroupAction));
  }

  std::vector<PreflightValue> preflightUpdatedValues;
  std::stringstream boxDimensions = std::stringstream();

  std::string lengthUnit = IGeometry::LengthUnitToString(static_cast<IGeometry::LengthUnit>(pLengthUnitValue));

  boxDimensions << "X Range: " << std::setprecision(8) << std::noshowpoint << pOriginValue[0] << " to " << std::setprecision(8) << std::noshowpoint
                << (pOriginValue[0] + (pDimensionsValue[0] * pSpacingValue[0])) << " (Delta: " << std::setprecision(8) << std::noshowpoint << (pDimensionsValue[0] * pSpacingValue[0]) << ") "
                << lengthUnit << "\n";
  boxDimensions << "Y Range: " << std::setprecision(8) << std::noshowpoint << pOriginValue[1] << " to " << std::setprecision(8) << std::noshowpoint
                << (pOriginValue[1] + (pDimensionsValue[1] * pSpacingValue[1])) << " (Delta: " << std::setprecision(8) << std::noshowpoint << (pDimensionsValue[1] * pSpacingValue[1]) << ") "
                << lengthUnit << "\n";
  boxDimensions << "Z Range: " << std::setprecision(8) << std::noshowpoint << pOriginValue[2] << " to " << std::setprecision(8) << std::noshowpoint
                << (pOriginValue[2] + (pDimensionsValue[2] * pSpacingValue[2])) << " (Delta: " << std::setprecision(8) << std::noshowpoint << (pDimensionsValue[2] * pSpacingValue[2]) << ") "
                << lengthUnit << "\n";

  float32 vol = (pDimensionsValue[0] * pSpacingValue[0]) * (pDimensionsValue[1] * pSpacingValue[1]) * (pDimensionsValue[2] * pSpacingValue[2]);

  boxDimensions << "Volume: " << std::setprecision(8) << std::noshowpoint << vol << " " << lengthUnit << "s ^3"
                << "\n";

  preflightUpdatedValues.push_back({"BoxDimensions", boxDimensions.str()});

  /////////////////////////////////////////////////////////////////////////////
  // CREATE THE EDGE GEOMETRY THAT WILL BE USED FOR THE POLYGONS
  // This Geometry will be deleted after the filter is completed
  {
    DataPath pSliceDataContainerNameValue({fmt::format(".{}_sliced", triangleGeometryPath.getTargetName())});
    std::string pEdgeAttributeMatrixNameValue("EdgeAttributeMatrix");
    std::string pSliceIdArrayNameValue("SliceIds");
    DataPath pRegionIdArrayPathValue({"NOT USED"});
    std::string pSliceAttributeMatrixNameValue("SliceAttributeMatrix");
    // create the edge geometry
    {
      auto createGeometryAction = std::make_unique<CreateEdgeGeometryAction>(pSliceDataContainerNameValue, 1, 2, INodeGeometry0D::k_VertexAttributeMatrixName, pEdgeAttributeMatrixNameValue,
                                                                             EdgeGeom::k_SharedVertexListName, EdgeGeom::k_SharedEdgeListName);
      resultOutputActions.value().appendAction(std::move(createGeometryAction));
    }

    std::vector<size_t> tDims = {1};
    const std::vector<size_t> compDims = {1};
    {
      DataPath path = pSliceDataContainerNameValue.createChildPath(pEdgeAttributeMatrixNameValue).createChildPath(pSliceIdArrayNameValue);
      auto createArray = std::make_unique<CreateArrayAction>(DataType::int32, tDims, compDims, path);
      resultOutputActions.value().appendAction(std::move(createArray));
    }

    DataPath featureSliceAttrMatPath = pSliceDataContainerNameValue.createChildPath(pSliceAttributeMatrixNameValue);
    {
      auto createAttributeMatrixAction = std::make_unique<CreateAttributeMatrixAction>(featureSliceAttrMatPath, tDims);
      resultOutputActions.value().appendAction(std::move(createAttributeMatrixAction));
    }

    auto deferredDeleteGeometryAction = std::make_unique<DeleteDataAction>(pSliceDataContainerNameValue);
    resultOutputActions.value().appendDeferredAction(std::move(deferredDeleteGeometryAction));
  }
  /////////////////////////////////////////////////////////////////////////////

  // Return both the resultOutputActions and the preflightUpdatedValues via std::move()
  return {std::move(resultOutputActions), std::move(preflightUpdatedValues)};
}

//------------------------------------------------------------------------------
Result<> RegularGridSampleSurfaceMeshFilter::executeImpl(DataStructure& dataStructure, const Arguments& filterArgs, const PipelineFilter* pipelineNode, const MessageHandler& messageHandler,
                                                         const std::atomic_bool& shouldCancel) const
{
  RegularGridSampleSurfaceMeshInputValues inputValues;

  inputValues.Dimensions = filterArgs.value<VectorUInt64Parameter::ValueType>(k_Dimensions_Key);
  inputValues.Spacing = filterArgs.value<VectorFloat32Parameter::ValueType>(k_Spacing_Key);
  inputValues.Origin = filterArgs.value<VectorFloat32Parameter::ValueType>(k_Origin_Key);
  inputValues.TriangleGeometryPath = filterArgs.value<DataPath>(k_TriangleGeometryPath_Key);
  inputValues.SurfaceMeshFaceLabelsArrayPath = filterArgs.value<DataPath>(k_SurfaceMeshFaceLabelsArrayPath_Key);
  inputValues.FeatureIdsArrayPath =
      filterArgs.value<DataPath>(k_ImageGeomPath_Key).createChildPath(filterArgs.value<std::string>(k_CellAMName_Key)).createChildPath(filterArgs.value<std::string>(k_FeatureIdsArrayName_Key));
  inputValues.ImageGeometryOutputPath = filterArgs.value<DataPath>(k_ImageGeomPath_Key);

  return RegularGridSampleSurfaceMesh(dataStructure, messageHandler, shouldCancel, &inputValues)();
}

namespace
{
namespace SIMPL
{
constexpr StringLiteral k_SurfaceMeshFaceLabelsArrayPathKey = "SurfaceMeshFaceLabelsArrayPath";
constexpr StringLiteral k_DimensionsKey = "Dimensions";
constexpr StringLiteral k_SpacingKey = "Spacing";
constexpr StringLiteral k_OriginKey = "Origin";
constexpr StringLiteral k_LengthUnitKey = "LengthUnit";
constexpr StringLiteral k_DataContainerNameKey = "DataContainerName";
constexpr StringLiteral k_CellAttributeMatrixNameKey = "CellAttributeMatrixName";
constexpr StringLiteral k_FeatureIdsArrayNameKey = "FeatureIdsArrayName";
} // namespace SIMPL
} // namespace

Result<Arguments> RegularGridSampleSurfaceMeshFilter::FromSIMPLJson(const nlohmann::json& json)
{
  Arguments args = RegularGridSampleSurfaceMeshFilter().getDefaultArguments();

  std::vector<Result<>> results;

  results.push_back(
      SIMPLConversion::ConvertParameter<SIMPLConversion::DataContainerSelectionFilterParameterConverter>(args, json, SIMPL::k_SurfaceMeshFaceLabelsArrayPathKey, k_TriangleGeometryPath_Key));
  results.push_back(
      SIMPLConversion::ConvertParameter<SIMPLConversion::DataArraySelectionFilterParameterConverter>(args, json, SIMPL::k_SurfaceMeshFaceLabelsArrayPathKey, k_SurfaceMeshFaceLabelsArrayPath_Key));
  results.push_back(SIMPLConversion::ConvertParameter<SIMPLConversion::UInt64Vec3FilterParameterConverter>(args, json, SIMPL::k_DimensionsKey, k_Dimensions_Key));
  results.push_back(SIMPLConversion::ConvertParameter<SIMPLConversion::FloatVec3FilterParameterConverter>(args, json, SIMPL::k_SpacingKey, k_Spacing_Key));
  results.push_back(SIMPLConversion::ConvertParameter<SIMPLConversion::FloatVec3FilterParameterConverter>(args, json, SIMPL::k_OriginKey, k_Origin_Key));
  results.push_back(SIMPLConversion::ConvertParameter<SIMPLConversion::ChoiceFilterParameterConverter>(args, json, SIMPL::k_LengthUnitKey, k_LengthUnit_Key));
  results.push_back(SIMPLConversion::ConvertParameter<SIMPLConversion::DataContainerCreationFilterParameterConverter>(args, json, SIMPL::k_DataContainerNameKey, k_ImageGeomPath_Key));
  results.push_back(SIMPLConversion::ConvertParameter<SIMPLConversion::LinkedPathCreationFilterParameterConverter>(args, json, SIMPL::k_CellAttributeMatrixNameKey, k_CellAMName_Key));
  results.push_back(SIMPLConversion::ConvertParameter<SIMPLConversion::LinkedPathCreationFilterParameterConverter>(args, json, SIMPL::k_FeatureIdsArrayNameKey, k_FeatureIdsArrayName_Key));

  Result<> conversionResult = MergeResults(std::move(results));

  return ConvertResultTo<Arguments>(std::move(conversionResult), std::move(args));
}
} // namespace nx::core
