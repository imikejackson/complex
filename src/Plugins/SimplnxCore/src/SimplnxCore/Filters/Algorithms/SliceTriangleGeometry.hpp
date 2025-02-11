#pragma once

#include "SimplnxCore/SimplnxCore_export.hpp"

#include "simplnx/DataStructure/DataPath.hpp"
#include "simplnx/DataStructure/DataStructure.hpp"
#include "simplnx/DataStructure/Geometry/INodeGeometry2D.hpp"
#include "simplnx/Filter/IFilter.hpp"
#include "simplnx/Parameters/ChoicesParameter.hpp"
#include "simplnx/Parameters/DataGroupCreationParameter.hpp"
#include "simplnx/Parameters/DataGroupSelectionParameter.hpp"
#include "simplnx/Parameters/DataObjectNameParameter.hpp"

namespace nx::core
{

namespace slice_triangle_geometry::constants
{
constexpr ChoicesParameter::ValueType k_FullRange = 0;
constexpr ChoicesParameter::ValueType k_UserDefinedRange = 1;
} // namespace slice_triangle_geometry::constants

struct SIMPLNXCORE_EXPORT SliceTriangleGeometryInputValues
{
  ChoicesParameter::ValueType SliceRange;
  float32 Zstart;
  float32 Zend;
  float32 SliceResolution;
  bool HaveRegionIds;
  DataPath CADDataContainerName;
  DataPath RegionIdArrayPath;
  DataGroupCreationParameter::ValueType SliceDataContainerName;
  DataObjectNameParameter::ValueType EdgeAttributeMatrixName;
  DataObjectNameParameter::ValueType SliceIdArrayName;
  DataObjectNameParameter::ValueType SliceAttributeMatrixName;
};

/**
 * @class SliceTriangleGeometry
 * @brief This filter slices an input Triangle Geometry, producing an Edge Geometry
 */

class SIMPLNXCORE_EXPORT SliceTriangleGeometry
{
public:
  SliceTriangleGeometry(DataStructure& dataStructure, const IFilter::MessageHandler& mesgHandler, const std::atomic_bool& shouldCancel, SliceTriangleGeometryInputValues* inputValues);
  ~SliceTriangleGeometry() noexcept;

  SliceTriangleGeometry(const SliceTriangleGeometry&) = delete;
  SliceTriangleGeometry(SliceTriangleGeometry&&) noexcept = delete;
  SliceTriangleGeometry& operator=(const SliceTriangleGeometry&) = delete;
  SliceTriangleGeometry& operator=(SliceTriangleGeometry&&) noexcept = delete;

  Result<> operator()();

  const std::atomic_bool& getCancel();

protected:
  using TriStore = AbstractDataStore<INodeGeometry2D::SharedFaceList::value_type>;
  using VertsStore = AbstractDataStore<INodeGeometry0D::SharedVertexList::value_type>;

private:
  DataStructure& m_DataStructure;
  const SliceTriangleGeometryInputValues* m_InputValues = nullptr;
  const std::atomic_bool& m_ShouldCancel;
  const IFilter::MessageHandler& m_MessageHandler;
};

} // namespace nx::core
