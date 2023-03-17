#pragma once

#include "OrientationAnalysis/OrientationAnalysis_export.hpp"

#include "complex/Filter/FilterTraits.hpp"
#include "complex/Filter/IFilter.hpp"

namespace complex
{
/**
 * @class FindTwinBoundarySchmidFactors
 * @brief This filter will ....
 */
class ORIENTATIONANALYSIS_EXPORT FindTwinBoundarySchmidFactors : public IFilter
{
public:
  FindTwinBoundarySchmidFactors() = default;
  ~FindTwinBoundarySchmidFactors() noexcept override = default;

  FindTwinBoundarySchmidFactors(const FindTwinBoundarySchmidFactors&) = delete;
  FindTwinBoundarySchmidFactors(FindTwinBoundarySchmidFactors&&) noexcept = delete;

  FindTwinBoundarySchmidFactors& operator=(const FindTwinBoundarySchmidFactors&) = delete;
  FindTwinBoundarySchmidFactors& operator=(FindTwinBoundarySchmidFactors&&) noexcept = delete;

  // Parameter Keys
  static inline constexpr StringLiteral k_LoadingDir_Key = "loading_dir";
  static inline constexpr StringLiteral k_WriteFile_Key = "write_file";
  static inline constexpr StringLiteral k_TwinBoundarySchmidFactorsFile_Key = "twin_boundary_schmid_factors_file";
  static inline constexpr StringLiteral k_AvgQuatsArrayPath_Key = "avg_quats_array_path";
  static inline constexpr StringLiteral k_FeaturePhasesArrayPath_Key = "feature_phases_array_path";
  static inline constexpr StringLiteral k_CrystalStructuresArrayPath_Key = "crystal_structures_array_path";
  static inline constexpr StringLiteral k_SurfaceMeshFaceLabelsArrayPath_Key = "surface_mesh_face_labels_array_path";
  static inline constexpr StringLiteral k_SurfaceMeshFaceNormalsArrayPath_Key = "surface_mesh_face_normals_array_path";
  static inline constexpr StringLiteral k_SurfaceMeshTwinBoundaryArrayPath_Key = "surface_mesh_twin_boundary_array_path";
  static inline constexpr StringLiteral k_SurfaceMeshTwinBoundarySchmidFactorsArrayName_Key = "surface_mesh_twin_boundary_schmid_factors_array_name";

  /**
   * @brief Returns the name of the filter.
   * @return
   */
  std::string name() const override;

  /**
   * @brief Returns the C++ classname of this filter.
   * @return
   */
  std::string className() const override;

  /**
   * @brief Returns the uuid of the filter.
   * @return
   */
  Uuid uuid() const override;

  /**
   * @brief Returns the human readable name of the filter.
   * @return
   */
  std::string humanName() const override;

  /**
   * @brief Returns the default tags for this filter.
   * @return
   */
  std::vector<std::string> defaultTags() const override;

  /**
   * @brief Returns the parameters of the filter (i.e. its inputs)
   * @return
   */
  Parameters parameters() const override;

  /**
   * @brief Returns a copy of the filter.
   * @return
   */
  UniquePointer clone() const override;

protected:
  /**
   * @brief Takes in a DataStructure and checks that the filter can be run on it with the given arguments.
   * Returns any warnings/errors. Also returns the changes that would be applied to the DataStructure.
   * Some parts of the actions may not be completely filled out if all the required information is not available at preflight time.
   * @param dataStructure The input DataStructure instance
   * @param filterArgs These are the input values for each parameter that is required for the filter
   * @param messageHandler The MessageHandler object
   * @return Returns a Result object with error or warning values if any of those occurred during execution of this function
   */
  PreflightResult preflightImpl(const DataStructure& dataStructure, const Arguments& filterArgs, const MessageHandler& messageHandler, const std::atomic_bool& shouldCancel) const override;

  /**
   * @brief Applies the filter's algorithm to the DataStructure with the given arguments. Returns any warnings/errors.
   * On failure, there is no guarantee that the DataStructure is in a correct state.
   * @param dataStructure The input DataStructure instance
   * @param filterArgs These are the input values for each parameter that is required for the filter
   * @param messageHandler The MessageHandler object
   * @return Returns a Result object with error or warning values if any of those occurred during execution of this function
   */
  Result<> executeImpl(DataStructure& data, const Arguments& filterArgs, const PipelineFilter* pipelineNode, const MessageHandler& messageHandler, const std::atomic_bool& shouldCancel) const override;
};
} // namespace complex

COMPLEX_DEF_FILTER_TRAITS(complex, FindTwinBoundarySchmidFactors, "19ad8363-40b5-45aa-829a-4406f35b00ce");
