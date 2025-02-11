#pragma once

#include "ITKImageProcessing/ITKImageProcessing_export.hpp"

#include "simplnx/Filter/FilterTraits.hpp"
#include "simplnx/Filter/IFilter.hpp"

namespace nx::core
{
/**
 * @class ITKImageReaderFilter
 * @brief This filter will ....
 */
class ITKIMAGEPROCESSING_EXPORT ITKImageReaderFilter : public IFilter
{
public:
  ITKImageReaderFilter() = default;
  ~ITKImageReaderFilter() noexcept override = default;

  ITKImageReaderFilter(const ITKImageReaderFilter&) = delete;
  ITKImageReaderFilter(ITKImageReaderFilter&&) noexcept = delete;

  ITKImageReaderFilter& operator=(const ITKImageReaderFilter&) = delete;
  ITKImageReaderFilter& operator=(ITKImageReaderFilter&&) noexcept = delete;

  // Parameter Keys
  static inline constexpr StringLiteral k_FileName_Key = "file_name";
  static inline constexpr StringLiteral k_ImageGeometryPath_Key = "output_geometry_path";
  static inline constexpr StringLiteral k_ImageDataArrayPath_Key = "image_data_array_name";
  static inline constexpr StringLiteral k_CellDataName_Key = "cell_attribute_matrix_name";

  static inline constexpr StringLiteral k_LengthUnit_Key = "length_unit_index";

  static inline constexpr StringLiteral k_ChangeOrigin_Key = "change_origin";
  static inline constexpr StringLiteral k_CenterOrigin_Key = "center_origin";
  static inline constexpr StringLiteral k_Origin_Key = "origin";

  static inline constexpr StringLiteral k_ChangeSpacing_Key = "change_spacing";
  static inline constexpr StringLiteral k_Spacing_Key = "spacing";

  static inline constexpr StringLiteral k_ChangeDataType_Key = "change_image_data_type";
  static inline constexpr StringLiteral k_ImageDataType_Key = "image_data_type_index";

  /**
   * @brief Reads SIMPL json and converts it simplnx Arguments.
   * @param json
   * @return Result<Arguments>
   */
  static Result<Arguments> FromSIMPLJson(const nlohmann::json& json);

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
   * @brief Returns parameters version integer.
   * Initial version should always be 1.
   * Should be incremented everytime the parameters change.
   * @return VersionType
   */
  VersionType parametersVersion() const override;

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
  Result<> executeImpl(DataStructure& dataStructure, const Arguments& filterArgs, const PipelineFilter* pipelineNode, const MessageHandler& messageHandler,
                       const std::atomic_bool& shouldCancel) const override;
};
} // namespace nx::core

SIMPLNX_DEF_FILTER_TRAITS(nx::core, ITKImageReaderFilter, "d72eaf98-9b1d-44c9-88f2-a5c3cf57b4f2");
