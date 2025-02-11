#pragma once

#include "OrientationAnalysis/OrientationAnalysis_export.hpp"

#include "simplnx/DataStructure/DataPath.hpp"
#include "simplnx/DataStructure/DataStructure.hpp"
#include "simplnx/Filter/IFilter.hpp"
#include "simplnx/Parameters/VectorParameter.hpp"

/**
* This is example code to put in the Execute Method of the filter.
  ComputeSchmidsInputValues inputValues;

  inputValues.LoadingDirection = filterArgs.value<VectorFloat32Parameter::ValueType>(k_LoadingDirection_Key);
  inputValues.StoreAngleComponents = filterArgs.value<bool>(k_StoreAngleComponents_Key);
  inputValues.OverrideSystem = filterArgs.value<bool>(k_OverrideSystem_Key);
  inputValues.SlipPlane = filterArgs.value<VectorFloat32Parameter::ValueType>(k_SlipPlane_Key);
  inputValues.SlipDirection = filterArgs.value<VectorFloat32Parameter::ValueType>(k_SlipDirection_Key);
  inputValues.FeaturePhasesArrayPath = filterArgs.value<DataPath>(k_FeaturePhasesArrayPath_Key);
  inputValues.AvgQuatsArrayPath = filterArgs.value<DataPath>(k_AvgQuatsArrayPath_Key);
  inputValues.CrystalStructuresArrayPath = filterArgs.value<DataPath>(k_CrystalStructuresArrayPath_Key);
  inputValues.SchmidsArrayName = filterArgs.value<DataPath>(k_SchmidsArrayName_Key);
  inputValues.SlipSystemsArrayName = filterArgs.value<DataPath>(k_SlipSystemsArrayName_Key);
  inputValues.PolesArrayName = filterArgs.value<DataPath>(k_PolesArrayName_Key);
  inputValues.PhisArrayName = filterArgs.value<DataPath>(k_PhisArrayName_Key);
  inputValues.LambdasArrayName = filterArgs.value<DataPath>(k_LambdasArrayName_Key);

  return ComputeSchmids(dataStructure, messageHandler, shouldCancel, &inputValues)();
*/

namespace nx::core
{

struct ORIENTATIONANALYSIS_EXPORT ComputeSchmidsInputValues
{
  VectorFloat32Parameter::ValueType LoadingDirection;
  bool StoreAngleComponents;
  bool OverrideSystem;
  VectorFloat32Parameter::ValueType SlipPlane;
  VectorFloat32Parameter::ValueType SlipDirection;
  DataPath FeaturePhasesArrayPath;
  DataPath AvgQuatsArrayPath;
  DataPath CrystalStructuresArrayPath;
  DataPath SchmidsArrayName;
  DataPath SlipSystemsArrayName;
  DataPath PolesArrayName;
  DataPath PhisArrayName;
  DataPath LambdasArrayName;
};

/**
 * @class
 */
class ORIENTATIONANALYSIS_EXPORT ComputeSchmids
{
public:
  ComputeSchmids(DataStructure& dataStructure, const IFilter::MessageHandler& mesgHandler, const std::atomic_bool& shouldCancel, ComputeSchmidsInputValues* inputValues);
  ~ComputeSchmids() noexcept;

  ComputeSchmids(const ComputeSchmids&) = delete;
  ComputeSchmids(ComputeSchmids&&) noexcept = delete;
  ComputeSchmids& operator=(const ComputeSchmids&) = delete;
  ComputeSchmids& operator=(ComputeSchmids&&) noexcept = delete;

  Result<> operator()();

  const std::atomic_bool& getCancel();

private:
  DataStructure& m_DataStructure;
  const ComputeSchmidsInputValues* m_InputValues = nullptr;
  const std::atomic_bool& m_ShouldCancel;
  const IFilter::MessageHandler& m_MessageHandler;
};

} // namespace nx::core
