#include "FileSystemPathParameter.hpp"

#include "simplnx/Common/Any.hpp"
#include "simplnx/Common/StringLiteral.hpp"
#include "simplnx/Common/StringLiteralFormatting.hpp"
#include "simplnx/Utilities/FileUtilities.hpp"
#include "simplnx/Utilities/SIMPLConversion.hpp"
#include "simplnx/Utilities/StringUtilities.hpp"

#include <fmt/core.h>

#include <nlohmann/json.hpp>

#include <cctype>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace fs = std::filesystem;

using namespace nx::core;

namespace
{
//-----------------------------------------------------------------------------
Result<> ValidateInputFile(const FileSystemPathParameter::ValueType& path)
{
  try
  {
    if(!fs::exists(path))
    {
      return MakeErrorResult(-2, fmt::format("File System Path '{}' does not exist", path.string()));
    }

    if(!fs::is_regular_file(path))
    {
      return MakeErrorResult(-3, fmt::format("File System Path '{}' is not a file", path.string()));
    }
  } catch(const fs::filesystem_error& exception)
  {
    return MakeErrorResult(-9, fmt::format("Filesystem exception: {}", exception.what()));
  }
  return {};
}

//-----------------------------------------------------------------------------
Result<> ValidateInputDir(const FileSystemPathParameter::ValueType& path)
{

  try
  {
    if(!fs::exists(path))
    {
      return MakeErrorResult(-4, fmt::format("File System Path '{}' does not exist", path.string()));
    }
    if(!fs::is_directory(path))
    {
      return MakeErrorResult(-5, fmt::format("File System Path '{}' is not a file", path.string()));
    }
  } catch(const fs::filesystem_error& exception)
  {
    return MakeErrorResult(-10, fmt::format("Filesystem exception: {}", exception.what()));
  }

  return {};
}

//-----------------------------------------------------------------------------
Result<> ValidateOutputFile(const FileSystemPathParameter::ValueType& path)
{
  try
  {
    if(fs::exists(path) && fs::is_directory(path))
    {
      return MakeErrorResult(-8, fmt::format("File System Path '{}' exists AND is a directory. The Parameter is set to save a file.", path.string()));
    }

    auto result = FileUtilities::ValidateDirectoryWritePermission(path, true);
    if(result.invalid())
    {
      return result;
    }
    if(!fs::exists(path))
    {
      return MakeWarningVoidResult(-6, fmt::format("File System Path '{}' does not exist. It will be created during execution.", path.string()));
    }
  } catch(const fs::filesystem_error& exception)
  {
    return MakeErrorResult(-11, fmt::format("Filesystem exception: {}", exception.what()));
  }

  return {};
}

//-----------------------------------------------------------------------------
Result<> ValidateOutputDir(const FileSystemPathParameter::ValueType& path)
{
  try
  {
    auto result = FileUtilities::ValidateDirectoryWritePermission(path, false);
    if(result.invalid())
    {
      return result;
    }
    if(!fs::exists(path))
    {
      return MakeWarningVoidResult(-7, fmt::format("File System Path '{}' does not exist. It will be created during execution.", path.string()));
    }
  } catch(const fs::filesystem_error& exception)
  {
    return MakeErrorResult(-12, fmt::format("Filesystem exception: {}", exception.what()));
  }
  return {};
}
} // namespace

namespace nx::core
{
//-----------------------------------------------------------------------------
FileSystemPathParameter::FileSystemPathParameter(const std::string& name, const std::string& humanName, const std::string& helpText, const ValueType& defaultValue,
                                                 const ExtensionsType& extensionsType, PathType pathType, bool acceptAllExtensions)
: ValueParameter(name, humanName, helpText)
, m_DefaultValue(defaultValue)
, m_PathType(pathType)
, m_AvailableExtensions(extensionsType)
, m_acceptAllExtensions(acceptAllExtensions)
{
  ExtensionsType validatedExtensions;
  for(const auto& ext : m_AvailableExtensions)
  {
    if(ext.empty())
    {
      throw std::runtime_error("FileSystemPathParameter: One of the given extensions was empty. The filter is required to use non-emtpy extensions");
    }
    if(ext.at(0) != '.')
    {
      validatedExtensions.insert('.' + nx::core::StringUtilities::toLower(ext));
    }
    else
    {
      validatedExtensions.insert(nx::core::StringUtilities::toLower(ext));
    }
  }
  m_AvailableExtensions = validatedExtensions;
}

//-----------------------------------------------------------------------------
Uuid FileSystemPathParameter::uuid() const
{
  return ParameterTraits<FileSystemPathParameter>::uuid;
}

//-----------------------------------------------------------------------------
IParameter::AcceptedTypes FileSystemPathParameter::acceptedTypes() const
{
  return {typeid(ValueType)};
}

//------------------------------------------------------------------------------
IParameter::VersionType FileSystemPathParameter::getVersion() const
{
  return 1;
}

//-----------------------------------------------------------------------------
bool FileSystemPathParameter::acceptAllExtensions() const
{
  return m_acceptAllExtensions;
}

//-----------------------------------------------------------------------------
nlohmann::json FileSystemPathParameter::toJsonImpl(const std::any& value) const
{
  const auto& path = GetAnyRef<ValueType>(value);
  nlohmann::json json = path.string();
  return json;
}

//-----------------------------------------------------------------------------
Result<std::any> FileSystemPathParameter::fromJsonImpl(const nlohmann::json& json, VersionType version) const
{
  static constexpr StringLiteral prefix = "FilterParameter 'FileSystemPathParameter' JSON Error: ";

  if(!json.is_string())
  {
    return MakeErrorResult<std::any>(-2, fmt::format("{}JSON value for key '{}' is not a string", prefix, name()));
  }
  auto pathString = json.get<std::string>();
  std::filesystem::path path = pathString;
  return {{path}};
}

//-----------------------------------------------------------------------------
IParameter::UniquePointer FileSystemPathParameter::clone() const
{
  return std::make_unique<FileSystemPathParameter>(name(), humanName(), helpText(), m_DefaultValue, m_AvailableExtensions, m_PathType, m_acceptAllExtensions);
}

//-----------------------------------------------------------------------------
std::any FileSystemPathParameter::defaultValue() const
{
  return defaultPath();
}

//-----------------------------------------------------------------------------
typename FileSystemPathParameter::ValueType FileSystemPathParameter::defaultPath() const
{
  return m_DefaultValue;
}

//-----------------------------------------------------------------------------
FileSystemPathParameter::PathType FileSystemPathParameter::getPathType() const
{
  return m_PathType;
}

//-----------------------------------------------------------------------------
FileSystemPathParameter::ExtensionsType FileSystemPathParameter::getAvailableExtensions() const
{
  return m_AvailableExtensions;
}

//-----------------------------------------------------------------------------
Result<> FileSystemPathParameter::validate(const std::any& value) const
{
  const auto& path = GetAnyRef<ValueType>(value);
  return validatePath(path);
}

//-----------------------------------------------------------------------------
Result<> FileSystemPathParameter::validatePath(const ValueType& path) const
{
  try
  {
    const std::string prefix = fmt::format("Parameter Name: '{}'\n    Parameter Key: '{}'\n    Validation Error: ", humanName(), name());

    if(path.empty())
    {
      return nx::core::MakeErrorResult(-3001, fmt::format("{} File System Path must not be empty", prefix));
    }

    if(!m_acceptAllExtensions && (m_PathType == nx::core::FileSystemPathParameter::PathType::InputFile || m_PathType == nx::core::FileSystemPathParameter::PathType::OutputFile))
    {
      if(!path.has_extension())
      {
        return {nonstd::make_unexpected(std::vector<Error>{{-3002, fmt::format("{} File System Path must include a file extension.\n  FilePath: '{}'", prefix, path.string())}})};
      }
      std::string lowerExtension = nx::core::StringUtilities::toLower(path.extension().string());
      if(path.has_extension() && !m_AvailableExtensions.empty() && m_AvailableExtensions.find(lowerExtension) == m_AvailableExtensions.end())
      {
        return {nonstd::make_unexpected(std::vector<Error>{{-3003, fmt::format("{} File extension '{}' is not a valid file extension", prefix, path.extension().string())}})};
      }
    }

    switch(m_PathType)
    {
    case nx::core::FileSystemPathParameter::PathType::InputFile:
      return ValidateInputFile(path);
    case nx::core::FileSystemPathParameter::PathType::InputDir:
      return ValidateInputDir(path);
    case nx::core::FileSystemPathParameter::PathType::OutputFile:
      return ValidateOutputFile(path);
    case nx::core::FileSystemPathParameter::PathType::OutputDir:
      return ValidateOutputDir(path);
    }
  } catch(const fs::filesystem_error& exception)
  {
    return MakeErrorResult(-9, fmt::format("Filesystem exception: {}", exception.what()));
  }

  return {};
}

namespace SIMPLConversion
{
Result<InputFileFilterParameterConverter::ValueType> InputFileFilterParameterConverter::convert(const nlohmann::json& json)
{
  auto filePathReult = ReadInputFilePath(json, "InputFileFilterParameter");
  if(filePathReult.invalid())
  {
    return ConvertInvalidResult<ValueType>(std::move(filePathReult));
  }

  return {filePathReult.value()};
}

Result<OutputFileFilterParameterConverter::ValueType> OutputFileFilterParameterConverter::convert(const nlohmann::json& json)
{
  auto filePathReult = ReadInputFilePath(json, "OutputFileFilterParameter");
  if(filePathReult.invalid())
  {
    return ConvertInvalidResult<ValueType>(std::move(filePathReult));
  }

  return {filePathReult.value()};
}
} // namespace SIMPLConversion
} // namespace nx::core
