#!/bin/bash -e
# \author mjackson
#
# This script assists with using clang-tidy.
#
# STEP 0;
# I recommend downloading or building the latest
# version of llvm with clang and the clang-tidy
# tools in a private repo. NOTE:  Apple does
# not currently distribute these.
#
# STEP 1:
# Build your source tree with cmake and use
# export PATH=${MYCLANG_FROM_STEP0}:$PATH
# export CXX=${MYCLANG_FROM_STEP0}
# export CC=${MYCLANG_FROM_STEP0}
# cd ${MYBLD}
# cmake -DCMAKE_EXPORT_COMPILE_COMMANDS:BOOL=ON ${MYSRC}
# to generate the compile_commands.json
# comple options database files.
#
# STEP 2:
# Run clang-tidy with the "run-clang-tidy.py" wrapper script from the llvm soruce tree
# Minor modifications to the "run-clang-tidy.py" script can allow
# for automatically skipping the "ThirdParty" directory or other performance
# enhancments.

SRCDIR=$1
BLDDIR=$2

this_script_name=$(basename $0)
CMTMSG=$(mktemp -q /tmp/${this_script_name}.XXXXXX)
FILES_TO_CHECK=$(mktemp -q /tmp/${this_script_name}_files.XXXXXX)

cat > ${CMTMSG} << EOF
STYLE: Adhere to naming conventions for constants

This check will replace constant variables to ensure they have an k_ prefix and CamelCase

SRCDIR=${SRCDIR} #My local SRC
BLDDIR=${BLDDIR} #My local BLD

cd ${BLDDIR}
run-clang-tidy -config=  -header-filter=.* -fix

EOF


PATH=/opt/local/clang+llvm-15.0.4-arm64-apple-darwin21.0/bin:$PATH
cd ${BLDDIR}
run-clang-tidy -header-filter=.* -fix -checks='-*,readability-identifier-naming' \
    -config="{CheckOptions:
[ { key: readability-identifier-naming.NamespaceCase,       value: lower_case },
{ key: readability-identifier-naming.ClassCase,           value: CamelCase  },
{ key: readability-identifier-naming.ClassMethodCase,     value: camelBack  },
{ key: readability-identifier-naming.PrivateMemberPrefix, value: m_         },
{ key: readability-identifier-naming.PrivateMember,       value: CamelCase  },
{ key: readability-identifier-naming.StructCase,          value: CamelCase  },
{ key: readability-identifier-naming.FunctionCase,        value: camelBack },
{ key: readability-identifier-naming.VariableCase,        value: camelBack },
{ key: readability-identifier-naming.GlobalConstantCase,  value: CamelCase },
{ key: readability-identifier-naming.GlobalConstantPrefix,value: s_ },
{ key: readability-identifier-naming.GlobalFunctionCase,  value: CamelCase },
{ key: readability-identifier-naming.LocalPointerCase,    value: camelBack },
{ key: readability-identifier-naming.LocalPointerSuffix,  value: Ptr },
{ key: readability-identifier-naming.TypeAliasCase,       value: CamelCase },
{ key: readability-identifier-naming.TypeAliasSuffix,     value: Type },
{ key: readability-identifier-naming.MacroDefinitionCase,  value: UPPER_CASE } ]
}"

cd ${SRCDIR}
#git add -A && git commit --file ${CMTMSG}
rm ${CMTMSG}
