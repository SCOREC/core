# - Try to find zoltan
# Once done this will define
#  TEUCHOS_FOUND - System has Teuchos
#  TEUCHOS_INCLUDE_DIRS - The Teuchos include directories
#  TEUCHOS_LIBRARIES - The libraries needed to use Teuchos
#  TEUCHOS_DEFINITIONS - Compiler switches required for using Teuchos


find_path(TEUCHOS_INCLUDE_DIR Teuchos_ParameterList.hpp)

find_library(TEUCHOS_LIBRARY teuchosparameterlist)

set(TEUCHOS_LIBRARIES ${TEUCHOS_LIBRARY})
set(TEUCHOS_INCLUDE_DIRS ${TEUCHOS_INCLUDE_DIR}) 

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set TEUCHOS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(
    TEUCHOS  
    DEFAULT_MSG
    TEUCHOS_LIBRARY TEUCHOS_INCLUDE_DIR
)

mark_as_advanced(TEUCHOS_INCLUDE_DIR TEUCHOS_LIBRARY )
