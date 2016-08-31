INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)

#
# A) Define the package
#

TRIBITS_PACKAGE(SCOREC)

#
# B) Set up package-specific options
#

# Tell SCOREC it is building inside Trilinos
SET(BUILD_IN_TRILINOS ON)

#
# C) Add the libraries, tests, and examples
#

TRIBITS_PROCESS_SUBPACKAGES()

TRIBITS_PACKAGE_DEF()

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
