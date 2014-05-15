
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

# Disable PCU threads
SET(ENABLE_THREADS OFF)

IF(TPL_ENABLE_PARASOLID)
  SET(Trilinos_ENABLE_SCORECpumi_geom_parasolid ON)
ENDIF()
IF(TPL_ENABLE_ACIS)
  SET(Trilinos_ENABLE_SCORECpumi_geom_acis ON)
ENDIF()

#
# C) Add the libraries, tests, and examples
#

TRIBITS_PROCESS_SUBPACKAGES()

TRIBITS_PACKAGE_DEF()


#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
