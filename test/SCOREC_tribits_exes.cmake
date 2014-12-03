include(TribitsSubPackageMacros) 
tribits_subpackage(exe)

tribits_add_executable(
  convert
  SOURCES convert.cc
  INSTALLABLE
)

tribits_subpackage_postprocess()
