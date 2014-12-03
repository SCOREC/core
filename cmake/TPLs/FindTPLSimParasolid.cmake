INCLUDE(TribitsTplDeclareLibraries)

TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( SimParasolid
  REQUIRED_HEADERS SimParasolidKrnl.h
  REQUIRED_LIBS_NAMES SimParasolid260 pskernel
)
