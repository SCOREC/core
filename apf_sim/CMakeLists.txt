if(DEFINED TRIBITS_PACKAGE)
  include(pkg_tribits.cmake)
  return()
endif()

#Sources & Headers
set(SOURCES apfSIM.cc)
set(HEADERS apfSIM.h)

#Library
add_library(apf_sim ${SOURCES})
target_link_libraries(apf_sim PUBLIC apf gmi_sim)
target_link_libraries(apf_sim PUBLIC ${SIMMODSUITE_LIBS})
target_include_directories(gmi_sim PUBLIC ${SIMMODSUITE_INCLUDE_DIR})

if(SIMMODSUITE_SimField_FOUND)
  target_compile_definitions(apf_sim PRIVATE -DUSE_FIELDSIM)
endif()

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/apf_simConfig.h.in"
    "${CMAKE_CURRENT_BINARY_DIR}/apf_simConfig.h")
#directory containing apf_simConfig.h
target_include_directories(apf_sim "${CMAKE_CURRENT_BINARY_DIR}")

scorec_export_library(apf_sim)

bob_end_subdir()
