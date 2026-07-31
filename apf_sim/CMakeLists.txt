if(DEFINED TRIBITS_PACKAGE)
  include(pkg_tribits.cmake)
  return()
endif()

if(NOT ENABLE_SIMMETRIX)
  return()
endif()


#Sources & Headers
set(SOURCES apfSIM.cc)
set(HEADERS apfSIM.h apf_simConfig.h)

#Library
add_library(apf_sim ${SOURCES})
target_link_libraries(apf_sim PUBLIC apf gmi_sim)
target_link_libraries(apf_sim PUBLIC ${SIMMODSUITE_LIBS})
target_include_directories(gmi_sim PUBLIC ${SIMMODSUITE_INCLUDE_DIR})
target_include_directories(apf_sim PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    )

scorec_export_library(apf_sim)

bob_end_subdir()
