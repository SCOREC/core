if(DEFINED TRIBITS_PACKAGE)
  include(pkg_tribits.cmake)
  return()
endif()

if(NOT ENABLE_SIMMETRIX)
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
target_include_directories(apf_sim PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    )

option(ENABLE_FIELDSIM "Enable use of FieldSim from Simmetrix SimModSuite" FALSE)
message(STATUS "ENABLE_FIELDSIM: ${ENABLE_FIELDSIM}")
set(USE_FIELDSIM FALSE)
if( ${SIMMODSUITE_SimField_FOUND} AND ENABLE_FIELDSIM )
  set(USE_FIELDSIM TRUE)
endif()
set(USE_SIM_ADVMESHING ${HAVE_SIMADVMESHING})

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/apf_simConfig.h.in"
    "${CMAKE_CURRENT_BINARY_DIR}/apf_simConfig.h")
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/apf_simConfig.h
        DESTINATION include)
#directory containing apf_simConfig.h
target_include_directories(apf_sim PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)


scorec_export_library(apf_sim)

bob_end_subdir()
