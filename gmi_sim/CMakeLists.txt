if(DEFINED TRIBITS_PACKAGE)
  include(pkg_tribits.cmake)
  return()
endif()

if(NOT ENABLE_SIMMETRIX)
  return()
endif()

option(SIM_PARASOLID "Enable Parasolid interface" OFF)
option(SIM_ACIS "Enable Acis interface" OFF)
option(SIM_DISCRETE "Enable Simmetrix discrete modeler interface" ON)
# this file brings SIM_PARASOLID, SIM_ACIS, SIM_DISCRETE from CMake to C++
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/gmi_sim_config.h.in"
               "${CMAKE_CURRENT_BINARY_DIR}/gmi_sim_config.h")

#Sources & Headers
set(SOURCES gmi_sim.cc)
set(HEADERS gmi_sim.h)

add_library(gmi_sim ${SOURCES})

target_include_directories(gmi_sim PUBLIC ${SIMMODSUITE_INCLUDE_DIR})
target_link_libraries(gmi_sim PUBLIC gmi pcu ${SIMMODSUITE_LIBS})

# Include directories
target_include_directories(gmi_sim PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    )
# make sure the compiler can find the config header
target_include_directories(gmi_sim PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)

scorec_export_library(gmi_sim)

bob_end_subdir()
