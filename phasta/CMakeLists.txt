project(ph)

find_package(ma PATHS ${CMAKE_BINARY_DIR})
find_package(mds PATHS ${CMAKE_BINARY_DIR})
find_package(gmi_sim PATHS ${CMAKE_BINARY_DIR})
find_package(sam PATHS ${CMAKE_BINARY_DIR})
set(PH_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${SAM_INCLUDE_DIRS}
  ${MA_INCLUDE_DIRS}
  ${MDS_INCLUDE_DIRS}
  ${GMI_SIM_INCLUDE_DIRS})
set(DEP_LIBS ${SAM_LIBS} ${MA_LIBS} ${MDS_LIBS} ${GMI_SIM_LIBS})
set(PH_LIBS ph ${DEP_LIBS})

set(HEADERS chef.h)

set(SOURCES
  phCook.cc
  phKitchen.cc
  ph.cc
  phAdjacent.cc
  phOutput.cc
  phLinks.cc
  phGeomBC.cc
  phBlock.cc
  phAdapt.cc
  phRestart.cc
  phInput.cc
  phPartition.cc
  phIO.c
  phBC.cc
  phAttrib.cc
  phBubble.cc
  phConstraint.cc
  phModelGeometry.cc
  phFilterMatching.cc
  phAxisymmetry.cc
  phInterfaceCutter.cc)

if(PHASTA_CHEF_ENABLED)
  set(DEP_LIBS ${DEP_LIBS} phstream)
  set(PH_INCLUDE_DIRS ${PH_INCLUDE_DIRS} ${PHASTACHEF_SOURCE_DIR}/phstream)
else()
  set(SOURCES ${SOURCES} phstreamEmpty/phstreamEmpty.cc)
  set(PH_INCLUDE_DIRS ${PH_INCLUDE_DIRS} phstreamEmpty)
endif()

include_directories(${PH_INCLUDE_DIRS})

add_library(ph ${SOURCES})
target_link_libraries(ph ${DEP_LIBS})

INSTALL(FILES ${HEADERS} DESTINATION include)
INSTALL(FILES "${PROJECT_BINARY_DIR}/libph.pc" DESTINATION lib/pkgconfig)
INSTALL(TARGETS 
  ph 
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/phConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/phConfig.cmake")

#pkgconfig {
set(prefix "${CMAKE_INSTALL_PREFIX}")
set(exec_prefix "${CMAKE_INSTALL_PREFIX}/bin")
set(libdir "${CMAKE_INSTALL_PREFIX}/lib")
set(includedir "${CMAKE_INSTALL_PREFIX}/include")
set(PACKAGE_VERSION "1.0")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/libph.pc.in"
    "${PROJECT_BINARY_DIR}/libph.pc"
    @ONLY)
#pkgconfig }

