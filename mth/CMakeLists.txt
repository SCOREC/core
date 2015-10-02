cmake_minimum_required(VERSION 2.8.6)
set(PACKAGE_VERSION "1.0")

project(mth)

set(MTH_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR})

find_package(can PATHS ${CMAKE_BINARY_DIR})
set(MTH_INCLUDE_DIRS
  ${MTH_INCLUDE_DIRS}
  ${CAN_INCLUDE_DIRS}
)
set(DEP_LIBS ${CAN_LIBS})
set(MTH_LIBS mth ${DEP_LIBS})

set(MTH_LIBS mth)

set(MTH_SOURCES
  mthQR.cc
)

set(MTH_HEADERS
  mth.h
  mth_def.h
  mthVector.h
  mthMatrix.h
  mthTensor.h
  mthQR.h
  mthAD.h
)

# THIS IS WHERE TRIBITS GETS HEADERS
include_directories(${MTH_INCLUDE_DIRS})

#Library
add_library(mth ${MTH_SOURCES})

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/mthConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/mthConfig.cmake")

#pkgconfig {
set(prefix "${CMAKE_INSTALL_PREFIX}")
set(exec_prefix "${CMAKE_INSTALL_PREFIX}/bin")
set(libdir "${CMAKE_INSTALL_PREFIX}/lib")
set(includedir "${CMAKE_INSTALL_PREFIX}/include")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/libmth.pc.in"
    "${PROJECT_BINARY_DIR}/libmth.pc"
    @ONLY)
#pkgconfig }

INSTALL(FILES ${MTH_HEADERS} DESTINATION include)
INSTALL(FILES "${PROJECT_BINARY_DIR}/libmth.pc" DESTINATION lib/pkgconfig)
