project(sam)

find_package(ma PATHS ${CMAKE_BINARY_DIR})
find_package(mds PATHS ${CMAKE_BINARY_DIR})
set(SAM_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${MA_INCLUDE_DIRS}
  ${MDS_INCLUDE_DIRS})
set(DEP_LIBS ${MA_LIBS} ${MDS_LIBS})
set(SAM_LIBS sam ${DEP_LIBS})

set(HEADERS
  sam.h
  samSz.h
  samElementCount.h
  )

set(SOURCES
  sam.cc
  samSz.cc
  samElementCount.cc
  )

include_directories(${SAM_INCLUDE_DIRS})

add_library(sam ${SOURCES})
target_link_libraries(sam ${DEP_LIBS})

INSTALL(FILES ${HEADERS} DESTINATION include)
INSTALL(FILES "${PROJECT_BINARY_DIR}/libsam.pc" DESTINATION lib/pkgconfig)
INSTALL(TARGETS
  sam
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/samConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/samConfig.cmake")

#pkgconfig {
set(prefix "${CMAKE_INSTALL_PREFIX}")
set(exec_prefix "${CMAKE_INSTALL_PREFIX}/bin")
set(libdir "${CMAKE_INSTALL_PREFIX}/lib")
set(includedir "${CMAKE_INSTALL_PREFIX}/include")
set(PACKAGE_VERSION "1.0")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/libsam.pc.in"
    "${PROJECT_BINARY_DIR}/libsam.pc"
    @ONLY)
#pkgconfig }

