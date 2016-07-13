# Let CMake know where to find package implementations
set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/")

# Gets C99 support
find_package(C99 REQUIRED)
set(CMAKE_C_FLAGS "${C99_C_FLAGS} ${CMAKE_C_FLAGS}")
message(STATUS "CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")

# Include directories
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

# Package sources
set(SOURCES
  gmi.c
  agm.c
  gmi_base.c
  gmi_file.c
  gmi_lookup.c
  gmi_mesh.c
  gmi_null.c
  gmi_analytic.c
)

# Package headers
set(HEADERS
  gmi.h
  agm.h
  gmi_base.h
  gmi_lookup.h
  gmi_mesh.h
  gmi_null.h
  gmi_analytic.h
)

# Add the gmi library
add_library(gmi ${SOURCES})

# Tell Cmake to add this package to the SCOREC project
set_property(
  TARGET gmi
  APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)
install(TARGETS gmi EXPORT gmi-targets DESTINATION lib)
install(FILES ${HEADERS} DESTINATION include)
install(EXPORT gmi-targets DESTINATION lib/scorec)
