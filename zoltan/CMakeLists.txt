# Package options
option(ENABLE_ZOLTAN "Enable Zoltan interface [ON|OFF]" OFF)
message(STATUS "ENABLE_ZOLTAN: " ${ENABLE_ZOLTAN})

# Let CMake know where to find custom find package implementations
set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/")

# Check for and enable zoltan support
if(ENABLE_ZOLTAN)
  find_package(Zoltan REQUIRED)
endif()

# Package sources
if(ENABLE_ZOLTAN)
  set(SOURCES
    apfInterElement.cc
    apfZoltan.cc
    apfZoltanMesh.cc
    apfZoltanCallbacks.cc
  )
else()
  set(SOURCES
    apfInterElement.cc
    apfZoltanEmpty.cc
  )
endif()

# Package headers
set(HEADERS apfZoltan.h)

# Add the apf_zoltan library
add_library(apf_zoltan ${SOURCES})

# Include directories
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link this package to these libraries
target_link_libraries(apf_zoltan apf)

# Do extra work if zoltan is enabled
if(ENABLE_ZOLTAN)
  include_directories(
    ${ZOLTAN_INCLUDE_DIRS}
    ${PARMETIS_INCLUDE_DIRS}
  )
  target_link_libraries(
    apf_zoltan zoltan
    ${ZOLTAN_LIBRARIES}
    ${PARMETIS_LIBRARIES}
  )
endif()

# Tell CMake to add this package to the SCOREC project
set_property(
  TARGET apf_zoltan
  APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)
install(TARGETS apf_zoltan EXPORT apf_zoltan-targets DESTINATION lib)
install(FILES ${HEADERS} DESTINATION include)
install(EXPORT apf_zoltan-targets DESTINATION lib/scorec)
