# Package options
option(PHASTA_CHEF_ENABLED "Enable the chef package [ON|OFF]" OFF)
message(STATUS "PHASTA_CHEF_ENABLED: " ${PHASTA_CHEF_ENABLED})

# Package sources
set(SOURCES
  phCook.cc
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
  phBubble.cc
  phConstraint.cc
  phModelGeometry.cc
  phFilterMatching.cc
  phAxisymmetry.cc
  phInterfaceCutter.cc
)
if(NOT PHASTA_CHEF_ENABLED)
  set(SOURCES
    ${SOURCES}
    phstreamEmpty/phstreamEmpty.cc
  )
endif()

# Package headers
set(HEADERS
  chef.h
)

# Add the ph library
add_library(ph ${SOURCES})

# Handle the complexity of the chef enabled flag
if(PHASTA_CHEF_ENABLED)
  include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${PHASTACHEF_SOURCE_DIR}/phstream
  )
  target_link_libraries(ph ma mds sam phstream)
else()
  include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/phstreamEmpty
  )
  target_link_libraries(ph ma mds sam)
endif()

# Tell CMake to add this package to the SCOREC project
set_property(
  TARGET ph
  APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)
install(TARGETS ph EXPORT ph-targets DESTINATION lib)
install(FILES ${HEADERS} DESTINATION include)
install(EXPORT ph-targets DESTINATION lib/scorec)
