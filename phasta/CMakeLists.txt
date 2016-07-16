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
  set(SOURCES ${SOURCES} phstreamEmpty/phstreamEmpty.cc)
endif()

# Package headers
set(HEADERS
  chef.h
)

# Add the ph library
add_library(ph ${SOURCES})

# Include directories
target_include_directories(ph INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    )
target_include_directories(ph PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    )
target_link_libraries(ph PRIVATE
    ma
    mds
    sam
    parma
    pcu
    apf_zoltan
    )

# Handle the complexity of the chef enabled flag
if(PHASTA_CHEF_ENABLED)
  target_link_libraries(ph PRIVATE phstream)
else()
  target_include_directories(ph PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/phstreamEmpty
    )
endif()

# Tell CMake to add this package to the SCOREC project
install(TARGETS ph EXPORT ph-targets DESTINATION lib)
install(FILES ${HEADERS} DESTINATION include)
install(EXPORT ph-targets DESTINATION lib/scorec)
