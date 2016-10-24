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
  phstream.cc
)

if(ENABLE_SIMMETRIX)
  set(SOURCES ${SOURCES}
      phAttrib.cc
     )
endif()

# Package headers
set(HEADERS
  chef.h
  phInput.h
  phstream.h
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
target_link_libraries(ph
    PUBLIC
      ma
      mds
      sam
      parma
      pcu
      apf_zoltan
    )
if(ENABLE_SIMMETRIX)
  target_link_libraries(ph PUBLIC
      gmi_sim
      apf_sim
    )
endif()

scorec_export_library(ph)

bob_end_subdir()
