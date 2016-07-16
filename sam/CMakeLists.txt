# Package sources
set(SOURCES
  samSz.cc
  sam.cc
)

# Package headers
set(HEADERS
  sam.h
)

# Add the sam library
add_library(sam ${SOURCES})

# Include directories
target_include_directories(sam INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    )

# Link this package to these libraries
target_link_libraries(sam PUBLIC apf)

# Tell CMake to add this package to the SCOREC project
install(TARGETS sam EXPORT sam-targets DESTINATION lib)
install(FILES ${HEADERS} DESTINATION include)
install(EXPORT sam-targets DESTINATION lib/scorec)
