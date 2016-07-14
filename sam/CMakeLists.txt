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
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link this package to these libraries
target_link_libraries(sam mds ma)

# Tell CMake to add this package to the SCOREC project
set_property(
  TARGET sam
  APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)
install(TARGETS sam EXPORT sam-targets DESTINATION lib)
install(FILES ${HEADERS} DESTINATION include)
install(EXPORT sam-targets DESTINATION lib/scorec)
