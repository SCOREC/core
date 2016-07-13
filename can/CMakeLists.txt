# Package headers
set(HEADERS
  canArray.h
  canNewArray.h
)

# Add the can 'library' as a CMake Interface
add_library(can INTERFACE)

# Tell CMake to add this package nicely to the SCOREC project
set_property(
  TARGET can
  APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)
install(TARGETS can EXPORT can-targets DESTINATION lib)
install(FILES ${HEADERS} DESTINATION include)
install(EXPORT can-targets DESTINATION lib/scorec)
