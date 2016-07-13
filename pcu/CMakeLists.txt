# Package options
option(PCU_COMPRESS "Enable SMB compression using libbzip2 [ON|OFF]" OFF)
message(STATUS "PCU_COMPRESS: " ${PCU_COMPRESS})

# Let CMake know where to find custom find package implementations
set(CMAKE_MODULE_PATH
  ${CMAKE_MODULE_PATH}
  "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/")

# Gets C99 support
find_package(C99 REQUIRED)
set(CMAKE_C_FLAGS "${C99_C_FLAGS} ${CMAKE_C_FLAGS}")
message(STATUS "CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")

# Check for and enable compression support
if(PCU_COMPRESS)
  find_package(BZip2 REQUIRED)
endif(PCU_COMPRESS)

# Package sources
set(SOURCES
  pcu.c
  pcu_aa.c
  pcu_coll.c
  pcu_io.c
  pcu_buffer.c
  pcu_mpi.c
  pcu_msg.c
  pcu_order.c
  pcu_pmpi.c
  noto/noto_malloc.c
  reel/reel.c
)

# Package headers
set(HEADERS
  PCU.h
  pcu_io.h
  noto/noto_malloc.h
  reel/reel.h
)

# Add the pcu library
add_library(pcu ${SOURCES})

# Include directories
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/noto
  ${CMAKE_CURRENT_SOURCE_DIR}/reel
)

# Do extra work if compression is enabled
if(PCU_COMPRESS)
 include_directories(${BZIP2_INCLUDE_DIR})
 target_link_libraries(pcu ${BZIP2_LIBRARIES})
 add_definitions(-DPCU_BZIP)
endif()

# Tell CMake to add this package to the SCOREC project
set_property(
  TARGET pcu
  APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include>
)
install(TARGETS pcu EXPORT pcu-targets DESTINATION lib)
install(FILES ${HEADERS} DESTINATION include)
install(EXPORT pcu-targets DESTINATION lib/scorec)
