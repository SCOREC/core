if(DEFINED TRIBITS_PACKAGE)
  include(pkg_tribits.cmake)
  return()
endif()

# Package options
option(PCU_COMPRESS "Enable SMB compression using libbzip2 [ON|OFF]" OFF)
message(STATUS "PCU_COMPRESS: " ${PCU_COMPRESS})

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
  pcu_util.c
  noto/noto_malloc.c
  reel/reel.c
)

# Package headers
set(HEADERS
  PCU.h
  pcu_util.h
  pcu_io.h
  reel/reel.h
)

# Add the pcu library
add_library(pcu ${SOURCES})

# Include directories
target_include_directories(pcu
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/reel>
    $<INSTALL_INTERFACE:include>
    PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/noto>
    )

# Check for and enable compression support
if(PCU_COMPRESS)
  find_package(BZip2 REQUIRED)
  target_include_directories(pcu PRIVATE ${BZIP2_INCLUDE_DIR})
  target_link_libraries(pcu PRIVATE ${BZIP2_LIBRARIES})
  target_compile_definitions(pcu PRIVATE "-DPCU_BZIP")
endif()

scorec_export_library(pcu)

bob_end_subdir()
