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

# Include directories
set(INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/noto
  ${CMAKE_CURRENT_SOURCE_DIR}/reel
)

# Add the pcu library
add_library(pcu ${SOURCES})
target_include_directories(pcu INTERFACE
    $<BUILD_INTERFACE:${INCLUDE_DIRS}>
    $<INSTALL_INTERFACE:include>
    )

# Check for and enable compression support
if(PCU_COMPRESS)
  find_package(BZip2 REQUIRED)
  target_include_directories(pcu PRIVATE ${BZIP2_INCLUDE_DIR})
  target_link_libraries(pcu PRIVATE ${BZIP2_LIBRARIES})
  target_add_definitions(pcu PRIVATE "-DPCU_BZIP")
endif()

install(TARGETS pcu EXPORT pcu-targets DESTINATION lib)
install(FILES ${HEADERS} DESTINATION include)
install(EXPORT pcu-targets DESTINATION lib/scorec)
