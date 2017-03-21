tribits_package(SCORECpcu)

option(PCU_COMPRESS "Enable SMB compression using libbzip2 [ON|OFF]" OFF)

set(CMAKE_MODULE_PATH
   ${CMAKE_MODULE_PATH}
   "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/")

#Gets C99 support
find_package(C99 REQUIRED)
set(CMAKE_C_FLAGS "${C99_C_FLAGS} ${CMAKE_C_FLAGS}")

if (PCU_COMPRESS)
find_package(BZip2 REQUIRED)
if (NOT BZIP2_INCLUDE_DIR STREQUAL "")
include_directories(${BZIP2_INCLUDE_DIR})
set(BZ_INCLUDE "-I${BZIP2_INCLUDE_DIR}")
string(REGEX REPLACE "libbz2.*" " " BZ_LIB_DIR "${BZIP2_LIBRARIES}")
set(BZ_LINK "-L${BZ_LIB_DIR} -lbz2")
endif()
endif(PCU_COMPRESS)

set(PCU_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/noto
  ${CMAKE_CURRENT_SOURCE_DIR}/reel
)

#directory containing pcu header files
include_directories(${PCU_INCLUDE_DIRS})
#directory containing reel_config.h
include_directories("${PROJECT_BINARY_DIR}")

#Sources & Headers
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

set(HEADERS
   PCU.h
   pcu_io.h
   pcu_util.h
   noto/noto_malloc.h
   reel/reel.h)

tribits_add_library(
   pcu
   HEADERS ${HEADERS}
   SOURCES ${SOURCES})

if (PCU_COMPRESS)
  include_directories(${BZIP_INCLUDE_DIR})
  target_link_libraries(pcu ${BZIP2_LIBRARIES})
  add_definitions(-DPCU_BZIP)
endif (PCU_COMPRESS)

tribits_package_postprocess()
