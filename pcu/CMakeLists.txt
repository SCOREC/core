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
  pcu_io.h
  pcu_util.h
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
  xsdk_add_tpl(BZIP2)
  find_package(BZip2 REQUIRED)
  target_include_directories(pcu PRIVATE ${BZIP2_INCLUDE_DIR})
  target_link_libraries(pcu PRIVATE ${BZIP2_LIBRARIES})
  target_compile_definitions(pcu PRIVATE "-DPCU_BZIP")
endif()


#######################
## Python Packaging  ##
#######################
if(ENABLE_PYTHON)
	set(MPICH_LNK "/usr/lib/x86_64-linux-gnu/libmpich.so")
	set(MPI_INC "/usr/include/mpich")
  find_package(SWIG REQUIRED)
  include(${SWIG_USE_FILE})
	set(CMAKE_FIND_LIBRARY_SUFFIXES ".so")
	include_directories(${PYTHON_INCLUDE_DIRS})
	include_directories(${MPI4PY_INCLUDE_DIR})
  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/noto)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/reel)
  set(CMAKE_SWIG_FLAGS "")

	set_source_files_properties(PCU.i PROPERTIES CPLUSPLUS ON)
	set_source_files_properties(PCU.i PROPERTIES SWIG_FLAGS "-includeall")
	swig_add_library(PCU LANGUAGE python SOURCES PCU.i pcu.c)
	swig_link_libraries(PCU ${PYTHON_LIBRARIES} pcu)
endif()

scorec_export_library(pcu)

bob_end_subdir()
