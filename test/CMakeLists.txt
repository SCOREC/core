find_package(mds REQUIRED PATHS ${CMAKE_BINARY_DIR})

include_directories(${MDS_INCLUDE_DIRS})

add_executable(gmi_test EXCLUDE_FROM_ALL gmi_test.c)
add_executable(mds_test EXCLUDE_FROM_ALL mds_test.cc)

set(tests
  gmi_test
  mds_test)

foreach(test ${tests})
  target_link_libraries(${test} ${MDS_LIBS})
endforeach()
