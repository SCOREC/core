find_package(mds REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(ma REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(gmi_sim REQUIRED PATHS ${CMAKE_BINARY_DIR})

include_directories(
  ${MA_INCLUDE_DIRS}
  ${MDS_INCLUDE_DIRS}
  ${GMI_SIM_INCLUDE_DIRS})

add_executable(gmi_test EXCLUDE_FROM_ALL gmi_test.cc)
add_executable(mds_test EXCLUDE_FROM_ALL mds_test.cc)
add_executable(ma_test EXCLUDE_FROM_ALL ma_test.cc)

set(tests
  gmi_test
  mds_test
  ma_test)

foreach(test ${tests})
  target_link_libraries(${test}
    ${MA_LIBS} 
    ${MDS_LIBS}
    ${GMI_SIM_LIBS})
endforeach()
