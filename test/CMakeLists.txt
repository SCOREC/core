find_package(mds REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(ma REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(gmi_sim REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(apf_sim REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(spr REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(parma REQUIRED PATHS ${CMAKE_BINARY_DIR})

include_directories(
  ${MA_INCLUDE_DIRS}
  ${MDS_INCLUDE_DIRS}
  ${GMI_SIM_INCLUDE_DIRS}
  ${APF_SIM_INCLUDE_DIRS}
  ${SPR_INCLUDE_DIRS}
  ${PARMA_INCLUDE_DIRS})

add_executable(gmi_test EXCLUDE_FROM_ALL gmi_test.cc)
add_executable(mds_test EXCLUDE_FROM_ALL mds_test.cc)
add_executable(ma_test EXCLUDE_FROM_ALL ma_test.cc)
add_executable(convert EXCLUDE_FROM_ALL convert.cc)
add_executable(uniform EXCLUDE_FROM_ALL uniform.cc)
add_executable(render EXCLUDE_FROM_ALL render.cc)
add_executable(verify EXCLUDE_FROM_ALL verify.cc)
add_executable(spr_test EXCLUDE_FROM_ALL spr_test.cc)
add_executable(separate_bl EXCLUDE_FROM_ALL separate_bl.cc)
add_executable(describe EXCLUDE_FROM_ALL describe.cc)
add_executable(split EXCLUDE_FROM_ALL split.cc)

set(tests
  gmi_test
  mds_test
  ma_test
  convert
  uniform
  render
  verify
  spr_test
  separate_bl
  describe
  split)

foreach(test ${tests})
  target_link_libraries(${test}
    ${MA_LIBS}
    ${MDS_LIBS}
    ${GMI_SIM_LIBS}
    ${APF_SIM_LIBS}
    ${SPR_LIBS}
    ${PARMA_LIBS})
endforeach()
