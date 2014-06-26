find_package(mds REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(ma REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(spr REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(parma REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(ph REQUIRED PATHS ${CMAKE_BINARY_DIR})

include_directories(
  ${MA_INCLUDE_DIRS}
  ${MDS_INCLUDE_DIRS}
  ${SPR_INCLUDE_DIRS}
  ${PARMA_INCLUDE_DIRS}
  ${PH_INCLUDE_DIRS})

add_executable(mdlConvert EXCLUDE_FROM_ALL mdlConvert.c)
add_executable(mkmodel EXCLUDE_FROM_ALL mkmodel.cc)
add_executable(mds_test EXCLUDE_FROM_ALL mds_test.cc)
add_executable(ma_test EXCLUDE_FROM_ALL ma_test.cc)
add_executable(ph_test EXCLUDE_FROM_ALL ph_test.cc)
add_executable(uniform EXCLUDE_FROM_ALL uniform.cc)
add_executable(fixshape EXCLUDE_FROM_ALL fixshape.cc)
add_executable(tetrahedronize EXCLUDE_FROM_ALL tetrahedronize.cc)
add_executable(render EXCLUDE_FROM_ALL render.cc)
add_executable(verify EXCLUDE_FROM_ALL verify.cc)
add_executable(prismCodeMatch EXCLUDE_FROM_ALL ../ma/prismCodeMatch.cc)
add_executable(spr_test EXCLUDE_FROM_ALL spr_test.cc)
add_executable(describe EXCLUDE_FROM_ALL describe.cc)
add_executable(balance EXCLUDE_FROM_ALL balance.cc)
add_executable(split EXCLUDE_FROM_ALL split.cc)
add_executable(zsplit EXCLUDE_FROM_ALL zsplit.cc)
add_executable(ghost EXCLUDE_FROM_ALL ghost.cc)
add_executable(bfs EXCLUDE_FROM_ALL bfs.cc)
add_executable(hps EXCLUDE_FROM_ALL hps.cc)
add_executable(fan EXCLUDE_FROM_ALL fan.cc)

set(tests
  mdlConvert
  mkmodel
  mds_test
  ma_test
  ph_test
  uniform
  fixshape
  tetrahedronize
  render
  verify
  prismCodeMatch
  spr_test
  describe
  balance
  split
  zsplit
  ghost
  bfs
  hps
  fan)

foreach(test ${tests})
  target_link_libraries(${test}
    ${MA_LIBS}
    ${MDS_LIBS}
    ${SPR_LIBS}
    ${PARMA_LIBS}
    ${PH_LIBS})
endforeach()
