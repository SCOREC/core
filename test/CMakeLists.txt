find_package(mds REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(ma REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(gmi_sim REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(apf_sim REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(spr REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(parma REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(ph REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(apf_mpas REQUIRED PATHS ${CMAKE_BINARY_DIR})

include_directories(
  ${MA_INCLUDE_DIRS}
  ${MDS_INCLUDE_DIRS}
  ${GMI_SIM_INCLUDE_DIRS}
  ${APF_SIM_INCLUDE_DIRS}
  ${SPR_INCLUDE_DIRS}
  ${PARMA_INCLUDE_DIRS}
  ${PH_INCLUDE_DIRS}
  ${APF_MPAS_INCLUDE_DIRS})

macro(buildExe exe src)
  add_executable(${exe} EXCLUDE_FROM_ALL ${src})
  set(tests ${tests} ${exe})
endmacro(buildExe)

set(tests)
buildExe(mdlConvert mdlConvert.cc)
buildExe(mkmodel mkmodel.cc)
buildExe(mds_test mds_test.cc)
buildExe(ma_test ma_test.cc)
buildExe(convert convert.cc)
if(SIMMODSUITE_SimAdvMeshing_FOUND)
  buildExe(generate generate.cc)
endif()
buildExe(snap snap.cc)
buildExe(ph_test ph_test.cc)
buildExe(uniform uniform.cc)
buildExe(fixshape fixshape.cc)
buildExe(tetrahedronize tetrahedronize.cc)
buildExe(fixlayer fixlayer.cc)
buildExe(render render.cc)
buildExe(verify verify.cc)
buildExe(prismCodeMatch ../ma/prismCodeMatch.cc)
buildExe(pyramidCodeMatch ../ma/pyramidCodeMatch.cc)
buildExe(spr_test spr_test.cc)
buildExe(describe describe.cc)
buildExe(balance balance.cc)
buildExe(split split.cc)
buildExe(zsplit zsplit.cc)
buildExe(ghost ghost.cc)
buildExe(mpas_read mpas_read.cc)
buildExe(mpas_write mpas_write.cc)
buildExe(hps hps.cc)
buildExe(fan fan.cc)

foreach(test ${tests})
  target_link_libraries(${test}
    ${MA_LIBS}
    ${MDS_LIBS}
    ${GMI_SIM_LIBS}
    ${APF_SIM_LIBS}
    ${SPR_LIBS}
    ${PARMA_LIBS}
    ${PH_LIBS}
    ${APF_MPAS_LIBS})
endforeach()
