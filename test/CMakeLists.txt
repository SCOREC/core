find_package(mds REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(ma REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(spr REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(parma REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(ph REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(apf_mpas REQUIRED PATHS ${CMAKE_BINARY_DIR})

include_directories(
  ${MA_INCLUDE_DIRS}
  ${MDS_INCLUDE_DIRS}
  ${SPR_INCLUDE_DIRS}
  ${PARMA_INCLUDE_DIRS}
  ${PH_INCLUDE_DIRS}
  ${APF_MPAS_INCLUDE_DIRS})
set(ALL_LIBS
  ${MA_LIBS}
  ${MDS_LIBS}
  ${SPR_LIBS}
  ${PARMA_LIBS}
  ${PH_LIBS}
  ${APF_MPAS_LIBS})

if(NOT IS_TESTING)
  set(MYEXCLUDE EXCLUDE_FROM_ALL)
else()
  set(MYEXCLUDE)
endif()

macro(setup_exe exename srcname)
  add_executable(${exename} ${MYEXCLUDE} ${srcname})
  target_link_libraries(${exename} ${ALL_LIBS})
  if(IS_TESTING)
    install(TARGETS ${exename} DESTINATION bin)
  endif()
endmacro(setup_exe)

setup_exe(mdlConvert mdlConvert.c)
setup_exe(mkmodel mkmodel.cc)
setup_exe(ma_insphere ma_insphere.cc)
setup_exe(ma_test ma_test.cc)
setup_exe(chef ../phasta/chef.cc)
setup_exe(uniform uniform.cc)
setup_exe(refine2x refine2x.cc)
setup_exe(fixshape fixshape.cc)
setup_exe(tetrahedronize tetrahedronize.cc)
setup_exe(fixlayer fixlayer.cc)
setup_exe(render render.cc)
setup_exe(verify verify.cc)
setup_exe(prismCodeMatch ../ma/prismCodeMatch.cc)
setup_exe(pyramidCodeMatch ../ma/pyramidCodeMatch.cc)
setup_exe(spr_test spr_test.cc)
setup_exe(describe describe.cc)
setup_exe(balance balance.cc)
setup_exe(gap gap.cc)
setup_exe(elmBalance elmBalance.cc)
setup_exe(vtxBalance vtxBalance.cc)
setup_exe(edgeBalance edgeBalance.cc)
setup_exe(vtxElmBalance vtxElmBalance.cc)
setup_exe(vtxEdgeElmBalance vtxEdgeElmBalance.cc)
setup_exe(split split.cc)
setup_exe(zsplit zsplit.cc)
setup_exe(ghost ghost.cc)
if(ENABLE_MPAS)
  setup_exe(mpas_read ../mpas/mpas_read.cc)
  setup_exe(mpas_write ../mpas/mpas_write.cc)
endif()
setup_exe(hps hps.cc)
setup_exe(fusion fusion.cc)
setup_exe(fusion2 fusion2.cc)
setup_exe(newdim newdim.cc)
setup_exe(construct construct.cc)
setup_exe(shapefun shapefun.cc)
#setup_exe(viz viz.cc)

include_directories("/lore/diamog/core/prefix/include")
setup_exe(viz viz.cc)
target_link_libraries(viz "/lore/diamog/core/prefix/lib/libmilo.a")

if(IS_TESTING)
  include(testing.cmake)
endif()
