find_package(mds REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(ma REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(spr REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(parma REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(ph REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(apf_mpas REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(viz REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(dsp REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(crv REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(mth REQUIRED PATHS ${CMAKE_BINARY_DIR})

include_directories(
  ${MA_INCLUDE_DIRS}
  ${MDS_INCLUDE_DIRS}
  ${SPR_INCLUDE_DIRS}
  ${PARMA_INCLUDE_DIRS}
  ${PH_INCLUDE_DIRS}
  ${APF_MPAS_INCLUDE_DIRS}
  ${VIZ_INCLUDE_DIRS}
  ${DSP_INCLUDE_DIRS}
  ${CRV_INCLUDE_DIRS}
  ${MTH_INCLUDE_DIRS}
)
set(ALL_LIBS
  ${MA_LIBS}
  ${MDS_LIBS}
  ${SPR_LIBS}
  ${PARMA_LIBS}
  ${PH_LIBS}
  ${APF_MPAS_LIBS}
  ${VIZ_LIBS}
  ${DSP_LIBS}
  ${CRV_LIBS}
  ${MTH_LIBS}
)

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
setup_exe(threshold ../phasta/threshold.cc)
setup_exe(uniform uniform.cc)
setup_exe(refine2x refine2x.cc)
setup_exe(fixshape fixshape.cc)
setup_exe(fixDisconnected fixDisconnected.cc)
setup_exe(tetrahedronize tetrahedronize.cc)
setup_exe(fixlayer fixlayer.cc)
setup_exe(render render.cc)
setup_exe(renderClass renderClass.cc)
setup_exe(verify verify.cc)
setup_exe(prismCodeMatch ../ma/prismCodeMatch.cc)
setup_exe(pyramidCodeMatch ../ma/pyramidCodeMatch.cc)
setup_exe(spr_test spr_test.cc)
setup_exe(describe describe.cc)
setup_exe(balance balance.cc)
if(ENABLE_ZOLTAN)
  setup_exe(zbalance zbalance.cc)
endif()
setup_exe(gap gap.cc)
setup_exe(elmBalance elmBalance.cc)
setup_exe(vtxBalance vtxBalance.cc)
setup_exe(vtxElmBalance vtxElmBalance.cc)
setup_exe(vtxEdgeElmBalance vtxEdgeElmBalance.cc)
if(ENABLE_THREADS)
  setup_exe(split split.cc)
  setup_exe(zsplit zsplit.cc)
  setup_exe(chef ../phasta/chef.cc)
  if(ENABLE_ZOLTAN)
    setup_exe(ptnParma ptnParma.cc)
  endif()
else()
  setup_exe(split split_nothread.cc)
  setup_exe(zsplit zsplit_nothread.cc)
  setup_exe(chef ../phasta/chef_nothread.cc)
  if(ENABLE_ZOLTAN)
    setup_exe(ptnParma_nothread ptnParma_nothread.cc)
  endif()
  if(PHASTA_CHEF_ENABLED)
    setup_exe(chefStream ../phasta/chefStream.cc)
  endif()
endif()
setup_exe(ghost ghost.cc)
if(ENABLE_MPAS)
  setup_exe(mpas_read ../mpas/mpas_read.cc)
  setup_exe(mpas_write ../mpas/mpas_write.cc)
endif()
setup_exe(1d 1d.cc)
setup_exe(serialize serialize.cc)
setup_exe(fusion fusion.cc)
setup_exe(fusion2 fusion2.cc)
setup_exe(fusion3 fusion3.cc)
setup_exe(newdim newdim.cc)
setup_exe(construct construct.cc)
setup_exe(shapefun shapefun.cc)
setup_exe(shapefun2 shapefun2.cc)
setup_exe(bezier1 bezier1.cc)
setup_exe(bezier2 bezier2.cc)
setup_exe(bezier3 bezier3.cc)

if (ENABLE_VIZ)
  setup_exe(viz_test viz.cc)
endif()
setup_exe(qr_test qr_test.cc)
setup_exe(from_neper neper.cc)
setup_exe(from_ansys ansys.cc)
setup_exe(eigen_test eigen_test.cc)
setup_exe(integrate integrate.cc)
setup_exe(align align.cc)
setup_exe(moving moving.cc)
setup_exe(graphdist graphdist.cc)
setup_exe(nektar_align nektar_align.cc)
setup_exe(reposition reposition.cc)
setup_exe(qr qr.cc)

if(IS_TESTING)
  include(testing.cmake)
endif()
