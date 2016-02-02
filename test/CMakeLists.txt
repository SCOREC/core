if(BUILD_IN_TRILINOS)
  include(SCOREC_tribits_exes.cmake)
  return()
endif()

find_package(mds REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(ma REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(gmi_sim REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(apf_sim REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(spr REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(parma REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(ph REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(viz REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(dsp REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(crv REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(lion REQUIRED PATHS ${CMAKE_BINARY_DIR})
find_package(apf_omega_h QUIET PATHS ${CMAKE_BINARY_DIR})

include_directories(
  ${MA_INCLUDE_DIRS}
  ${MDS_INCLUDE_DIRS}
  ${GMI_SIM_INCLUDE_DIRS}
  ${APF_SIM_INCLUDE_DIRS}
  ${SPR_INCLUDE_DIRS}
  ${PARMA_INCLUDE_DIRS}
  ${PH_INCLUDE_DIRS}
  ${VIZ_INCLUDE_DIRS}
  ${DSP_INCLUDE_DIRS}
  ${CRV_INCLUDE_DIRS}
  ${LION_INCLUDE_DIRS}
  ${APF_OMEGA_H_INCLUDE_DIRS}
)
set(ALL_LIBS
  ${MA_LIBS}
  ${MDS_LIBS}
  ${GMI_SIM_LIBS}
  ${APF_SIM_LIBS}
  ${SPR_LIBS}
  ${PARMA_LIBS}
  ${PH_LIBS}
  ${VIZ_LIBS}
  ${DSP_LIBS}
  ${CRV_LIBS}
  ${LION_LIBS}
  ${APF_OMEGA_H_LIBS}
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

setup_exe(mdlConvert mdlConvert.cc)
setup_exe(mkmodel mkmodel.cc)
setup_exe(ma_insphere ma_insphere.cc)
setup_exe(ma_test ma_test.cc)
setup_exe(convert convert.cc)
setup_exe(ph_convert ../phasta/ph_convert.cc)
if(SIMMODSUITE_SimAdvMeshing_FOUND)
  setup_exe(generate generate.cc)
endif()
setup_exe(snap snap.cc)
setup_exe(threshold ../phasta/threshold.cc)
setup_exe(uniform uniform.cc)
setup_exe(refine2x refine2x.cc)
setup_exe(fixshape fixshape.cc)
setup_exe(fixDisconnected fixDisconnected.cc)
setup_exe(tetrahedronize tetrahedronize.cc)
setup_exe(fixlayer fixlayer.cc)
setup_exe(render render.cc)
setup_exe(renderClass renderClass.cc)
setup_exe(render_ascii render_ascii.cc)
setup_exe(render_binary render_binary.cc)
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
setup_exe(vtxElmMixedBalance vtxElmMixedBalance.cc)
setup_exe(vtxEdgeElmBalance vtxEdgeElmBalance.cc)
setup_exe(split split.cc)
setup_exe(zsplit zsplit.cc)
setup_exe(repartition repartition.cc)
setup_exe(chef ../phasta/chef.cc)
if(ENABLE_ZOLTAN)
  setup_exe(ptnParma ptnParma.cc)
endif()
if(PHASTA_CHEF_ENABLED)
  setup_exe(chefStream ../phasta/chefStream.cc)
endif()
setup_exe(cut_interface ../phasta/cut_interface.cc)
setup_exe(chefReadUrPrep ../phasta/readUrPrep.cc)
setup_exe(ghost ghost.cc)
setup_exe(1d 1d.cc)
setup_exe(serialize serialize.cc)
setup_exe(fusion fusion.cc)
setup_exe(fusion2 fusion2.cc)
setup_exe(fusion3 fusion3.cc)
setup_exe(newdim newdim.cc)
setup_exe(construct construct.cc)
setup_exe(shapefun shapefun.cc)
setup_exe(shapefun2 shapefun2.cc)
setup_exe(bezierElevation bezierElevation.cc)
setup_exe(bezierMesh bezierMesh.cc)
setup_exe(bezierMisc bezierMisc.cc)
setup_exe(bezierRefine bezierRefine.cc)
setup_exe(bezierSubdivision bezierSubdivision.cc)
setup_exe(bezierValidity bezierValidity.cc)

if (ENABLE_VIZ)
  setup_exe(viz_test viz.cc)
endif()
setup_exe(from_neper neper.cc)
setup_exe(from_gmsh gmsh.cc)
setup_exe(from_ugrid ugrid.cc)
setup_exe(from_ansys ansys.cc)
setup_exe(eigen_test eigen_test.cc)
setup_exe(curvetest curvetest.cc)
setup_exe(curve_to_bezier curve_to_bezier.cc)
setup_exe(integrate integrate.cc)
setup_exe(align align.cc)
setup_exe(moving moving.cc)
setup_exe(graphdist graphdist.cc)
setup_exe(nektar_align nektar_align.cc)
setup_exe(reposition reposition.cc)
setup_exe(qr qr.cc)
setup_exe(base64 base64.cc)
setup_exe(test_AD test_AD.cc)
if(ENABLE_OMEGA_H)
  setup_exe(from_omega_h from_omega_h.cc)
  setup_exe(to_omega_h to_omega_h.cc)
endif()

if(IS_TESTING)
  include(testing.cmake)
endif()

#executables which call only Simmetrix functions
include_directories(${PROJECT_BINARY_DIR}/gmi_sim)
setup_exe(simTranslate simTranslate.cc)
if(SIM_PARASOLID)
setup_exe(simSpjToSmd simSpjToSmd.cc)
endif()
