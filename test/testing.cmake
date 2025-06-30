set(MESHES ""
    CACHE STRING
    "Extracted http://scorec.rpi.edu/pumi/pumi_test_meshes.tar.gz")
function(mpi_test TESTNAME PROCS EXE)
  if(SCOREC_NO_MPI)
    if(${PROCS} EQUAL "1")
      add_test(
        NAME ${TESTNAME}
        COMMAND ${VALGRIND} ${VALGRIND_ARGS} ${EXE} ${ARGN}
      )
    endif()
  else()
    add_test(
      NAME ${TESTNAME}
      COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} ${PROCS}
        ${VALGRIND} ${VALGRIND_ARGS} ${EXE} ${ARGN}
    )
    set_tests_properties(${TESTNAME} PROPERTIES PROCESSORS ${PROCS})
  endif()
endfunction(mpi_test)

# USAGE: set_test_depends(TESTS test1 [test2...] DEPENDS test3 [test4...])
function(set_test_depends)
  cmake_parse_arguments(SET_TEST_DEPENDS "" "" "TESTS;DEPENDS" ${ARGN})
  if (NOT DEFINED SET_TEST_DEPENDS_TESTS OR NOT DEFINED SET_TEST_DEPENDS_DEPENDS)
    return()
  endif()
  if(SCOREC_NO_MPI)
    # Check for test existence as it may be a multiproc test.
    if(TEST "${TESTNAME}")
      set_tests_properties(${SET_TEST_DEPENDS_TESTS} PROPERTIES
        DEPENDS "${SET_TEST_DEPENDS_DEPENDS}")
    endif()
  else()
    # Don't check if test exists because it's more likely a typo.
    set_tests_properties(${SET_TEST_DEPENDS_TESTS} PROPERTIES
      DEPENDS "${SET_TEST_DEPENDS_DEPENDS}")
  endif()
endfunction()

mpi_test(shapefun 1 ./shapefun)
mpi_test(shapefun2 1 ./shapefun2)
mpi_test(bezierElevation 1 ./bezierElevation)
mpi_test(bezierMesh 1 ./bezierMesh)
mpi_test(bezierMisc 1 ./bezierMisc)
mpi_test(bezierRefine 1 ./bezierRefine)
mpi_test(bezierSubdivision 1 ./bezierSubdivision)
mpi_test(bezierValidity 1 ./bezierValidity)
mpi_test(ma_analytic 1 ./ma_test_analytic_model)

if(ENABLE_ZOLTAN)
mpi_test(print_pumipic_partion 1
         ./print_pumipic_partition
         ${MESHES}/cube/cube.dmg
         ${MESHES}/cube/pumi11/cube.smb
         4
         pumipic_cube
         )
endif()

mpi_test(align 1 ./align)
mpi_test(eigen_test 1 ./eigen_test)
mpi_test(integrate 1 ./integrate)
mpi_test(qr_test 1 ./qr)
mpi_test(swapDoubles 1 ./swapDoubles)
mpi_test(base64 1 ./base64)
mpi_test(tensor_test 1 ./tensor)
mpi_test(verify_convert 1 ./verify_convert)
mpi_test(test_integrator 1
         ./test_integrator
         "${MESHES}/cube/cube.dmg"
         "${MESHES}/cube/pumi11/cube.smb"
         )
mpi_test(test_matrix_gradient 1
         ./test_matrix_gradient
         "${MESHES}/cube/cube.dmg"
         "${MESHES}/cube/pumi11/cube.smb"
         )

mpi_test(modelInfo_dmg 1
  ./modelInfo
  "${MESHES}/cube/cube.dmg")
if(ENABLE_SIMMETRIX)
  mpi_test(in_closure_of 1
    ./inClosureOf_test
    "${MESHES}/cube/cube.smd")
  mpi_test(modelInfo_smd 1
    ./modelInfo
    "${MESHES}/cube/cube.smd")
  mpi_test(highorder_sizefield 1
    ./highOrderSizeFields
    "${MESHES}/cube/cube.smd"
    "${MESHES}/cube/pumi11/cube.smb")
endif(ENABLE_SIMMETRIX)

if(ENABLE_SIMMETRIX)
  set(GXT smd)
else()
  set(GXT dmg)
endif()

set(MDIR ${MESHES}/phasta/dg)

if(ENABLE_SIMMETRIX AND SIM_PARASOLID AND SIMMODSUITE_SimAdvMeshing_FOUND)
  set(MDIR ${MESHES}/phasta/BL_query)
  mpi_test(chef-BL_query 4 ${CMAKE_CURRENT_BINARY_DIR}/chef
    WORKING_DIRECTORY ${MDIR}/run_case)
  add_test(NAME chef-BL_query-diff
    COMMAND diff -r run_case/4-procs_case/ good_case/4-procs_case
    WORKING_DIRECTORY ${MDIR})
  set_test_depends(TESTS chef-BL_query-diff DEPENDS chef-BL_query)
  set(MDIR ${MESHES}/simExtrusionInfo)
  mpi_test(convertExtrudedRoots 1 ${CMAKE_CURRENT_BINARY_DIR}/convert
    --model-face-root=${MDIR}/ExtruRootID.txt
    --native-model=${MDIR}/geom.xmt_txt
    ${MDIR}/geom.smd ${MDIR}/geom.sms ${MDIR}/mdsMesh.smb
    WORKING_DIRECTORY ${MDIR})
  add_test(NAME convertExtrudedRoots_diff_cnn
    COMMAND diff -r geom.cnn geom_expected.cnn
    WORKING_DIRECTORY ${MDIR})
  add_test(NAME convertExtrudedRoots_diff_crd
    COMMAND diff -r geom.crd geom_expected.crd
    WORKING_DIRECTORY ${MDIR})
  set_test_depends(TESTS
    convertExtrudedRoots_diff_cnn convertExtrudedRoots_diff_crd
    DEPENDS convertExtrudedRoots
  )
endif()

if(ENABLE_SIMMETRIX AND SIM_PARASOLID AND SIMMODSUITE_SimAdvMeshing_FOUND)
  if(SIM_DOT_VERSION VERSION_GREATER 12.0.171000)
    set(MDIR ${MESHES}/faceExtrusion)
    mpi_test(rm_extrusion 1
      ${CMAKE_CURRENT_BINARY_DIR}/rm_extrusion
      ${MDIR}/plate.x_t
      ${MDIR}/extrusion.sms
      ${MDIR}/extrusion_noatts.sms)
  endif()
  set(MDIR ${MESHES}/phasta/BL_query/cut_and_partition)
  add_test(NAME cut_interface_sim
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/cut_interface
   "${MDIR}/model_nat.x_t"
   "${MDIR}/model.smd"
   "${MDIR}/mesh.sms"
   "${MDIR}/mesh_cut.sms"
   WORKING_DIRECTORY ${MDIR})
  if(SIM_DOT_VERSION VERSION_GREATER 11.0.170826)
    mpi_test(partition_sim 4
     ${CMAKE_CURRENT_BINARY_DIR}/sim_part
     "${MDIR}/model_nat.x_t"
     "${MDIR}/model.smd"
     "${MDIR}/mesh_cut.sms"
     4
     WORKING_DIRECTORY ${MDIR})
    set_test_depends(TESTS partition_sim DEPENDS cut_interface_sim)
  endif()
  if(SIMMODSUITE_SimAdvMeshing_FOUND)
    add_test(NAME countBL_cut_mesh
      COMMAND ${CMAKE_CURRENT_BINARY_DIR}/sim_countBL
     "${MDIR}/model_nat.x_t"
     "${MDIR}/model.smd"
     "${MDIR}/mesh_cut.sms"
     3504
     WORKING_DIRECTORY ${MDIR})
    set_test_depends(TESTS countBL_cut_mesh DEPENDS cut_interface_sim)
    if(SIM_DOT_VERSION VERSION_GREATER 11.0.170826)
      mpi_test(countBL_part_mesh 4
       ${CMAKE_CURRENT_BINARY_DIR}/sim_countBL
       "${MDIR}/model_nat.x_t"
       "${MDIR}/model.smd"
       "${MDIR}/outmesh_4_parts.sms"
       3504
       WORKING_DIRECTORY ${MDIR})
      set_test_depends(TESTS countBL_part_mesh DEPENDS partition_sim)
    endif()
  endif()
endif(ENABLE_SIMMETRIX AND SIM_PARASOLID AND SIMMODSUITE_SimAdvMeshing_FOUND)

set(MDIR ${MESHES}/phasta/loopDriver)
if(ENABLE_ZOLTAN AND ENABLE_SIMMETRIX AND PCU_COMPRESS AND SIM_PARASOLID
    AND SIMMODSUITE_SimAdvMeshing_FOUND)
  mpi_test(ph_adapt 1
    ${CMAKE_CURRENT_BINARY_DIR}/ph_adapt
    "${MDIR}/model.smd"
    "${MDIR}/mesh_.smb"
    WORKING_DIRECTORY ${MDIR})
endif()


if(ENABLE_ZOLTAN)
  mpi_test(pumi3d-1p 4
    ./test_pumi
    ${MESHES}/pumi/3d-1p/model.dmg
    ${MESHES}/pumi/3d-1p/part.smb
    out.smb 1 0)
endif()
if(ENABLE_OMEGA_H)
  mpi_test(mdsToOmega 1
    ./smb2osh
    ${MESHES}/cube/cube.dmg
    ${MESHES}/cube/pumi670/cube.smb
    cube.osh)
  mpi_test(omegaToMds 1
    ./osh2smb
    cube.osh
    ${MESHES}/cube/cube.dmg
    converted.smb)
endif()
mpi_test(test_scaling 1
  ./test_scaling
  ${MESHES}/cube/cube.dmg
  ${MESHES}/cube/pumi670/cube.smb)
mpi_test(render 1
  ./render
  ${MESHES}/cube/cube.dmg
  ${MESHES}/cube/pumi11/cube.smb
  render_test)
mpi_test(render_ascii 1
  ./render_ascii
  ${MESHES}/cube/cube.dmg
  ${MESHES}/cube/pumi11/cube.smb
  render_ascii_test)
mpi_test(field_io 1
  ./field_io
  ${MESHES}/cube/cube.dmg
  ${MESHES}/cube/pumi11/cube.smb)
mpi_test(reorder_serial 1
  ./reorder
  ${MESHES}/cube/cube.dmg
  ${MESHES}/cube/pumi7k/cube.smb
  cube_bfs.smb)
mpi_test(create_misCube 1
  ./create_mis
  ${MESHES}/cube/cube.dmg
  ${MESHES}/cube/pumi670/cube.smb
  mis_test)
mpi_test(create_misSquare 1
  ./create_mis
  ${MESHES}/square/square.dmg
  ${MESHES}/square/square.smb
  mis_test)

set(MDIR ${MESHES}/matchedNodeElementReader)
mpi_test(matchedNodeElementReader_p1 1
  ./matchedNodeElmReader
  "${MDIR}/model.dmg"
  "${MDIR}/1part/geom3D.cnndt"
  "${MDIR}/1part/geom3D.coord"
  "${MDIR}/1part/geom3D.match"
  "${MDIR}/1part/geom3D.class"
  "${MDIR}/1part/geom3D.fathr"
  "NULL"
  "${MDIR}/1part/geom3DHead.cnn"
  "geom.dmg" "geom.smb")

mpi_test(matchedNodeElementReader_p4 4
  ./matchedNodeElmReader
  "${MDIR}/model.dmg"
  "${MDIR}/4part/geom3D.cnndt"
  "${MDIR}/4part/geom3D.coord"
  "${MDIR}/4part/geom3D.match"
  "${MDIR}/4part/geom3D.class"
  "${MDIR}/4part/geom3D.fathr"
  "NULL"
  "${MDIR}/4part/geom3DHead.cnn"
  "geom.dmg" "geom.smb")

set(MDIR ${MESHES}/gmsh)
mpi_test(gmshv2TwoQuads 1
  ./from_gmsh
  "none"
  "${MDIR}/twoQuads.msh"
  "${MDIR}/twoQuads.smb"
  "${MDIR}/twoQuads.dmg")

set(MDIR ${MESHES}/gmsh/v4)
mpi_test(gmshV4AirFoil 1
  ./from_gmsh
  "none"
  "${MDIR}/AirfoilDemo.msh"
  "${MDIR}/AirfoilDemo.smb"
  "${MDIR}/AirfoilDemo.dmg")
add_test(NAME gmshV4AirFoil_dmgDiff
  COMMAND diff -r ${MDIR}/AirfoilDemo.dmg AirfoilDemo_gold.dmg
  WORKING_DIRECTORY ${MDIR})
set_test_depends(TESTS gmshV4AirFoil_dmgDiff DEPENDS gmshV4AirFoil)

set(MDIR ${MESHES}/ugrid)
mpi_test(naca_ugrid 2
  ./from_ugrid
  "${MDIR}/inviscid_egg.b8.ugrid"
  "${MDIR}/naca.dmg"
  "${MDIR}/2/"
  "2")
mpi_test(inviscid_ugrid 4
  ./from_ugrid
  "${MDIR}/inviscid_egg.b8.ugrid"
  "${MDIR}/inviscid_egg.dmg"
  "${MDIR}/4/"
  "4")
mpi_test(inviscid_ghost 4
  ./ghost
  "${MDIR}/inviscid_egg.dmg"
  "${MDIR}/4/"
  "${MDIR}/vis")
set_test_depends(TESTS inviscid_ghost DEPENDS inviscid_ugrid)

set(MDIR ${SMOKE_TEST_MESHES}/pipe)
if(ENABLE_SIMMETRIX)
  mpi_test(convert 1
    ./convert
    "${MDIR}/pipe.smd"
    "${MDIR}/pipe.sms"
    "pipe.smb")
else()
  file(COPY "${MDIR}/pipe0.smb" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()
mpi_test(verify_serial 1
  ./verify
  "${MDIR}/pipe.${GXT}"
  "pipe.smb")
if(ENABLE_SIMMETRIX)
  set_test_depends(TESTS verify_serial DEPENDS convert)
endif()
if(ENABLE_SIMMETRIX AND SIM_PARASOLID)
  mpi_test(convert_2d_quads 1
    ./convert
    "${MESHES}/disk/disk.smd"
    "${MESHES}/disk/disk_quad_mesh.sms"
    "disk_quad_mesh.smb")
  mpi_test(convert_2d_tris 1
    ./convert
    "${MESHES}/disk/disk.smd"
    "${MESHES}/disk/disk_tri_mesh.sms"
    "disk_tri_mesh.smb")
endif()
mpi_test(verify_2nd_order_shape_quads 1
  ./verify_2nd_order_shapes
  "disk_quad_mesh.smb")
if(ENABLE_SIMMETRIX AND SIM_PARASOLID)
  set_test_depends(TESTS verify_2nd_order_shape_quads DEPENDS convert_2d_quads)
else()
  file(COPY "${MESHES}/disk/disk_quad_mesh0.smb" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()
mpi_test(verify_2nd_order_shape_tris 1
  ./verify_2nd_order_shapes
  "disk_tri_mesh.smb")
if(ENABLE_SIMMETRIX AND SIM_PARASOLID)
  set_test_depends(TESTS verify_2nd_order_shape_tris DEPENDS convert_2d_tris)
else()
  file(COPY "${MESHES}/disk/disk_tri_mesh0.smb" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()
mpi_test(uniform_serial 1
  ./uniform
  "${MDIR}/pipe.${GXT}"
  "pipe.smb"
  "pipe_unif.smb")
mpi_test(classifyThenAdapt 1 ./classifyThenAdapt)
if(ENABLE_SIMMETRIX)
  mpi_test(snap_serial 1
    ./snap
    "${MDIR}/pipe.${GXT}"
    "pipe_unif.smb"
    "pipe.smb")
  set_test_depends(TESTS snap_serial DEPENDS uniform_serial)
endif()
if(ENABLE_ZOLTAN)
  mpi_test(ma_serial 1
    ./ma_test
    "${MDIR}/pipe.${GXT}"
    "pipe.smb")
  mpi_test(aniso_ma_serial 1
    ./aniso_ma_test
    "${MESHES}/cube/cube.dmg"
    "${MESHES}/cube/pumi670/cube.smb"
    "0")
  mpi_test(aniso_ma_serial_log_interpolation 1
    ./aniso_ma_test
    "${MESHES}/cube/cube.dmg"
    "${MESHES}/cube/pumi670/cube.smb"
    "1")
  mpi_test(torus_ma_parallel 4
    ./torus_ma_test
    "${MESHES}/torus/torus.dmg"
    "${MESHES}/torus/4imb/torus.smb")
  mpi_test(ma_2dLayersOff 1
    ./ma_test
    "${MESHES}/2dlayersNoAdapt/model.dmg"
    "${MESHES}/2dlayersNoAdapt/mesh.smb"
    "doNotAdapt" "8")
  set_tests_properties(ma_2dLayersOff PROPERTIES
    PASS_REGULAR_EXPRESSION "number of triangle 18698")
endif()
mpi_test(tet_serial 1
  ./tetrahedronize
  "${MDIR}/pipe.${GXT}"
  "pipe.smb"
  "tet.smb")
if(ENABLE_SIMMETRIX AND SIM_PARASOLID)
  mpi_test(test_residual_error_estimate 1
    ./residualErrorEstimation_test
    "${MESHES}/electromagnetic/fichera_geomSim.smd"
    "${MESHES}/electromagnetic/fichera_1k.smb")
endif()
if(PCU_COMPRESS)
  set(MESHFILE "bz2:pipe_2_.smb")
else()
  set(MESHFILE "pipe_2_.smb")
endif()
mpi_test(split_2 2
  ./split
  "${MDIR}/pipe.${GXT}"
  "pipe.smb"
  ${MESHFILE}
  2)
mpi_test(collapse_2 2
  ./collapse
  "${MDIR}/pipe.${GXT}"
  ${MESHFILE}
  pipe_p1_.smb
  2)
if(ENABLE_SIMMETRIX)
  set_test_depends(TESTS split_2 collapse_2 tet_serial DEPENDS convert)
endif()
if(ENABLE_METIS)
  mpi_test(msplit_2 2
    ./msplit
    "${MDIR}/pipe.${GXT}" "pipe.smb"
    "pipe_m2_.smb"
    2
  )
  mpi_test(msplit_3 3
    ./msplit
    "${MDIR}/pipe.${GXT}" "pipe.smb"
    "pipe_m3_.smb"
    3
  )
  mpi_test(msplit_6 6
    ./msplit
    "${MDIR}/pipe.${GXT}" "pipe_m2_.smb"
    "pipe_m6_.smb"
    3
  )
  set_test_depends(TESTS msplit_6 DEPENDS msplit_2)
  mpi_test(mbalanceEmpty 4
    ./mbalanceEmpty
    "${MDIR}/pipe.dmg" "pipe.smb" 1
    "pipe_mbe_.smb"
  )
  if(ENABLE_SIMMETRIX)
    set_test_depends(
      TESTS msplit_2 msplit_3 msplit_6 mbalanceEmpty
      DEPENDS convert
    )
  endif()
endif()
if(ENABLE_ZOLTAN)
  mpi_test(refineX 2
    ./refine2x
    "${MDIR}/pipe.${GXT}"
    ${MESHFILE}
    0
    "refXpipe/")
  set_test_depends(TESTS refineX DEPENDS split_2)
  mpi_test(split_4 4
    ./zsplit
    "${MDIR}/pipe.${GXT}"
    ${MESHFILE}
    "pipe_4_.smb"
    2)
else()
  mpi_test(split_4 4
    ./split
    "${MDIR}/pipe.${GXT}"
    ${MESHFILE}
    "pipe_4_.smb"
    2)
endif()
set_test_depends(TESTS split_4 DEPENDS split_2)
mpi_test(pipe_condense 4
  ./serialize
  "${MDIR}/pipe.${GXT}"
  "pipe_4_.smb"
  "pipe_2.smb"
  2)
mpi_test(verify_parallel 4
  ./verify
  "${MDIR}/pipe.${GXT}"
  "pipe_4_.smb")
mpi_test(vtxElmMixedBalance 4
  ./vtxElmMixedBalance
  "${MDIR}/pipe.${GXT}"
  "pipe_4_.smb")
if(ENABLE_ZOLTAN)
  mpi_test(ma_parallel 4
    ./ma_test
    "${MDIR}/pipe.${GXT}"
    "pipe_4_.smb")
  mpi_test(tet_parallel 4
    ./tetrahedronize
    "${MDIR}/pipe.${GXT}"
    "pipe_4_.smb"
    "tet.smb")
  set_test_depends(TESTS ma_parallel tet_parallel DEPENDS split_4)
endif()
mpi_test(fieldReduce 4
  ./fieldReduce
  "${MDIR}/pipe.${GXT}"
  "pipe_4_.smb")
set_test_depends(TESTS pipe_condense verify_parallel fieldReduce
  verify_parallel vtxElmMixedBalance DEPENDS split_4)

set(MDIR ${MESHES}/torus)
mpi_test(reorder 4
  ./reorder
  "${MDIR}/torus.dmg"
  "${MDIR}/4imb/torus.smb"
  "torusBfs4p/")
mpi_test(balance 4
  ./balance
  "${MDIR}/torus.dmg"
  "${MDIR}/4imb/torus.smb"
  "${MDIR}/torusBal4p/")
mpi_test(gap 4
  ./gap
  "${MDIR}/torus.dmg"
  "${MDIR}/torusBal4p/"
  "1.08"
  "${MDIR}/torusOpt4p/")
set_test_depends(TESTS gap DEPENDS balance)
mpi_test(applyMatrixFunc 1
  ./applyMatrixFunc)
if(ENABLE_ZOLTAN)
  mpi_test(zbalance 4
    ./zbalance
    "${MDIR}/torus.dmg"
    "${MDIR}/4imb/torus.smb"
    "torusZbal4p/")
endif()
if(ENABLE_METIS)
  mpi_test(mbalance 4
    ./mbalance
    "${MDIR}/torus.dmg"
    "${MDIR}/4imb/torus.smb"
    "torusMbal4p/")
endif()
mpi_test(ghostMPAS 4
  ./ghostMPAS
  "${MDIR}/torus.dmg"
  "${MDIR}/4imb/torus.smb"
  "torusGhostEle4p/")
mpi_test(ghostEdge 4
  ./ghostEdge
  "${MDIR}/torus.dmg"
  "${MDIR}/4imb/torus.smb"
  "torusGhostEle4p/")
mpi_test(fixDisconnected 4
  ./fixDisconnected
  "${MDIR}/torus.dmg"
  "${MDIR}/4imb/torus.smb"
  "torusDcFix4p/")
mpi_test(outputcontrol 1
  ./outputcontrol 2)
mpi_test(quality 4
  ./quality
  "${MDIR}/torus.dmg"
  "${MDIR}/4imb/torus.smb"
  .3)
set(MDIR ${MESHES}/airFoilAfosr)
mpi_test(elmBalance 4
  ./elmBalance
  "${MDIR}/afosr.dmg"
  "${MDIR}/4imb/"
  "afosrBal4p/")
mpi_test(vtxBalance 4
  ./vtxBalance
  "${MDIR}/afosr.${GXT}"
  "${MDIR}/4imb/"
  "afosrBal4p/")
mpi_test(vtxEdgeElmBalance 4
  ./vtxEdgeElmBalance
  "${MDIR}/afosr.${GXT}"
  "${MDIR}/4imb/"
  "afosrBal4p/"
  "2"
  "1.10")
mpi_test(vtxElmBalance 4
  ./vtxElmBalance
  "${MDIR}/afosr.dmg"
  "${MDIR}/4imb/"
  "afosrBal4p/")
mpi_test(parmaSerial 1
  ./vtxElmBalance
  "${MESHES}/cube/cube.dmg"
  "${MESHES}/cube/pumi670/cube.smb"
  "cubeBal.smb/")
if(ENABLE_ZOLTAN AND ENABLE_SIMMETRIX AND SIM_PARASOLID)
  set(MDIR ${MESHES}/annular)
  mpi_test(simZBalance_4 4
    ./simZBalance
    "${MDIR}/annular.smd"
    "${MDIR}/annular_4_part.sms")
endif()
set(MDIR ${MESHES}/cube)
if(ENABLE_ZOLTAN)
  mpi_test(ptnParma_cube 4
    ./ptnParma
    "${MDIR}/cube.dmg"
    "${MDIR}/pumi670/cube.smb"
    "ptnParmaCube/"
    "4" "rib" "reptn" "1"
  )
endif()

if(ENABLE_CGNS AND ENABLE_ZOLTAN)
#
# sort of an arbitrary choice
set(numProcs 4)
#
set(CGNSDIR ${MESHES}/cgns/basic)
#
# 2D tests including for mixed cells
#
mpi_test(cgns_2d_1 ${numProcs}
  ./from_cgns
  "${CGNSDIR}/2D/4quads.cgns"
  4quads.smb
  additional)
mpi_test(cgns_2d_2 ${numProcs}
  ./from_cgns
  "${CGNSDIR}/2D/5quad1Tri.cgns"
  5quad1Tri.smb
  additional)
mpi_test(cgns_2d_3 ${numProcs}
  ./from_cgns
  "${CGNSDIR}/2D/5quad2Tri.cgns"
  5quad2Tri.smb
  additional)
mpi_test(cgns_2d_4 ${numProcs}
  ./from_cgns
  "${CGNSDIR}/2D/9tris.cgns"
  9tris.smb
  additional)
#
# 3D tests including for mixed cells
#
mpi_test(cgns_3d_1 ${numProcs}
  ./from_cgns
  "${CGNSDIR}/3D/tets_pyra.cgns"
  tets_pyra.smb
  additional)
mpi_test(cgns_3d_2 ${numProcs}
  ./from_cgns
  "${CGNSDIR}/3D/hexs.cgns"
  hexs.smb
  additional)
#
# 3D BCS tests
#
set(numProcs 5)
#
set(CGNSDIR ${MESHES}/cgns/withBCS/3D)
#
mpi_test(cgns_bcs_1 ${numProcs}
  ./from_cgns
  "${CGNSDIR}/mixed.cgns"
  bcs1.smb
  additional)

mpi_test(cgns_bcs_hex ${numProcs}
  ./from_cgns
  "${CGNSDIR}/8hexs.cgns"
  bcshex.smb
  additional)
#
# 2D BCS tests
#
set(numProcs 4)
#
set(CGNSDIR ${MESHES}/cgns/withBCS/2D)
#
mpi_test(cgns_bcs_2 ${numProcs}
  ./from_cgns
  "${CGNSDIR}/4quads.cgns"
  bcs2.smb
  additional)
#
# 1D BCS tests
#
set(numProcs 3)
#
set(CGNSDIR ${MESHES}/cgns/withBCS/1D)
#
mpi_test(cgns_bcs_3 ${numProcs}
  ./from_cgns
  "${CGNSDIR}/edges.cgns"
  bcs3.smb
  additional)

endif(ENABLE_CGNS AND ENABLE_ZOLTAN)

mpi_test(construct 4
  ./construct
  "${MDIR}/cube.dmg"
  "${MDIR}/pumi7k/4/cube.smb")
mpi_test(constructThenGhost 4
  ./constructThenGhost
  "${MDIR}/cube.dmg"
  "${MDIR}/pumi7k/4/cube.smb")
mpi_test(construct_bottom_up 1
  ./construct_bottom_up
  "${MDIR}/bottom_up_constructed_cube.smb"
  "${MDIR}/cube.dmg")
set(MDIR ${MESHES}/embeddedEdges)
mpi_test(embedded_edges 1
  ./embedded_edges
  "${MDIR}/edges-embedded-in-3D.smb")
set(MDIR ${MESHES}/spr)
mpi_test(spr_3D 4
  ./spr_test
  "${MDIR}/spr.dmg"
  "${MDIR}/quadspr.smb"
  spr3D
  2)
mpi_test(spr_2D 4
  ./spr_test
  "${MDIR}/square.dmg"
  "${MDIR}/square.smb"
  spr2D
  1)
mpi_test(mixedNumbering 4
  ./mixedNumbering
  "${MDIR}/square.dmg"
  "${MDIR}/square.smb"
  out)
set(MDIR ${MESHES}/square)
mpi_test(hierarchic_2p_2D 1
  ./hierarchic
  "${MDIR}/square.dmg"
  "${MDIR}/square.smb"
  2)
mpi_test(hierarchic_3p_2D 1
  ./hierarchic
  "${MDIR}/square.dmg"
  "${MDIR}/square.smb"
  3)
set(MDIR ${MESHES}/cube/pumi24)
mpi_test(hierarchic_2p_3D 1
  ./hierarchic
  "${MDIR}/cube.dmg"
  "${MDIR}/cube.smb"
  2)
set(MDIR ${MESHES}/cube/pumi24)
mpi_test(nedelec 1
  ./nedelecShapes
  "${MDIR}/cube.dmg"
  "${MDIR}/cube.smb")
set(MDIR ${MESHES}/cube/pumi670)
mpi_test(l2_shape_tet_serial 1
  ./L2Shapes
  "${MDIR}/../cube.dmg"
  "${MDIR}/cube.smb")
mpi_test(l2_shape_tet_parallel 4
  ./L2Shapes
  "${MDIR}/../cube.dmg"
  "${MDIR}/4/cube.smb")
mpi_test(h1_shape_serial 1
  ./H1Shapes
  "${MDIR}/../cube.dmg"
  "${MDIR}/cube.smb")
mpi_test(h1_shape_parallel 4
  ./H1Shapes
  "${MDIR}/../cube.dmg"
  "${MDIR}/4/cube.smb")
set(MDIR ${MESHES}/cube/pumi24)
mpi_test(pumiLoadMesh-1p 1
  ./pumiLoadMesh
  ${MDIR}/cube.dmg
  ${MDIR}/cube.smb)
set(MDIR ${MESHES}/cube)
mpi_test(test_verify 4
  ./test_verify
  "${MDIR}/cube.dmg"
  "${MDIR}/pumi7k/4/cube.smb")
set_test_depends(TESTS l2_shape_tet_parallel h1_shape_parallel test_verify
  DEPENDS "construct;constructThenGhost")
set(MDIR ${MESHES}/nonmanifold)
mpi_test(nonmanif_verify 1
  ./verify
  "${MDIR}/nonmanifold.dmg"
  "${MDIR}/nonmanifold.smb")
mpi_test(nonmanif_split 2
  ./split
  "${MDIR}/nonmanifold.dmg"
  "${MDIR}/nonmanifold.smb"
  "nonmanifold_2_.smb"
  2)
mpi_test(nonmanif_verify2 2
  ./verify
  "${MDIR}/nonmanifold.dmg"
  "nonmanifold_2_.smb")
set_test_depends(TESTS nonmanif_verify2 DEPENDS nonmanif_split2)
set(MDIR ${MESHES}/fusion)
mpi_test(mkmodel_fusion 1
  ./mkmodel
  "${MDIR}/fusion.smb"
  "fusionNull.dmg")
mpi_test(mktopomodel_fusion 1
  ./mktopomodel
  "${MDIR}/fusion.smb"
  "fusion.dmg")
mpi_test(split_fusion 2
  ./split
  "fusion.dmg"
  "${MDIR}/fusion.smb"
  "fusion_2_.smb"
  2)
set_test_depends(TESTS split_fusion DEPENDS mktopomodel_fusion)
# the part count mismatch is intentional,
# this test runs on half its procs
if(ENABLE_ZOLTAN)
  mpi_test(adapt_fusion 4
    ./fusion
    "fusion_2_.smb")
  set_test_depends(TESTS adapt_fusion DEPENDS split_fusion)
endif()
mpi_test(fusion_field 2
  ./fusion2)
mpi_test(change_dim 1
  ./newdim)
mpi_test(ma_insphere 1
  ./ma_insphere)
if(ENABLE_SIMMETRIX)
  set(MDIR ${MESHES}/upright)
  if(SIMMODSUITE_SimAdvMeshing_FOUND)
    mpi_test(parallel_meshgen 4
      ./generate
      "${MDIR}/upright.smd"
      "67k")
    mpi_test(parallel_meshgen_surf 4
      ./generate
      "--disable-volume"
      "--surface-mesh=${MDIR}/67k_surf.sms"
      "${MDIR}/upright.smd"
      "67k")
    mpi_test(parallel_meshgen_vol 4
      ./generate
      "--disable-surface"
      "--surface-mesh=${MDIR}/67k_surf_ref.sms"
      "${MDIR}/upright.smd"
      "67k")
    if(SIM_PARASOLID)
      mpi_test(parallel_meshgen_para 4
      ./generate
      "--native-model=${MDIR}/upright.x_t"
      "${MDIR}/upright.smd"
      "67k")
    endif()
    if(ENABLE_ZOLTAN)
      # adapt_meshgen uses the output of parallel_meshgen
      mpi_test(adapt_meshgen 4
        ./ma_test
        "${MDIR}/upright.smd"
        "67k/")
      set_test_depends(TESTS adapt_meshgen DEPENDS parallel_meshgen)
    endif()
  endif()
  if(SIM_PARASOLID)
    mpi_test(convert_para 1
      ./convert
      "--enable-log"
      "--native-model=${MDIR}/upright.x_t"
      "${MDIR}/upright.smd"
      "${MDIR}/67k.sms"
      "67k.smb")
    mpi_test(ph_convert_para 1
      ./ph_convert
      "--enable-log"
      "--attach-order"
      "--native-model=${MDIR}/upright.x_t"
      "${MDIR}/upright.smd"
      "${MDIR}/67k.sms"
      "67k.smb")
    set(MDIR ${MESHES}/curved)
    mpi_test(curvedSphere 1
      ./curvetest
      "${MDIR}/sphere1_geomSim.smd"
      "${MDIR}/sphere1_4.smb")
    mpi_test(highOrderSolutionTransfer 1
      ./highOrderSolutionTransfer
      "${MDIR}/sphere1_geomSim.smd"
      "${MDIR}/sphere1_4.smb")
    mpi_test(curvedKova 1
      ./curvetest
      "${MDIR}/Kova_geomSim.smd"
      "${MDIR}/Kova.smb")
    mpi_test(degen_shpere_full 1
      ./degenerate_test
      "${MDIR}/sph_full_geomSim.smd"
      "${MDIR}/sph_full.smb"
      "sph_full_refine"
      "3")
    mpi_test(degen_shpere_no_north 1
      ./degenerate_test
      "${MDIR}/sph_no_north_geomSim.smd"
      "${MDIR}/sph_no_north.smb"
      "sph_no_north_refine"
      "3")
    mpi_test(degen_shpere_vertical_slice 1
      ./degenerate_test
      "${MDIR}/sph_vertical_slice_geomSim.smd"
      "${MDIR}/sph_vertical_slice.smb"
      "sph_vertical_slice_refine"
      "3")
    mpi_test(crack_test 1
      ./crack_test
      "${MDIR}/crack_geomSim.smd")
  endif(SIM_PARASOLID)
endif()
if (PCU_COMPRESS)
  if(ENABLE_SIMMETRIX)
    set(RUNDIR run_sim)
  else()
    set(RUNDIR run)
  endif()
  set(MDIR ${MESHES}/phasta/1-1-Chef-Tet-Part/${RUNDIR})
  if(NOT APPLE)
    mpi_test(chefStream 1 ${CMAKE_CURRENT_BINARY_DIR}/chefStream
      WORKING_DIRECTORY ${MDIR})
  endif()
  mpi_test(chef0 1 ${CMAKE_CURRENT_BINARY_DIR}/chef
    WORKING_DIRECTORY ${MDIR})
  set(MDIR ${MESHES}/phasta/1-1-Chef-Tet-Part)
  if(ENABLE_SIMMETRIX)
    add_test(NAME chef1
      COMMAND diff -r ${RUNDIR}/1-procs_case/ good_phasta/
      WORKING_DIRECTORY ${MDIR})
    set_test_depends(TESTS chef1 DEPENDS chef0)
  endif()
  add_test(NAME chef2
    COMMAND diff -r out_mesh/ good_mesh/
    WORKING_DIRECTORY ${MDIR})
  if(ENABLE_ZOLTAN)
    mpi_test(chef3 2 ${CMAKE_CURRENT_BINARY_DIR}/chef
      WORKING_DIRECTORY ${MESHES}/phasta/2-1-Chef-Tet-Part/${RUNDIR})
    mpi_test(chef4 4 ${CMAKE_CURRENT_BINARY_DIR}/chef
      WORKING_DIRECTORY ${MESHES}/phasta/2-1-Chef-Tet-Part/4-2-Chef-Part/${RUNDIR})
    mpi_test(chef5 4 ${CMAKE_CURRENT_BINARY_DIR}/chef
      WORKING_DIRECTORY ${MESHES}/phasta/4-1-Chef-Tet-Part/${RUNDIR})
  endif()
  set(MDIR ${MESHES}/phasta/4-1-Chef-Tet-Part/4-4-Chef-Part-ts20)
  mpi_test(chef6 4 ${CMAKE_CURRENT_BINARY_DIR}/chef
    WORKING_DIRECTORY ${MDIR}/${RUNDIR})
  if(ENABLE_SIMMETRIX)
    add_test(NAME chef7
      COMMAND diff -r ${RUNDIR}/4-procs_case/ good_phasta/
      WORKING_DIRECTORY ${MDIR})
    set_test_depends(TESTS chef7 DEPENDS chef6)
  endif()
  add_test(NAME chef8
    COMMAND diff -r out_mesh/ good_mesh/
    WORKING_DIRECTORY ${MDIR})
  set_test_depends(TESTS chef8 DEPENDS chef6)
  if(ENABLE_ZOLTAN AND ENABLE_SIMMETRIX)
    mpi_test(chef9 2 ${CMAKE_CURRENT_BINARY_DIR}/chef
      WORKING_DIRECTORY ${MESHES}/phasta/simModelAndAttributes)
  endif()
  mpi_test(chefReadUrPrep 4 ${CMAKE_CURRENT_BINARY_DIR}/chefReadUrPrep
    ../../../model.dmg bz2:../good_mesh/ adapt.ur.inp
    WORKING_DIRECTORY ${MESHES}/phasta/4-1-Chef-Tet-Part/4-4-Chef-Part-ts20/run)
  if(ENABLE_ZOLTAN)
    mpi_test(chefReadRibUrPrep 4 ${CMAKE_CURRENT_BINARY_DIR}/chefReadUrPrep
      ../../../model.dmg bz2:../good_mesh/ adapt.prerib.inp
      WORKING_DIRECTORY ${MESHES}/phasta/4-1-Chef-Tet-Part/4-4-Chef-Part-ts20/run)
  endif()
  set(MDIR ${MESHES}/phasta/4-1-Chef-Tet-Part/4-4-Chef-Part-ts20)
  mpi_test(chef10 4 ${CMAKE_CURRENT_BINARY_DIR}/chef adapt.reducePartCount.inp
    WORKING_DIRECTORY ${MDIR}/run)
endif()

if(PUMI_ENABLE_CAPSTONE)
  set(MDIR "${MESHES}/cap")
  mpi_test(capVol_Wing 1 ./capVol 2 ${MDIR}/wing_surf_only.cre capVol_Wing.cre
  )
  mpi_test(capVol_WingMds 1
    ./capVol -m 2 ${MDIR}/wing_surf_only.cre capVol_WingMds.cre
  )
  mpi_test(capVol_BLCylMds3D 1
    ./capVol -gm 4 ${MDIR}/cyl_surf_only.cre capVol_BLCylMds3D.cre
  )
  mpi_test(capVol_BLCylMds3D-4 4
    ./capVol -gm 4 ${MDIR}/cyl_surf_only.cre capVol_BLCylMds3D.cre
  )
  mpi_test(capVol_CubeMds3D 1
    ./capVol -gm 3 ${MDIR}/cube_surf_only.cre capVol_CubeMds3D.cre
  )
  mpi_test(capVol_CubeMds3D-4 4
    ./capVol -gm 3 ${MDIR}/cube_surf_only.cre capVol_CubeMds3D.cre
  )
  if(PUMI_TEST_CAPVOL_EXTRA)
    mpi_test(capVol_BLCyl3D 1
      ./capVol -g 4 ${MDIR}/cyl_surf_only.cre capVol_BLCyl3D.cre
    )
    mpi_test(capVol_Cube3D 1
      ./capVol -g 3 ${MDIR}/cube_surf_only.cre capVol_Cube3D.cre
    )
    mpi_test(capVol_Cyl3D 1 ./capVol -g 1 ${MDIR}/cyl_surf_only.cre out.cre)
    mpi_test(capVol_Cube3D 1 ./capVol -g 3 ${MDIR}/cube_surf_only.cre out.cre)
    mpi_test(capVol_BLCyl3D 1 ./capVol -g 4 ${MDIR}/cyl_surf_only.cre out.cre)
    mpi_test(capVol_CylMds 1 ./capVol -m 1 ${MDIR}/cyl_surf_only.cre out.cre)
    mpi_test(capVol_CubeMds 1 ./capVol -m 3 ${MDIR}/cube_surf_only.cre out.cre)
    mpi_test(capVol_BLCylMds 1 ./capVol -m 4 ${MDIR}/cyl_surf_only.cre out.cre)
    mpi_test(capVol_CylMds3D 1 ./capVol -gm 1 ${MDIR}/cyl_surf_only.cre out.cre)
    mpi_test(capVol_Wing3D 1 ./capVol -g 2 ${MDIR}/wing_surf_only.cre out.cre)
    mpi_test(capVol_WingMds3D 1 ./capVol -gm 2 ${MDIR}/wing_surf_only.cre out.cre)
  endif()
  mpi_test(cap_inClosureOf 1 ./cap_inClosureOf ${MESHES}/cap/cube_surf_only.cre)
  mpi_test(cap_closestPoint 1 ./cap_closestPoint ${MESHES}/cap/cube_surf_only.cre)
  if(PUMI_HAS_CAPSTONE_SIZINGMETRICTOOL)
    mpi_test(cap_smooth 1 ./cap_smooth)
  endif()
endif()
