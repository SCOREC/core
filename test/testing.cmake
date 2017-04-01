set(MESHES ""
    CACHE STRING
    "Extracted http://scorec.rpi.edu/pumi/pumi_test_meshes.tar.gz")
function(mpi_test TESTNAME PROCS EXE)
  add_test(
    NAME ${TESTNAME}
    COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} ${PROCS} ${VALGRIND} ${VALGRIND_ARGS} ${EXE} ${ARGN}
  )
endfunction(mpi_test)
mpi_test(shapefun 1 ./shapefun)
mpi_test(shapefun2 1 ./shapefun2)
mpi_test(bezierElevation 1 ./bezierElevation)
mpi_test(bezierMesh 1 ./bezierMesh)
mpi_test(bezierMisc 1 ./bezierMisc)
mpi_test(bezierRefine 1 ./bezierRefine)
mpi_test(bezierSubdivision 1 ./bezierSubdivision)
mpi_test(bezierValidity 1 ./bezierValidity)

mpi_test(align 1 ./align)
mpi_test(eigen_test 1 ./eigen_test)
mpi_test(integrate 1 ./integrate)
mpi_test(qr_test 1 ./qr)
mpi_test(base64 1 ./base64)
mpi_test(tensor_test 1 ./tensor)


if(ENABLE_SIMMETRIX)
  mpi_test(in_closure_of 1
    ./inClosureOf_test
    "${MESHES}/cube/cube.smd")
endif(ENABLE_SIMMETRIX)

if(ENABLE_SIMMETRIX)
  set(GXT smd)
else()
  set(GXT dmg)
endif()

set(MDIR ${MESHES}/phasta/dg)
if(ENABLE_SIMMETRIX)
  mpi_test(migrate_interface 4
    ./migrate_interface
    "${MDIR}/box.smd"
    "${MDIR}/box.smb"
    "${MDIR}/4/")
  mpi_test(dg_ma_test 4
    ./dg_ma_test
    "${MDIR}/box.smd"
    "${MDIR}/4/")
endif(ENABLE_SIMMETRIX)

set(MDIR ${MESHES}/phasta/loopDriver)
if(ENABLE_SIMMETRIX AND PCU_COMPRESS)
  mpi_test(ph_adapt 1
    ${CMAKE_CURRENT_BINARY_DIR}/ph_adapt
    "${MDIR}/model.smd"
    "${MDIR}/mesh_.smb"
    WORKING_DIRECTORY ${MDIR})
endif(ENABLE_SIMMETRIX AND PCU_COMPRESS)

if(ENABLE_ZOLTAN)
  mpi_test(pumi3d-1p 4
    ./test_pumi
    ${MESHES}/pumi/3d-1p/model.dmg
    ${MESHES}/pumi/3d-1p/part.smb
    out.smb 1 0)
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

set(MDIR ${MESHES}/fun3d)
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
set(MDIR ${MESHES}/pipe)
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
  mpi_test(convert_2d_quads 1
    ./convert
    "${MESHES}/disk/disk.smd"
    "${MESHES}/disk/disk_quad_mesh.sms"
    "disk_quad_mesh.smb")
else()
  file(COPY "${MESHES}/disk/disk_quad_mesh0.smb" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()
if(ENABLE_SIMMETRIX)
  mpi_test(convert_2d_tris 1
    ./convert
    "${MESHES}/disk/disk.smd"
    "${MESHES}/disk/disk_tri_mesh.sms"
    "disk_tri_mesh.smb")
else()
  file(COPY "${MESHES}/disk/disk_tri_mesh0.smb" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()
mpi_test(verify_2nd_order_shape_quads 1
  ./verify_2nd_order_shapes
  "disk_quad_mesh.smb")
mpi_test(verify_2nd_order_shape_tris 1
  ./verify_2nd_order_shapes
  "disk_tri_mesh.smb")
mpi_test(uniform_serial 1
  ./uniform
  "${MDIR}/pipe.${GXT}"
  "pipe.smb"
  "pipe_unif.smb")
if(ENABLE_SIMMETRIX)
  mpi_test(snap_serial 1
    ./snap
    "${MDIR}/pipe.${GXT}"
    "pipe_unif.smb"
    "pipe.smb")
endif()
mpi_test(ma_serial 1
  ./ma_test
  "${MDIR}/pipe.${GXT}"
  "pipe.smb")
mpi_test(aniso_ma_serial 1
  ./aniso_ma_test
  "${MESHES}/cube/cube.dmg"
  "${MESHES}/cube/pumi670/cube.smb")
if(ENABLE_ZOLTAN)
  mpi_test(torus_ma_parallel 4
    ./torus_ma_test
    "${MESHES}/torus/torus.dmg"
    "${MESHES}/torus/4imb/torus.smb")
endif()
mpi_test(tet_serial 1
  ./tetrahedronize
  "${MDIR}/pipe.${GXT}"
  "pipe.smb"
  "tet.smb")
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
if(ENABLE_ZOLTAN)
  mpi_test(refineX 2
    ./refine2x
    "${MDIR}/pipe.${GXT}"
    ${MESHFILE}
    0
    "refXpipe/")
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
endif()
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
if(ENABLE_ZOLTAN)
  mpi_test(zbalance 4
    ./zbalance
    "${MDIR}/torus.dmg"
    "${MDIR}/4imb/torus.smb"
    "torusZbal4p/")
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
if(ENABLE_ZOLTAN AND ENABLE_SIMMETRIX)
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
mpi_test(construct 4
  ./construct
  "${MDIR}/cube.dmg"
  "${MDIR}/pumi7k/4/cube.smb")
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
add_test(hierarchic_2p_2D
  ./hierarchic
  "${MDIR}/square.dmg"
  "${MDIR}/square.smb"
  2)
add_test(hierarchic_3p_2D
  ./hierarchic
  "${MDIR}/square.dmg"
  "${MDIR}/square.smb"
  3)
set(MDIR ${MESHES}/cube/pumi24)
add_test(hierarchic_2p_3D
  ./hierarchic
  "${MDIR}/cube.dmg"
  "${MDIR}/cube.smb"
  2)
set(MDIR ${MESHES}/cube)
mpi_test(test_verify 4
  ./test_verify
  "${MDIR}/cube.dmg"
  "${MDIR}/pumi7k/4/cube.smb")
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
set(MDIR ${MESHES}/fusion)
mpi_test(mkmodel_fusion 1
  ./mkmodel
  "${MDIR}/fusion.smb"
  "fusion.dmg")
mpi_test(split_fusion 2
  ./split
  "fusion.dmg"
  "${MDIR}/fusion.smb"
  "fusion_2_.smb"
  2)
# the part count mismatch is intentional,
# this test runs on half its procs
if(ENABLE_ZOLTAN)
  mpi_test(adapt_fusion 4
    ./fusion
    "fusion_2_.smb")
endif()
mpi_test(fusion_field 2
  ./fusion2)
mpi_test(change_dim 1
  ./newdim)
mpi_test(ma_insphere 1
  ./ma_insphere)
if(ENABLE_SIMMETRIX)
  set(MDIR ${MESHES}/upright)
  mpi_test(parallel_meshgen 4
    ./generate
    "${MDIR}/upright.smd"
    "67k")
  mpi_test(adapt_meshgen 4
    ./ma_test
    "${MDIR}/upright.smd"
    "67k/")
  if(SIM_PARASOLID)
    set(MDIR ${MESHES}/curved)
    mpi_test(curvedSphere 1
      ./curvetest
      "${MDIR}/sphere1.xmt_txt"
      "${MDIR}/sphere1_4.smb")
    mpi_test(curvedKova 1
      ./curvetest
      "${MDIR}/Kova.xmt_txt"
      "${MDIR}/Kova.smb")
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
  endif()
  add_test(NAME chef8
    COMMAND diff -r out_mesh/ good_mesh/
    WORKING_DIRECTORY ${MDIR})
  if(ENABLE_SIMMETRIX)
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
endif()
