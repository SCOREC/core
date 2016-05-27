set(MESHES "/lore/dibanez/meshes"
    CACHE string
    "path to the meshes svn repo")
macro(mpi_test TESTNAME PROCS EXE)
  add_test(
    NAME ${TESTNAME}
    COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} ${PROCS} ${VALGRIND} ${VALGRIND_ARGS} ${EXE} ${ARGN}
  )
endmacro(mpi_test)
add_test(shapefun shapefun)
add_test(shapefun2 shapefun2)
add_test(bezierElevation bezierElevation)
add_test(bezierMesh bezierMesh)
add_test(bezierMisc bezierMisc)
add_test(bezierRefine bezierRefine)
add_test(bezierSubdivision bezierSubdivision)
add_test(bezierValidity bezierValidity)

add_test(align align)
add_test(eigen_test eigen_test)
add_test(integrate integrate)
add_test(qr_test qr)
add_test(base64 base64)
add_test(tensor_test tensor)

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

set(MDIR ${MESHES}/fun3d)
mpi_test(inviscid_ugrid 1
  ./from_ugrid
  "${MDIR}/inviscid_egg.b8.ugrid"
  "${MDIR}/inviscid_egg.dmg"
  "${MDIR}/inviscid_egg.smb")
mpi_test(inviscid_split 4
  ./split
  "${MDIR}/inviscid_egg.dmg"
  "${MDIR}/inviscid_egg.smb"
  "${MDIR}/4/"
  4)
mpi_test(inviscid_ghost 4
  ./ghost
  "${MDIR}/inviscid_egg.dmg"
  "${MDIR}/4/"
  "${MDIR}/vis")
set(MDIR ${MESHES}/pipe)
mpi_test(convert 1
  ./convert
  "${MDIR}/pipe.smd"
  "${MDIR}/pipe.sms"
  "pipe.smb")
mpi_test(verify_serial 1
  ./verify
  "${MDIR}/pipe.smd"
  "pipe.smb")
mpi_test(uniform_serial 1
  ./uniform
  "${MDIR}/pipe.smd"
  "pipe.smb"
  "pipe_unif.smb")
mpi_test(snap_serial 1
  ./snap
  "${MDIR}/pipe.smd"
  "pipe_unif.smb"
  "pipe.smb")
mpi_test(ma_serial 1
  ./ma_test
  "${MDIR}/pipe.smd"
  "pipe.smb")
mpi_test(tet_serial 1
  ./tetrahedronize
  "${MDIR}/pipe.smd"
  "pipe.smb"
  "tet.smb")
if (PCU_COMPRESS)
  set(MESHFILE "bz2:pipe_2_.smb")
else()
  set(MESHFILE "pipe_2_.smb")
endif()
mpi_test(split_2 2
  ./split
  "${MDIR}/pipe.smd"
  "pipe.smb"
  ${MESHFILE}
  2)
if(ENABLE_ZOLTAN)
  mpi_test(refineX 2
    ./refine2x
    "${MDIR}/pipe.dmg"
    ${MESHFILE}
    0
    "refXpipe/")
  mpi_test(split_4 4
    ./zsplit
    "${MDIR}/pipe.smd"
    ${MESHFILE}
    "pipe_4_.smb"
    2)
else()
  mpi_test(split_4 4
    ./split
    "${MDIR}/pipe.smd"
    ${MESHFILE}
    "pipe_4_.smb"
    2)
endif()
mpi_test(verify_parallel 4
  ./verify
  "${MDIR}/pipe.smd"
  "pipe_4_.smb")
mpi_test(vtxElmMixedBalance 4
  ./vtxElmMixedBalance
  "${MDIR}/pipe.dmg"
  "pipe_4_.smb")
if(ENABLE_ZOLTAN)
  mpi_test(ma_parallel 4
    ./ma_test
    "${MDIR}/pipe.smd"
    "pipe_4_.smb")
  mpi_test(tet_parallel 4
    ./tetrahedronize
    "${MDIR}/pipe.smd"
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
  "${MDIR}/afosr.smd"
  "${MDIR}/4imb/"
  "afosrBal4p/")
mpi_test(vtxEdgeElmBalance 4
  ./vtxEdgeElmBalance
  "${MDIR}/afosr.smd"
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
set(MDIR ${MESHES}/upright)
mpi_test(parallel_meshgen 4
  ./generate
  "${MDIR}/upright.smd"
  "67k")
mpi_test(adapt_meshgen 4
  ./ma_test
  "${MDIR}/upright.smd"
  "67k/")
mpi_test(ma_insphere 1
  ./ma_insphere)
set(MDIR ${MESHES}/curved)
mpi_test(curvedSphere 1
  ./curvetest
  "${MDIR}/sphere1.xmt_txt"
  "${MDIR}/sphere1_4.smb")
mpi_test(curvedKova 1
  ./curvetest
  "${MDIR}/Kova.xmt_txt"
  "${MDIR}/Kova.smb")
if (PCU_COMPRESS)
  set(MDIR ${MESHES}/phasta/1-1-Chef-Tet-Part/run_sim)
  if (PHASTA_CHEF_ENABLED)
    mpi_test(chefStream 1 ${CMAKE_CURRENT_BINARY_DIR}/chefStream
      WORKING_DIRECTORY ${MDIR})
    set(cmd 
      ${CMAKE_BINARY_DIR}/phasta/bin/checkphasta 
      ${MDIR}/1-procs_case/ 
      ${MESHES}/phasta/1-1-Chef-Tet-Part/good_phasta/
      0 1e-6)
    add_test(
      NAME compareChefStream
      COMMAND ${cmd}
      WORKING_DIRECTORY ${MDIR}
    )
  endif()
  mpi_test(chef0 1 ${CMAKE_CURRENT_BINARY_DIR}/chef
    WORKING_DIRECTORY ${MDIR})
  set(MDIR ${MESHES}/phasta/1-1-Chef-Tet-Part)
  add_test(NAME chef1
    COMMAND diff -r -x .svn run_sim/1-procs_case/ good_phasta/
    WORKING_DIRECTORY ${MDIR})
  add_test(NAME chef2
    COMMAND diff -r -x .svn out_mesh/ good_mesh/
    WORKING_DIRECTORY ${MDIR})
  if(ENABLE_ZOLTAN)
    mpi_test(chef3 2 ${CMAKE_CURRENT_BINARY_DIR}/chef
      WORKING_DIRECTORY ${MESHES}/phasta/2-1-Chef-Tet-Part/run_sim)
    mpi_test(chef4 4 ${CMAKE_CURRENT_BINARY_DIR}/chef
      WORKING_DIRECTORY ${MESHES}/phasta/2-1-Chef-Tet-Part/4-2-Chef-Part/run_sim)
    mpi_test(chef5 4 ${CMAKE_CURRENT_BINARY_DIR}/chef
      WORKING_DIRECTORY ${MESHES}/phasta/4-1-Chef-Tet-Part/run_sim)
  endif()
  mpi_test(chef6 4 ${CMAKE_CURRENT_BINARY_DIR}/chef
    WORKING_DIRECTORY ${MESHES}/phasta/4-1-Chef-Tet-Part/4-4-Chef-Part-ts20/run_sim)
  set(MDIR ${MESHES}/phasta/4-1-Chef-Tet-Part/4-4-Chef-Part-ts20)
  add_test(NAME chef7
    COMMAND diff -r -x .svn run_sim/4-procs_case/ good_phasta/
    WORKING_DIRECTORY ${MDIR})
  add_test(NAME chef8
    COMMAND diff -r -x .svn out_mesh/ good_mesh/
    WORKING_DIRECTORY ${MDIR})
  mpi_test(chef9 2 ${CMAKE_CURRENT_BINARY_DIR}/chef
    WORKING_DIRECTORY ${MESHES}/phasta/simModelAndAttributes)
  mpi_test(chefReadUrPrep 4 ${CMAKE_CURRENT_BINARY_DIR}/chefReadUrPrep
    ../../../model.dmg bz2:../good_mesh/ adapt.inp
    WORKING_DIRECTORY ${MESHES}/phasta/4-1-Chef-Tet-Part/4-4-Chef-Part-ts20/run)
  if(ENABLE_ZOLTAN)
    mpi_test(chefReadRibUrPrep 4 ${CMAKE_CURRENT_BINARY_DIR}/chefReadUrPrep
      ../../../model.dmg bz2:../good_mesh/ adapt.prerib.inp
      WORKING_DIRECTORY ${MESHES}/phasta/4-1-Chef-Tet-Part/4-4-Chef-Part-ts20/run)
  endif()
endif()
