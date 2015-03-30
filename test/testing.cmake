set(MESHES "/lore/dibanez/meshes"
    CACHE string
    "path to the meshes svn repo")
if(ENABLE_THREADS)
  macro(splitfun TESTNAME PROG MODEL IN OUT PARTS FACTOR)
    add_test("${TESTNAME}"
      ${MPIRUN} ${MPIRUN_PROCFLAG} ${PARTS}
      "${PROG}"
      "${MODEL}"
      "${IN}"
      "${OUT}"
      ${FACTOR})
  endmacro()
  macro(cook TESTNAME PROG PARTS FACTOR WORKDIR)
    add_test(NAME "${TESTNAME}"
      COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} ${PARTS} "${PROG}"
      WORKING_DIRECTORY "${WORKDIR}")
  endmacro()
else()
  macro(splitfun TESTNAME PROG MODEL IN OUT PARTS FACTOR)
    math(EXPR OUTPARTS "${PARTS} * ${FACTOR}")
    add_test("${TESTNAME}"
      ${MPIRUN} ${MPIRUN_PROCFLAG} ${OUTPARTS}
      "${PROG}"
      "${MODEL}"
      "${IN}"
      "${OUT}"
      ${FACTOR})
  endmacro()
  macro(cook TESTNAME PROG PARTS FACTOR WORKDIR)
    math(EXPR OUTPARTS "${PARTS} * ${FACTOR}")
    add_test(NAME "${TESTNAME}"
      COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} ${OUTPARTS} "${PROG}"
      WORKING_DIRECTORY "${WORKDIR}")
  endmacro()
endif()
add_test(shapefun shapefun)
add_test(eigen_test eigen_test)
add_test(integrate integrate)
add_test(qr_test qr_test)
set(MDIR ${MESHES}/pipe)
add_test(verify_serial
  verify
  "${MDIR}/pipe.dmg"
  "${MDIR}/pipe.smb")
add_test(uniform_serial
  uniform
  "${MDIR}/pipe.dmg"
  "${MDIR}/pipe.smb"
  "pipe.smb")
add_test(ma_serial
  ma_test
  "${MDIR}/pipe.dmg"
  "pipe.smb")
add_test(tet_serial
  tetrahedronize
  "${MDIR}/pipe.dmg"
  "pipe.smb"
  "tet.smb")
if (PCU_COMPRESS)
  set(MESHFILE "bz2:pipe_2_.smb")
else()
  set(MESHFILE "pipe_2_.smb")
endif()
splitfun(split_2
  ./split
  "${MDIR}/pipe.dmg"
  "pipe.smb"
  ${MESHFILE}
  1 2)
add_test(refineX
  ${MPIRUN} ${MPIRUN_PROCFLAG} 2
  ./refine2x
  "${MDIR}/pipe.dmg"
  ${MESHFILE}
  0
  "refXpipe/")
splitfun(split_4
  ./zsplit
  "${MDIR}/pipe.dmg"
  ${MESHFILE}
  "pipe_4_.smb"
  2 2)
add_test(verify_parallel
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  ./verify
  "${MDIR}/pipe.dmg"
  "pipe_4_.smb")
add_test(ma_parallel
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  ./ma_test
  "${MDIR}/pipe.dmg"
  "pipe_4_.smb")
add_test(tet_parallel
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  ./tetrahedronize
  "${MDIR}/pipe.dmg"
  "pipe_4_.smb"
  "tet.smb")
set(MDIR ${MESHES}/torus)
add_test(balance
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  ./balance
  "${MDIR}/torus.dmg"
  "${MDIR}/4imb/torus.smb"
  "torusBal4p/")
add_test(zbalance
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  ./zbalance
  "${MDIR}/torus.dmg"
  "${MDIR}/4imb/torus.smb"
  "torusZbal4p/")
add_test(gap
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  ./gap
  "${MDIR}/torus.dmg"
  "${MDIR}/4imb/torus.smb"
  "torusOpt4p/")
add_test(hps
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  ./hps
  "${MDIR}/torus.dmg"
  "${MDIR}/4imb/torus.smb"
  "torusHps4p/")
add_test(fixDisconnected
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  ./fixDisconnected
  "${MDIR}/torus.dmg"
  "${MDIR}/4imb/torus.smb"
  "torusDcFix4p/")
set(MDIR ${MESHES}/airFoilAfosr)
add_test(elmBalance
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  ./elmBalance
  "${MDIR}/afosr.dmg"
  "${MDIR}/4imb/"
  "afosrBal4p/")
add_test(vtxBalance
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  ./vtxBalance
  "${MDIR}/afosr.dmg"
  "${MDIR}/4imb/"
  "afosrBal4p/")
add_test(edgeBalance
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  ./edgeBalance
  "${MDIR}/afosr.dmg"
  "${MDIR}/4imb/"
  "afosrBal4p/")
add_test(vtxEdgeElmBalance
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  ./vtxEdgeElmBalance
  "${MDIR}/afosr.dmg"
  "${MDIR}/4imb/"
  "afosrBal4p/")
add_test(vtxElmBalance
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  ./vtxElmBalance
  "${MDIR}/afosr.dmg"
  "${MDIR}/4imb/"
  "afosrBal4p/")
set(MDIR ${MESHES}/cube)
add_test(construct
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  ./construct
  "${MDIR}/cube.dmg"
  "${MDIR}/pumi7k/4/cube.smb")
set(MDIR ${MESHES}/spr)
add_test(spr_3D
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  ./spr_test
  "${MDIR}/spr.dmg"
  "${MDIR}/quadspr.smb"
  spr3D
  2)
add_test(spr_2D
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  ./spr_test
  "${MDIR}/square.dmg"
  "${MDIR}/square.smb"
  spr2D
  1)
set(MDIR ${MESHES}/nonmanifold)
add_test(nonmanif_verify
  ./verify
  "${MDIR}/nonmanifold.dmg"
  "${MDIR}/nonmanifold.smb")
splitfun(nonmanif_split
  ./split
  "${MDIR}/nonmanifold.dmg"
  "${MDIR}/nonmanifold.smb"
  "nonmanifold_2_.smb"
  1 2)
add_test(nonmanif_verify2
  ${MPIRUN} ${MPIRUN_PROCFLAG} 2
  ./verify
  "${MDIR}/nonmanifold.dmg"
  "nonmanifold_2_.smb")
if (ENABLE_MPAS)
  set(MDIR ${MESHES}/mpas)
  add_test(read_mpas
    mpas_read
    "${MDIR}/ocean_QU_240km.nc"
    "mpas.dmg"
    "mpas.smb")
  if(ENABLE_THREADS)
    splitfun(split_mpas
      ./split
      "mpas.dmg"
      "mpas.smb"
      "mpas_4_.smb"
      1 4)
    add_test(verify_mpas
      ${MPIRUN} ${MPIRUN_PROCFLAG} 4
      ./verify
      "mpas.dmg"
      "mpas_4_.smb")
    add_test(ghost_mpas
      ${MPIRUN} ${MPIRUN_PROCFLAG} 4
      ./ghost
      "mpas.dmg"
      "mpas_4_.smb"
      "ghost_4_.smb")
    add_test(write_mpas
      ${MPIRUN} ${MPIRUN_PROCFLAG} 4
      ./mpas_write
      "mpas.dmg"
      "ghost_4_.smb"
      "${MDIR}/ocean_QU_240km.nc"
      "mpas_part_")
  endif()
endif()
set(MDIR ${MESHES}/fusion)
add_test(mkmodel_fusion
  mkmodel
  "${MDIR}/fusion.smb"
  "fusion.dmg")
splitfun(split_fusion
  ./split
  "fusion.dmg"
  "${MDIR}/fusion.smb"
  "fusion_2_.smb"
  1 2)
# the part count mismatch is intentional,
# this test runs on half its procs
add_test(adapt_fusion
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  ./fusion
  "fusion_2_.smb")
add_test(fusion_field
  ${MPIRUN} ${MPIRUN_PROCFLAG} 2
  ./fusion2)
add_test(change_dim
  ./newdim)
add_test(ma_insphere
  ma_insphere)
if (PCU_COMPRESS)
  set(MDIR ${MESHES}/phasta/1-1-Chef-Tet-Part/run)
  cook(chef0 ${CMAKE_CURRENT_BINARY_DIR}/chef 1 1 ${MDIR})
  set(MDIR ${MESHES}/phasta/1-1-Chef-Tet-Part)
  add_test(NAME chef1
    COMMAND diff -r -x .svn out_mesh/ good_mesh/
    WORKING_DIRECTORY ${MDIR})
  set(MDIR ${MESHES}/phasta/2-1-Chef-Tet-Part/run)
  cook(chef2 ${CMAKE_CURRENT_BINARY_DIR}/chef 1 2 ${MDIR})
  set(MDIR ${MESHES}/phasta/2-1-Chef-Tet-Part/4-2-Chef-Part/run)
  cook(chef3 ${CMAKE_CURRENT_BINARY_DIR}/chef 2 2 ${MDIR})
  set(MDIR ${MESHES}/phasta/4-1-Chef-Tet-Part/run)
  cook(chef4 ${CMAKE_CURRENT_BINARY_DIR}/chef 1 4 ${MDIR})
  set(MDIR ${MESHES}/phasta/4-1-Chef-Tet-Part/4-4-Chef-Part-ts20/run)
  cook(chef5 ${CMAKE_CURRENT_BINARY_DIR}/chef 4 1 ${MDIR})
  set(MDIR ${MESHES}/phasta/4-1-Chef-Tet-Part/4-4-Chef-Part-ts20)
  add_test(NAME chef6
    COMMAND diff -r -x .svn out_mesh/ good_mesh/
    WORKING_DIRECTORY ${MDIR})
endif()
