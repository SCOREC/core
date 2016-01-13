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
  macro(parma TESTNAME MDL IN OUT FACTOR METHOD APPROACH ISLOCAL PARTS)
    add_test("${TESTNAME}"
      ${MPIRUN} ${MPIRUN_PROCFLAG} ${PARTS} "./ptnParma" 
      ${MDL} ${IN} ${OUT} ${FACTOR} ${METHOD} ${APPROACH} ${ISLOCAL}
    )
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
  macro(parma TESTNAME MDL IN OUT FACTOR METHOD APPROACH ISLOCAL PARTS)
    math(EXPR OUTPARTS "${PARTS} * ${FACTOR}")
    add_test("${TESTNAME}"
      ${MPIRUN} ${MPIRUN_PROCFLAG} ${OUTPARTS} "./ptnParma_nothread" 
      ${MDL} ${IN} ${OUT} ${FACTOR} ${METHOD} ${APPROACH} ${ISLOCAL}
    )
  endmacro()
endif()
add_test(shapefun shapefun)
add_test(shapefun2 shapefun2)
add_test(bezierElevation bezierElevation)
add_test(bezierExperimental bezierExperimental)
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
if(ENABLE_ZOLTAN)
  splitfun(split_4
    ./zsplit
    "${MDIR}/pipe.dmg"
    ${MESHFILE}
    "pipe_4_.smb"
    2 2)
else()
  splitfun(split_4
    ./split
    "${MDIR}/pipe.dmg"
    ${MESHFILE}
    "pipe_4_.smb"
    2 2)
endif()
add_test(verify_parallel
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  ./verify
  "${MDIR}/pipe.dmg"
  "pipe_4_.smb")
add_test(vtxElmMixedBalance
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  ./vtxElmMixedBalance
  "${MDIR}/pipe.dmg"
  "pipe_4_.smb")
if(ENABLE_ZOLTAN)
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
endif()
set(MDIR ${MESHES}/torus)
add_test(balance
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  ./balance
  "${MDIR}/torus.dmg"
  "${MDIR}/4imb/torus.smb"
  "torusBal4p/")
if(ENABLE_ZOLTAN)
  add_test(zbalance
    ${MPIRUN} ${MPIRUN_PROCFLAG} 4
    ./zbalance
    "${MDIR}/torus.dmg"
    "${MDIR}/4imb/torus.smb"
    "torusZbal4p/")
endif()
add_test(gap
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  ./gap
  "${MDIR}/torus.dmg"
  "${MDIR}/4imb/torus.smb"
  "torusOpt4p/")
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
add_test(vtxEdgeElmBalance
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  ./vtxEdgeElmBalance
  "${MDIR}/afosr.dmg"
  "${MDIR}/4imb/"
  "afosrBal4p/"
  "2"
  "1.10")
add_test(vtxElmBalance
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  ./vtxElmBalance
  "${MDIR}/afosr.dmg"
  "${MDIR}/4imb/"
  "afosrBal4p/")
add_test(parmaSerial
  ${MPIRUN} ${MPIRUN_PROCFLAG} 1
  ./vtxElmBalance
  "${MESHES}/cube/cube.dmg"
  "${MESHES}/cube/pumi670/cube.smb"
  "cubeBal.smb/")
set(MDIR ${MESHES}/cube)
if(ENABLE_ZOLTAN)
  parma(ptnParma_cube
    "${MDIR}/cube.dmg"
    "${MDIR}/pumi670/cube.smb"
    "ptnParmaCube/" 
    "4" "rib" "reptn" "1"
    1)
endif()
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
if(ENABLE_ZOLTAN)
  add_test(adapt_fusion
    ${MPIRUN} ${MPIRUN_PROCFLAG} 4
    ./fusion
    "fusion_2_.smb")
endif()
add_test(fusion_field
  ${MPIRUN} ${MPIRUN_PROCFLAG} 2
  ./fusion2)
add_test(change_dim
  ./newdim)
add_test(ma_insphere
  ma_insphere)
if (PCU_COMPRESS)
  set(MDIR ${MESHES}/phasta/1-1-Chef-Tet-Part/run)
  if (PHASTA_CHEF_ENABLED)
    cook(chefStream ${CMAKE_CURRENT_BINARY_DIR}/chefStream 1 1 ${MDIR})
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
  cook(chef0 ${CMAKE_CURRENT_BINARY_DIR}/chef 1 1 ${MDIR})
  set(MDIR ${MESHES}/phasta/1-1-Chef-Tet-Part)
  add_test(NAME chef1
    COMMAND diff -r -x .svn out_mesh/ good_mesh/
    WORKING_DIRECTORY ${MDIR})
  set(MDIR ${MESHES}/phasta/2-1-Chef-Tet-Part/run)
  if(ENABLE_ZOLTAN)
    cook(chef2 ${CMAKE_CURRENT_BINARY_DIR}/chef 1 2 ${MDIR})
    set(MDIR ${MESHES}/phasta/2-1-Chef-Tet-Part/4-2-Chef-Part/run)
    cook(chef3 ${CMAKE_CURRENT_BINARY_DIR}/chef 2 2 ${MDIR})
    set(MDIR ${MESHES}/phasta/4-1-Chef-Tet-Part/run)
    cook(chef4 ${CMAKE_CURRENT_BINARY_DIR}/chef 1 4 ${MDIR})
  endif()
  set(MDIR ${MESHES}/phasta/4-1-Chef-Tet-Part/4-4-Chef-Part-ts20/run)
  cook(chef5 ${CMAKE_CURRENT_BINARY_DIR}/chef 4 1 ${MDIR})
  set(MDIR ${MESHES}/phasta/4-1-Chef-Tet-Part/4-4-Chef-Part-ts20)
  add_test(NAME chef6
    COMMAND diff -r -x .svn out_mesh/ good_mesh/
    WORKING_DIRECTORY ${MDIR})
  set(MDIR ${MESHES}/phasta/4-1-Chef-Tet-Part/4-4-Chef-Part-ts20/run)
  add_test(NAME chefReadUrPrep
    COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} 4
    ${CMAKE_CURRENT_BINARY_DIR}/chefReadUrPrep ../../../model.dmg bz2:../good_mesh/ adapt.inp
    WORKING_DIRECTORY "${MDIR}")
endif()
