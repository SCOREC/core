set(MESHES "/lore/dibanez/meshes"
    CACHE string 
    "path to the meshes svn repo")
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
add_test(split_2
  split
  "${MDIR}/pipe.dmg"
  "pipe.smb"
  ${MESHFILE}
  2)
add_test(split_4
  ${MPIRUN} ${MPIRUN_PROCFLAG} 2
  ./zsplit
  "${MDIR}/pipe.dmg"
  ${MESHFILE}
  "pipe_4_.smb"
  2)
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
set(MDIR ${MESHES}/airFoilAfosr)
add_test(vtxBalance
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  ./vtxBalance
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
add_test(spr
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  ./spr_test
  "${MDIR}/spr.dmg"
  "${MDIR}/quadspr.smb")
set(MDIR ${MESHES}/nonmanifold)
add_test(nonmanif_verify
  ./verify
  "${MDIR}/nonmanifold.dmg"
  "${MDIR}/nonmanifold.smb")
add_test(nonmanif_split
  ./split
  "${MDIR}/nonmanifold.dmg"
  "${MDIR}/nonmanifold.smb"
  "nonmanifold_2_.smb"
  2)
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
  add_test(split_mpas
    split
    "mpas.dmg"
    "mpas.smb"
    "mpas_4_.smb"
    4)
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
set(MDIR ${MESHES}/fusion)
add_test(mkmodel_fusion
  mkmodel
  "${MDIR}/fusion.smb"
  "fusion.dmg")
add_test(split_fusion
  split
  "fusion.dmg"
  "${MDIR}/fusion.smb"
  "fusion_2_.smb"
  2)
# the part count mismatch is intentional,
# this test runs on half its procs
add_test(adapt_fusion
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  ./fusion
  "fusion_2_.smb")
add_test(change_dim
  ./newdim)
add_test(ma_insphere
  ma_insphere)
add_test(shapefun shapefun)
#todo - ph_test on crossflow ?
