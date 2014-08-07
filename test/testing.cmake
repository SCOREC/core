set(MESHES /lore/dibanez/meshes)
set(MDIR ${MESHES}/pipe)
add_test(verify_serial
  verify
  "${MDIR}/pipe.smd"
  "${MDIR}/pipe.smb")
add_test(uniform_serial
  uniform
  "${MDIR}/pipe.smd"
  "${MDIR}/pipe.smb"
  "pipe.smb")
add_test(ma_serial
  ma_test
  "${MDIR}/pipe.smd"
  "pipe.smb")
add_test(tet_serial
  tetrahedronize
  "${MDIR}/pipe.smd"
  "pipe.smb"
  "tet.smb")
add_test(split_2
  split
  "${MDIR}/pipe.smd"
  "pipe.smb"
  "pipe_2_.smb"
  2)
add_test(split_4
  ${MPIRUN} ${MPIRUN_PROCFLAG} 2
  zsplit
  "${MDIR}/pipe.smd"
  "pipe_2_.smb"
  "pipe_4_.smb"
  2)
add_test(verify_parallel
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  ma_test
  "${MDIR}/pipe.smd"
  "pipe_4_.smb")
add_test(ma_parallel
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  ma_test
  "${MDIR}/pipe.smd"
  "pipe_4_.smb")
add_test(tet_parallel
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  tetrahedronize
  "${MDIR}/pipe.smd"
  "pipe_4_.smb"
  "tet.smb")
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
  verify
  "mpas.dmg"
  "mpas_4_.smb")
add_test(ghost_mpas
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  ghost
  "mpas.dmg"
  "mpas_4_.smb"
  "ghost_4_.smb")
add_test(write_mpas
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  mpas_write
  "mpas.dmg"
  "ghost_4_.smb"
  "${MDIR}/ocean_QU_240km.nc"
  "mpas_part_")
set(MDIR ${MESHES}/fusion)
add_test(mkmodel_fusion
  mkmodel
  "${MDIR}/fusion.smb"
  "fusion.dmg")
add_test(split_fusion
  split
  "fusion.dmg"
  "${MDIR}/fusion.smb"
  "fusion_4_.smb"
  4)
add_test(adapt_fusion
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  fan
  "fusion_4_.smb")
set(MDIR ${MESHES}/upright)
add_test(parallel_meshgen
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  generate
  "${MDIR}/upright.smd"
  "67k")
add_test(adapt_meshgen
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  ma_test
  "${MDIR}/upright.smd"
  "67k/")
#todo - ph_test on crossflow, fusion (fan.cc), etc...
