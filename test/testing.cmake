set(MDIR ../../meshes/pipe)
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
add_test(split_2
  split
  "${MDIR}/pipe.dmg"
  "pipe.smb"
  "pipe_2_.smb"
  2)
add_test(split_4
  ${MPIRUN} ${MPIRUN_PROCFLAG} 2
  zsplit
  "${MDIR}/pipe.dmg"
  "pipe_2_.smb"
  "pipe_4_.smb"
  2)
add_test(verify_parallel
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  ma_test
  "${MDIR}/cube.dmg"
  "pipe_4_.smb")
add_test(ma_parallel
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  ma_test
  "${MDIR}/cube.dmg"
  "pipe_4_.smb")
add_test(tet_parallel
  ${MPIRUN} ${MPIRUN_PROCFLAG} 4
  tetrahedronize
  "${MDIR}/pipe.dmg"
  "pipe_4_.smb"
  "tet.smb")
#todo - mpas, ph_test on crossflow, ...
