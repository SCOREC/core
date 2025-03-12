function(smoke_test TESTNAME PROCS EXE)
  set(tname smoke_test_${TESTNAME})
  if(SCOREC_NO_MPI)
    if(PROCS EQUAL "1")
      add_test(NAME ${tname} COMMAND ${VALGRIND} ${VALGRIND_ARGS}
        ${EXE} ${ARGN})
    else()
      return()
    endif()
  else()
    add_test(
      NAME ${tname}
      COMMAND ${MPIRUN} ${MPIRUN_PROCFLAG} ${PROCS}
        ${VALGRIND} ${VALGRIND_ARGS} ${EXE} ${ARGN})
  endif()
  SET_TESTS_PROPERTIES(${tname} PROPERTIES LABELS "SMOKE_TEST" )
endfunction(smoke_test)

set(MDIR ${SMOKE_TEST_MESHES}/pipe)
smoke_test(uniform_serial 1
  ./uniform
  "${MDIR}/pipe.dmg"
  "${MDIR}/pipe.smb"
  "pipe_unif.smb")

smoke_test(split_2 2
  ./split
  "${MDIR}/pipe.dmg"
  "${MDIR}/pipe.smb"
  "pipe_2_.smb"
  2)

include(GNUInstallDirs)
# install the test input files for use in spack post-install tests
install(FILES "${MDIR}/pipe.dmg" "${MDIR}/pipe0.smb"
  DESTINATION ${CMAKE_INSTALL_DATADIR}/testdata)

