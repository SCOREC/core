function(util_exe_func exename srcname)
  add_executable(${exename} ${srcname})
  target_link_libraries(${exename} core apf_cap)
  bob_export_target(${exename})
  bob_end_subdir()
endfunction(util_exe_func)

function(test_exe_func exename srcname)
  if(IS_TESTING)
    add_executable(${exename} ${srcname})
  else()
    add_executable(${exename} EXCLUDE_FROM_ALL ${srcname})
  endif()
  target_link_libraries(${exename} core)
endfunction(test_exe_func)

# Capstone utilities
if(ENABLE_CAPSTONE)
  util_exe_func(capStoneAttachSolution capStoneAttachSolution.cc)
endif()

bob_end_subdir()
