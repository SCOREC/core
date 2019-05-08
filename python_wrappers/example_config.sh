cmake .. \
	-DSCOREC_PREFIX="/usr/local/core_install" \
	-DCMAKE_INSTALL_PREFIX="/usr/local/pyCore" \
  -DCMAKE_C_COMPILER="mpicc" \
  -DCMAKE_CXX_COMPILER="mpicxx" \
  -DENABLE_PYTHON=ON
