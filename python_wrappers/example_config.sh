cmake .. \
	-DSCOREC_PREFIX="/usr/local/core_install" \
	-DCMAKE_INSTALL_PREFIX="/usr/local/pcu_python" \
  -DCMAKE_C_COMPILER="mpicc" \
  -DCMAKE_CXX_COMPILER="mpicxx" \
  -DENABLE_PYTHON=ON
