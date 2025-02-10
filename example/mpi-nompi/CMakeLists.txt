cmake_minimum_required(VERSION 3.8)
project(mpi-nompi VERSION 1.0.0 LANGUAGES CXX)

# Allow the user to indicate where they installed SCOREC
# via "-DSCOREC_PREFIX=/home/somewhere" when calling `cmake`
set(SCOREC_PREFIX "" CACHE STRING "Directory where SCOREC is installed")

# If SCOREC_PREFIX was specified then use only that directory.
if (SCOREC_PREFIX)
  find_package(SCOREC 3 REQUIRED CONFIG PATHS ${SCOREC_PREFIX} NO_DEFAULT_PATH)
else()
  find_package(SCOREC 3 REQUIRED CONFIG)
endif()

add_executable(hello hello.cc)
target_link_libraries(hello PRIVATE SCOREC::core)

enable_testing()
add_test(NAME hello1 COMMAND mpirun -np 1 $<TARGET_FILE:hello>)
set_property(TEST hello1 PROPERTY PASS_REGULAR_EXPRESSION 
  "MPI size: 1\; PCU size: 1
Hello from MPI rank: 0\; PCU rank: 0
MPI\\(0\\) sum: 1\; PCU sum: 1"
)
add_test(NAME hello2 COMMAND mpirun -np 2 $<TARGET_FILE:hello>)
set_property(TEST hello2 PROPERTY PASS_REGULAR_EXPRESSION 
  "MPI size: 2\; PCU size: 1
Hello from MPI rank: 0\; PCU rank: 0
Hello from MPI rank: 1\; PCU rank: 0
MPI\\(0\\) sum: 3\; PCU sum: 1
MPI\\(1\\) sum: 3\; PCU sum: 2"
)
add_test(NAME hello4 COMMAND mpirun -np 4 $<TARGET_FILE:hello>)
set_property(TEST hello4 PROPERTY PASS_REGULAR_EXPRESSION 
  "MPI size: 4\; PCU size: 1
Hello from MPI rank: 0\; PCU rank: 0
Hello from MPI rank: 1\; PCU rank: 0
Hello from MPI rank: 2\; PCU rank: 0
Hello from MPI rank: 3\; PCU rank: 0
MPI\\(0\\) sum: 10\; PCU sum: 1
MPI\\(1\\) sum: 10\; PCU sum: 2
MPI\\(2\\) sum: 10\; PCU sum: 3
MPI\\(3\\) sum: 10\; PCU sum: 4"
)

