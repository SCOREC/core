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
add_executable(test_coll test_coll.cc)
target_link_libraries(test_coll PRIVATE SCOREC::core)

enable_testing()
add_test(NAME coll1 COMMAND mpirun -np 1 $<TARGET_FILE:test_coll>)
add_test(NAME coll2 COMMAND mpirun -np 2 $<TARGET_FILE:test_coll>)
add_test(NAME coll4 COMMAND mpirun -np 4 $<TARGET_FILE:test_coll>)

