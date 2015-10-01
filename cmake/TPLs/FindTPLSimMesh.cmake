INCLUDE(TribitsTplDeclareLibraries)

TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( SimMesh
  REQUIRED_HEADERS MeshSim.h SimPartitionedMesh.h
  REQUIRED_LIBS_NAMES SimPartitionedMesh-mpi${SIM_DEBUG} SimMeshing${SIM_DEBUG} SimMeshTools${SIM_DEBUG} SimPartitionWrapper-${SIM_MPI}
)

