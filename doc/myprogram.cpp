#include <PCU.h>
#include "mylibrary.h"
int main(int argc, char** argv) {
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  makeMesh();
  PCU_Comm_Free();
  MPI_Finalize();
  return 0;
}
