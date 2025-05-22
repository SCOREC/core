#include <PCU.h>
#include "mylibrary.h"
int main(int argc, char** argv) {
#ifndef SCOREC_NO_MPI
  MPI_Init(&argc,&argv);
#endif
  pcu::PCU PCUObj;
  makeMesh(&PCUObj);
#ifndef SCOREC_NO_MPI
  MPI_Finalize();
#endif
  return 0;
}
