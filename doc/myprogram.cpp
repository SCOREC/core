#include <PCU.h>
#include "mylibrary.h"
int main(int argc, char** argv) {
  MPI_Init(&argc,&argv);
  pcu::PCU PCUObj = pcu::PCU;
  makeMesh(&PCUObj);
  MPI_Finalize();
  return 0;
}
