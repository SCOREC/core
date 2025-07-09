#include <PCU.h>
#include "mylibrary.h"
int main(int argc, char** argv) {
  pcu::Init(&argc, &argv);
  pcu::PCU PCUObj;
  makeMesh(&PCUObj);
  pcu::Finalize();
  return 0;
}
