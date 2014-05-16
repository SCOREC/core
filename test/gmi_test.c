#include "gmi_mesh.h"

int main(int argc, char** argv)
{
  gmi_register_mesh();
  gmi_destroy(gmi_load(argv[1]));
  return 0;
}
