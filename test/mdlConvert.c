#include <gmi_mesh.h>
#include <assert.h>
#include <PCU.h>

int main(int argc, char** argv)
{
  assert(argc == 3);
  struct gmi_model* m;
  gmi_register_mesh();
  m = gmi_load(argv[1]);
  gmi_write_dmg(m, argv[2]);
  gmi_destroy(m);
  return 0;
}
