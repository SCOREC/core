#include <gmi_mesh.h>
#include <assert.h>
#include <PCU.h>

int main(int argc, char** argv)
{
  struct gmi_model* m;
  assert(argc == 3);
  gmi_register_mesh();
  m = gmi_load(argv[1]);
  gmi_write_dmg(m, argv[2]);
  gmi_destroy(m);
  return 0;
}
