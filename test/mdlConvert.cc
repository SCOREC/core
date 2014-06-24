#include <stdio.h>
#include <SimUtil.h>
#include <SimModel.h>
#include "gmi_mesh.h"
#include "gmi_sim.h"

int main(int argc, char** argv)
{
  if (argc != 3) {
    printf("Convert simmetrix smd model to a gmi dmg model\n");
    printf("Usage: %s <simmetrix smd model> <gmi dmg model>\n", argv[0]);
    return 0;
  }
  Sim_readLicenseFile(0);
  SimModel_start();
  gmi_register_mesh();
  gmi_register_sim();
  gmi_model* m = gmi_load(argv[1]);
  gmi_write_dmg(m, argv[2]);
  gmi_destroy(m);
  SimModel_stop();
  Sim_unregisterAllKeys();
  return 0;
}
