#include <SimUtil.h>
#include <SimModel.h>
#include "gmi_mesh.h"
#include "gmi_sim.h"

int main(int argc, char** argv)
{
  Sim_readLicenseFile(0);
  SimModel_start();
  gmi_register_mesh();
  gmi_register_sim();
  gmi_destroy(gmi_load(argv[1]));
  SimModel_stop();
  Sim_unregisterAllKeys();
  return 0;
}
