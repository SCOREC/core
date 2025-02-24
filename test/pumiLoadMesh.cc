#include "pumi.h"

int main(int argc, char** argv)
{
  pcu::PCU_Init(&argc,&argv);
  {
  pcu::PCU PCUObj;
  pumi_load_pcu(&PCUObj);
  pGeom g = pumi_geom_load(argv[1], "mesh");
  pMesh m = pumi_mesh_load(g, argv[2], 1);
  pumi_mesh_delete(m);
  pumi_geom_delete(g);
  }
  pcu::PCU_Finalize();
}

