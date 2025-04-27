#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <PCU.h>
#include <pcu_util.h>
#include <stdlib.h>

int main(int ac, char * av[])
{

  MPI_Init(&ac,&av);
  PCU_Comm_Init();
  if (ac < 2)
  {
    if (PCU_Comm_Self() == 0)
      printf("USAGE1: %s <mesh.smb>\n", av[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  gmi_register_mesh();
  gmi_register_null();

  apf::Mesh * msh = apf::loadMdsMesh(".null",av[1]);

  int lcl_vrt = apf::countEntitiesOfType(msh,apf::Mesh::VERTEX);
  int own_vrt = apf::countEntitiesOfType(msh,apf::Mesh::VERTEX,apf::getSharing(msh));

  int lcl_vrt2 = apf::countOwned(msh,0,apf::getNoSharing());
  int own_vrt2 = apf::countOwned(msh,0);

  int lcl_edg = apf::countEntitiesOfType(msh,apf::Mesh::EDGE);
  int own_edg = apf::countEntitiesOfType(msh,apf::Mesh::EDGE,apf::getSharing(msh));

  int lcl_edg2 = apf::countOwned(msh,1,apf::getNoSharing());
  int own_edg2 = apf::countOwned(msh,1);

  // higher dimension mesh entities is where counting specific types gets tricky
  // using the default interfaces,
  // which is the entire reason for this function in the first place.

  PCU_ALWAYS_ASSERT((lcl_vrt == lcl_vrt2) && (own_vrt == own_vrt2));
  PCU_ALWAYS_ASSERT((lcl_edg == lcl_edg2) && (own_edg == own_edg2));

  apf::Field * fld = apf::createPackedField(msh,"scalar",1);
  apf::Sharing * shr = apf::getSharing(msh);

  int lcl_nds = apf::countLocalNodes(fld);
  int own_nds = apf::countLocalNodes(fld,shr);

  PCU_ALWAYS_ASSERT((lcl_vrt == lcl_nds) && (own_vrt == own_nds));

  int gbl_nds = apf::countGlobalNodes(fld);
  int all_nds = apf::countGlobalNodes(fld,apf::getNoSharing());

  int gbl_nds2 = PCU_Add_Int(own_nds);
  int all_nds2 = PCU_Add_Int(lcl_nds);

  PCU_ALWAYS_ASSERT((gbl_nds == gbl_nds2) && (all_nds == all_nds2));

  delete shr;

  msh->destroyNative();
  apf::destroyMesh(msh);

  PCU_Comm_Free();
  MPI_Finalize();
}

