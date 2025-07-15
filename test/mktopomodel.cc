#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <gmi_null.h>
#include <gmi_mesh.h>
#include <lionPrint.h>
#include <pcu_util.h>
#include <cstdlib>

int main(int argc, char** argv)
{
  pcu::Init(&argc,&argv);
  {
  pcu::PCU pcu_obj;
  lion_set_verbosity(1);
  if ( argc != 3 ) {
    if ( !pcu_obj.Self() )
      printf("Create a topological geometric model from a mesh\n"
             "Usage: %s <mesh> <out model (.dmg)>\n", argv[0]);
    pcu::Finalize();
    exit(EXIT_FAILURE);
  }
  gmi_register_null();
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(".null", argv[1], &pcu_obj);
  apf::deriveMdsModel(m);
  gmi_model* g = m->getModel();
  gmi_write_dmg(g, argv[2]);
  m->destroyNative();
  apf::destroyMesh(m);
  }
  pcu::Finalize();
}

