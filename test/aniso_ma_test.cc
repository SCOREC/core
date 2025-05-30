#include "ma.h"
#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <lionPrint.h>
#include <pcu_util.h>

#include <stdlib.h>

int main(int argc, char** argv)
{
  PCU_ALWAYS_ASSERT(argc==4);
  const char* modelFile = argv[1];
  const char* meshFile = argv[2];
  pcu::Init(&argc,&argv);
  {
  pcu::PCU PCUObj;
  lion_set_verbosity(1);
  gmi_register_mesh();
  ma::Mesh* m = apf::loadMdsMesh(modelFile,meshFile,&PCUObj);
  m->verify();
  apf::writeVtkFiles("aniso_before",m);
  AnIso sf(m);
  ma::Input* in = ma::makeAdvanced(ma::configure(m, &sf, 0, logInterpolation));
  in->goodQuality = 0.2;
  ma::adapt(in);
  m->verify();
  if (logInterpolation)
    apf::writeVtkFiles("aniso_log_interpolation_after",m);
  else
    apf::writeVtkFiles("aniso_after",m);
  m->destroyNative();
  apf::destroyMesh(m);
  }
  pcu::Finalize();
}

