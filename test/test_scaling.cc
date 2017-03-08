#include <pcu_util.h>
#include <iostream>

#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apf.h>
#include <PCU.h>
#include <samSz.h>
#include <samElementCount.h>

int main(int argc, char** argv)
{
  PCU_ALWAYS_ASSERT(argc==3);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  apf::Field* identity_size = samSz::isoSize(m);
  double scaling_factor = sam::getIsoLengthScalar(identity_size,
      m->count(m->getDimension()));
  std::cout << "scaling factor " << scaling_factor << '\n';
  PCU_ALWAYS_ASSERT(scaling_factor < 2.0);
  PCU_ALWAYS_ASSERT(0.5 < scaling_factor);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}
