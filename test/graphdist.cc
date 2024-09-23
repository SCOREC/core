#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <lionPrint.h>
#include <dsp.h>
#include <dspGraphDistance.h>
#include <cstdlib>

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  {
  pcu::PCU PCUObj =  pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
  if ( argc != 4 ) {
    if ( !PCUObj.Self() )
      printf("Usage: %s <model> <mesh> <out prefix>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2],&PCUObj);
  dsp::Boundary moving;
  moving.insert(m->findModelEntity(2, 57));
  moving.insert(m->findModelEntity(2, 62));
  moving.insert(m->findModelEntity(2, 66));
  dsp::closeBoundary(m, moving);
  std::vector<apf::MeshEntity*> vs;
  dsp::getGraphDistance(m, moving, vs);
  apf::writeVtkFiles(argv[3], m);
  m->destroyNative();
  apf::destroyMesh(m);
  }
  MPI_Finalize();
}
