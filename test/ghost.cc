#include <gmi_mesh.h>
#include <ma.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <PCU.h>
#include <parma.h>
#include <cassert>

namespace {
  const char* modelFile = 0;
  const char* meshFile = 0;

  void freeMesh(apf::Mesh* m)
  {
    m->destroyNative();
    apf::destroyMesh(m);
  }

  void getConfig(int argc, char** argv)
  {
    assert(argc==4);
    modelFile = argv[1];
    meshFile = argv[2];
  }

  apf::MeshTag* applyUnitWeight(apf::Mesh* m) {
    apf::MeshTag* wtag = m->createDoubleTag("ghostUnitWeight",1);
    apf::MeshEntity* e;
    for(int d=0; d <= m->getDimension(); d++) {
      apf::MeshIterator* itr = m->begin(d);
      double w = 1;
      while( (e = m->iterate(itr)) )
        m->setDoubleTag(e, wtag, &w);
      m->end(itr);
    }
    return wtag;
  }

  void runParma(apf::Mesh* m) {
    apf::MeshTag* weights = applyUnitWeight(m);
    const int layers = 1;
    const int bridgeDim = 0;
    const double stepFactor = 0.5;
    const int verbosity = 2;
    apf::Balancer* ghost =
      Parma_MakeGhostDiffuser(m, layers, bridgeDim, stepFactor, verbosity);
    ghost->balance(weights, 1.05);
    m->destroyTag(weights);
    delete ghost;
  }
}

int main(int argc, char** argv)
{
  int provided;
  MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided);
  assert(provided==MPI_THREAD_MULTIPLE);
  PCU_Comm_Init();
  PCU_Debug_Open();
  gmi_register_mesh();
  getConfig(argc,argv);
  apf::Mesh2* m = apf::loadMdsMesh(modelFile,meshFile);
  runParma(m);
  apf::writeVtkFiles(argv[3],m);
  freeMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}
