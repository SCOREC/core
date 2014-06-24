#include <gmi_mesh.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <PCU.h>
#include <parma.h>

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
    assert(argc==3);
    modelFile = argv[1];
    meshFile = argv[2];
  }

  apf::MeshTag* applyUnitVtxWeight(apf::Mesh* m) {
    apf::MeshTag* wtag = m->createDoubleTag("hpsUnitWeight",1);
    apf::MeshEntity* e;
    apf::MeshIterator* itr = m->begin(m->getDimension());
    double w = 1;
    while( (e = m->iterate(itr)) )
      m->setDoubleTag(e, wtag, &w);
    m->end(itr);
    return wtag;
  }

  void runParma(apf::Mesh* m) {
    apf::MeshTag* weights = applyUnitVtxWeight(m);
    apf::Balancer* hps = Parma_MakeHpsBalancer(m);
    double ignored = 3.14;
    hps->balance(weights, ignored);
    delete hps;
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
  PCU_Protect();
  getConfig(argc,argv);
  apf::Mesh2* m = apf::loadMdsMesh(modelFile,meshFile);
  apf::writeVtkFiles("before", m);
  runParma(m);
  apf::writeVtkFiles("after", m);
  freeMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}


