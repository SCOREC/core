#include <gmi_mesh.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <PCU.h>
#include <parma.h>
#include <pcu_util.h>

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
    PCU_ALWAYS_ASSERT(argc==4);
    modelFile = argv[1];
    meshFile = argv[2];
  }

  apf::MeshTag* applyUnitWeight(apf::Mesh* m) {
    apf::MeshTag* wtag = m->createDoubleTag("ghostUnitWeight",1);
    apf::MeshEntity* e;
    double one = 1;
    for(int d=0; d<=m->getDimension(); d++) {
      apf::MeshIterator* itr = m->begin(d);
      while( (e = m->iterate(itr)) )
        m->setDoubleTag(e, wtag, &one);
      m->end(itr);
    }
    return wtag;
  }

  void runParma(apf::Mesh* m) {
    apf::MeshTag* weights = applyUnitWeight(m);
    double factor = 0.6;
    int verbose = 2;
    apf::Balancer* ghost = Parma_MakeGhostEdgeDiffuser(m,factor,verbose);
    ghost->balance(weights, 1.05);
    m->destroyTag(weights);
    delete ghost;
  }
}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_mesh();
  getConfig(argc,argv);
  apf::Mesh2* m = apf::loadMdsMesh(modelFile,meshFile);
  runParma(m);
  m->writeNative(argv[3]);
  freeMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}
