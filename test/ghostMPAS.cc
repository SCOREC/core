#include <gmi_mesh.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <lionPrint.h>
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
    double w = 1;
    int dims[2] = {0,m->getDimension()};
    const size_t len = sizeof(dims)/sizeof(int);
    for(size_t i=0; i<len; i++) {
      apf::MeshIterator* itr = m->begin(dims[i]);
      while( (e = m->iterate(itr)) )
        m->setDoubleTag(e, wtag, &w);
      m->end(itr);
    }
    return wtag;
  }

  void runParma(apf::Mesh* m) {
    apf::MeshTag* weights = applyUnitWeight(m);
    const int layers = 3;
    const int bridgeDim = 1;
    apf::Balancer* ghost = Parma_MakeMPASDiffuser(m, layers, bridgeDim);
    ghost->balance(weights, 1.01);
    m->destroyTag(weights);
    delete ghost;
  }
}

int main(int argc, char** argv)
{
  pcu::Init(&argc,&argv);
  {
  pcu::PCU PCUObj;
  lion_set_verbosity(1);
  gmi_register_mesh();
  getConfig(argc,argv);
  apf::Mesh2* m = apf::loadMdsMesh(modelFile,meshFile,&PCUObj);
  runParma(m);
  m->writeNative(argv[3]);
  freeMesh(m);
  }
  pcu::Finalize();
}
