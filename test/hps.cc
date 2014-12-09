#include <gmi_mesh.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <PCU.h>
#include <parma.h>

namespace {
  apf::MeshTag* applyUnitWeight(apf::Mesh* m) {
    apf::MeshTag* wtag = m->createDoubleTag("hpsUnitWeight",1);
    apf::MeshEntity* e;
    apf::MeshIterator* itr = m->begin(m->getDimension());
    double w = 1.0;
    //TODO Remove if not testing with torus/4imb
    if(PCU_Comm_Self() == 0) w = 1.8193;
    else if (PCU_Comm_Self() == 3) w = .804069;
    //END
    while( (e = m->iterate(itr)) )
      m->setDoubleTag(e, wtag, &w);
    m->end(itr);
    return wtag;
  }
}

int main(int argc, char** argv) {
  assert(argc == 4);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  PCU_Protect();
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  apf::MeshTag* weights = applyUnitWeight(m);
  apf::Balancer* hps = Parma_MakeHpsBalancer(m);
  double ignored = 3.14;
  hps->balance(weights, ignored);
  delete hps;
  removeTagFromDimension(m, weights, m->getDimension());
  m->destroyTag(weights);
  Parma_PrintPtnStats(m, "final");
  m->writeNative(argv[3]);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}
