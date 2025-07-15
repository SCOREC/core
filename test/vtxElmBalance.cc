#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <gmi_mesh.h>
#include <parma.h>
#include <lionPrint.h>
#include <pcu_util.h>
#include <cstdlib>

namespace {
  void setWeight(apf::Mesh* m, apf::MeshTag* tag, int dim) {
    double w = 1.0;
    apf::MeshEntity* e;
    apf::MeshIterator* it = m->begin(dim);
    while ((e = m->iterate(it)))
      m->setDoubleTag(e, tag, &w);
    m->end(it);
  }

  apf::MeshTag* setWeights(apf::Mesh* m) {
    apf::MeshTag* tag = m->createDoubleTag("parma_weight", 1);
    setWeight(m, tag, 0);
    setWeight(m, tag, m->getDimension());
    return tag;
  }

  void clearTags(apf::Mesh* m, apf::MeshTag* t) {
    apf::removeTagFromDimension(m, t, 0);
    apf::removeTagFromDimension(m, t, m->getDimension());
  }
}

int main(int argc, char** argv)
{
  PCU_ALWAYS_ASSERT(argc == 4);
  pcu::Init(&argc,&argv);
  {
  pcu::PCU PCUObj;
  lion_set_verbosity(1);
  if ( argc != 4 ) {
    if ( !PCUObj.Self() )
      printf("Usage: %s <model> <mesh> <out mesh>\n", argv[0]);
    pcu::Finalize();
    exit(EXIT_FAILURE);
  }
  gmi_register_mesh();
  //load model and mesh
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2],&PCUObj);
  Parma_PrintPtnStats(m, "initial");
  apf::MeshTag* weights = setWeights(m);
  const double step = 0.5; const int verbose = 2;
  apf::Balancer* balancer = Parma_MakeVtxElmBalancer(m, step, verbose);
  balancer->balance(weights, 1.05);
  delete balancer;
  clearTags(m, weights);
  m->destroyTag(weights);
  Parma_PrintPtnStats(m, "final");
  m->writeNative(argv[3]);
  // destroy mds
  m->destroyNative();
  apf::destroyMesh(m);
  }
  pcu::Finalize();
}
