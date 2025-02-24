#include <apf.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <lionPrint.h>
#include <parma.h>
#include <cstdlib>
#include <pcu_util.h>

apf::Migration* getPlan(apf::Mesh* m, const int partitionFactor)
{
  apf::Splitter* splitter = Parma_MakeRibSplitter(m);
  apf::MeshTag* weights = Parma_WeighByMemory(m);
  apf::Migration* plan = splitter->split(weights, 1.10, partitionFactor);
  apf::removeTagFromDimension(m, weights, m->getDimension());
  m->destroyTag(weights);
  delete splitter;
  return plan;
}

int main(int argc, char** argv)
{
  pcu::PCU_Init(&argc,&argv);
  {
  pcu::PCU PCUObj;
  lion_set_verbosity(1);
  if ( argc != 5 ) {
    if ( !PCUObj.Self() )
      printf("Usage: %s <in .[b8|lb8].ugrid> <out .dmg> <out .smb> <partition factor>\n", argv[0]);
    pcu::PCU_Finalize();
    exit(EXIT_FAILURE);
  }
  gmi_register_null();
  const int partitionFactor = atoi(argv[4]);
  PCU_ALWAYS_ASSERT(partitionFactor <= PCUObj.Peers());
  bool isOriginal = ((PCUObj.Self() % partitionFactor) == 0);
  gmi_model* g = gmi_load(".null");
  apf::Mesh2* m = 0;
  apf::Migration* plan = 0;
  auto groupedPCUObj = PCUObj.Split(
    PCUObj.Self() % partitionFactor, PCUObj.Self() / partitionFactor
  );
  if (isOriginal) {
    m = apf::loadMdsFromUgrid(g, argv[1], groupedPCUObj.get());
    apf::deriveMdsModel(m);
    m->verify();
    plan = getPlan(m, partitionFactor);
  }
  //used switchPCU here to load the mesh on the groupedPCU, perform tasks and then call repeatMdsMesh
  //on the globalPCU
  if(m != nullptr) m->switchPCU(&PCUObj);
  m = repeatMdsMesh(m, g, plan, partitionFactor, &PCUObj);
  Parma_PrintPtnStats(m, "");
  gmi_write_dmg(g,argv[2]);
  m->writeNative(argv[3]);
  m->destroyNative();
  apf::destroyMesh(m);
  }
  pcu::PCU_Finalize();
}

