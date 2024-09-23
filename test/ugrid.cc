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

pcu::PCU* getGroupedPCU(const int partitionFactor, pcu::PCU *PCUObj)
{
  int self = PCUObj->Self();
  int groupRank = self / partitionFactor;
  int group = self % partitionFactor;
  MPI_Comm groupComm;
  MPI_Comm_split(MPI_COMM_WORLD, group, groupRank, &groupComm);
  return new pcu::PCU(groupComm);
}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
  if ( argc != 5 ) {
    if ( !PCUObj.Self() )
      printf("Usage: %s <in .[b8|lb8].ugrid> <out .dmg> <out .smb> <partition factor>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  gmi_register_null();
  const int partitionFactor = atoi(argv[4]);
  PCU_ALWAYS_ASSERT(partitionFactor <= PCUObj.Peers());
  bool isOriginal = ((PCUObj.Self() % partitionFactor) == 0);
  gmi_model* g = gmi_load(".null");
  apf::Mesh2* m = 0;
  apf::Migration* plan = 0;
  pcu::PCU *groupedPCUObj = getGroupedPCU(partitionFactor, &PCUObj);
  if (isOriginal) {
    m = apf::loadMdsFromUgrid(g, argv[1], groupedPCUObj);
    apf::deriveMdsModel(m);
    m->verify();
    plan = getPlan(m, partitionFactor);
  }
  //used switchPCU here to load the mesh on the groupedPCU, perform tasks and then call repeatMdsMesh
  //on the globalPCU
  if(m != nullptr) m->switchPCU(&PCUObj);
  delete groupedPCUObj;
  m = repeatMdsMesh(m, g, plan, partitionFactor, &PCUObj);
  Parma_PrintPtnStats(m, "");
  gmi_write_dmg(g,argv[2]);
  m->writeNative(argv[3]);
  m->destroyNative();
  apf::destroyMesh(m);
  }
  MPI_Finalize();
}

