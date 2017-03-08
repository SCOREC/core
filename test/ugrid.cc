#include <apf.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
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

void switchToOriginals(const int partitionFactor)
{
  int self = PCU_Comm_Self();
  int groupRank = self / partitionFactor;
  int group = self % partitionFactor;
  MPI_Comm groupComm;
  MPI_Comm_split(MPI_COMM_WORLD, group, groupRank, &groupComm);
  PCU_Switch_Comm(groupComm);
}

void switchToAll()
{
  MPI_Comm prevComm = PCU_Get_Comm();
  PCU_Switch_Comm(MPI_COMM_WORLD);
  MPI_Comm_free(&prevComm);
  PCU_Barrier();
}
int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  if ( argc != 5 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <in .[b8|lb8].ugrid> <out .dmg> <out .smb> <partition factor>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  gmi_register_null();
  const int partitionFactor = atoi(argv[4]);
  PCU_ALWAYS_ASSERT(partitionFactor <= PCU_Comm_Peers());
  bool isOriginal = ((PCU_Comm_Self() % partitionFactor) == 0);
  gmi_model* g = gmi_load(".null");
  apf::Mesh2* m = 0;
  apf::Migration* plan = 0;
  switchToOriginals(partitionFactor);
  if (isOriginal) {
    m = apf::loadMdsFromUgrid(g, argv[1]);
    apf::deriveMdsModel(m);
    m->verify();
    plan = getPlan(m, partitionFactor);
  }
  switchToAll();
  m = repeatMdsMesh(m, g, plan, partitionFactor);
  Parma_PrintPtnStats(m, "");
  gmi_write_dmg(g,argv[2]);
  m->writeNative(argv[3]);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

