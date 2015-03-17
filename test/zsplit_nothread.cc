#include <gmi_mesh.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <PCU.h>
#include <parma.h>
#include <apfZoltan.h>

namespace {

const char* modelFile = 0;
const char* meshFile = 0;
const char* outFile = 0;
int partitionFactor = 1;

void freeMesh(apf::Mesh* m)
{
  m->destroyNative();
  apf::destroyMesh(m);
}

apf::Migration* getPlan(apf::Mesh* m)
{
  apf::Splitter* splitter = apf::makeZoltanSplitter(
      m, apf::GRAPH, apf::PARTITION, false);
  apf::MeshTag* weights = Parma_WeighByMemory(m);
  apf::Migration* plan = splitter->split(weights, 1.05, partitionFactor);
  apf::removeTagFromDimension(m, weights, m->getDimension());
  m->destroyTag(weights);
  delete splitter;
  return plan;
}

void switchToMasters()
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

void remapMesh(apf::Mesh2* m)
{
  apf::Multiply remap(partitionFactor);
  apf::remapPartition(m, remap);
}

void mymain(bool ismaster)
{
  gmi_model* g = 0;
  apf::Mesh2* m = 0;
  apf::Migration* plan = 0;
  int dim, matched;
  switchToMasters();
  if (ismaster) {
    m = apf::loadMdsMesh(modelFile,meshFile);
    dim = m->getDimension();
    matched = m->hasMatching();
    plan = getPlan(m);
  } else {
    g = gmi_load(modelFile);
  }
  switchToAll();
  PCU_Comm_Begin();
  if (ismaster)
    for (int i = 1; i < partitionFactor; ++i) {
      PCU_COMM_PACK(PCU_Comm_Self() + i, dim);
      PCU_COMM_PACK(PCU_Comm_Self() + i, matched);
    }
  PCU_Comm_Send();
  while (PCU_Comm_Receive()) {
    PCU_COMM_UNPACK(dim);
    PCU_COMM_UNPACK(matched);
  }
  if (!ismaster) {
    m = apf::makeEmptyMdsMesh(g, dim, matched);
    plan = new apf::Migration(m);
  }
  remapMesh(m);
  m->migrate(plan);
  m->writeNative(outFile);
  freeMesh(m);
}

void getConfig(int argc, char** argv)
{
  if ( argc != 5 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh> <outMesh> <factor>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  modelFile = argv[1];
  meshFile = argv[2];
  outFile = argv[3];
  partitionFactor = atoi(argv[4]);
  assert(partitionFactor <= PCU_Comm_Peers());
}

}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_mesh();
  getConfig(argc,argv);
  if (PCU_Comm_Self() % partitionFactor)
    mymain(false);
  else
    mymain(true);
  PCU_Comm_Free();
  MPI_Finalize();
}

