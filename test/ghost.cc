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

apf::MeshTag* applyUnitElmWeight(apf::Mesh* m) {
  apf::MeshTag* wtag = m->createDoubleTag("ghostUnitWeight",1);
  apf::MeshEntity* e;
  apf::MeshIterator* itr = m->begin(m->getDimension());
  double w = 1;
  while( (e = m->iterate(itr)) )
    m->setDoubleTag(e, wtag, &w);
  m->end(itr);
  return wtag;
}

void runParma(apf::Mesh* m) {
  apf::MeshTag* weights = applyUnitElmWeight(m);
  const int layers = 3;
  const int bridgeDim = 1;
  Parma_RunGhostPtnImprovement(m, weights, 1.05, layers, bridgeDim, 1);
}

}

int main(int argc, char** argv)
{
  int provided;
  MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided);
  assert(provided==MPI_THREAD_MULTIPLE);
  PCU_Comm_Init();
  gmi_register_mesh();
  PCU_Protect();
  getConfig(argc,argv);
  apf::Mesh2* m = apf::loadMdsMesh(modelFile,meshFile);
  runParma(m);
  freeMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}


