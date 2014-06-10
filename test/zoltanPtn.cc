#include <gmi_mesh.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <parma.h>
#include <PCU.h>
#include <apfPartition.h>
#include <apfZoltan.h>
#include <stdlib.h>
#include <stdio.h>

void applyElmWeights(apf::Mesh* m, apf::MeshTag* wtag) { 
  assert(wtag); 
  apf::MeshEntity* e; 
  apf::MeshIterator* itr = m->begin(m->getDimension()); 
  double w = 1; 
  while( (e = m->iterate(itr)) )  
    m->setDoubleTag(e, wtag, &w); 
  m->end(itr); 
}

static std::string outMeshFN;

void getStats(apf::Mesh2* m) {
  m->verify();
  Parma_PrintPtnStats(m, "split");
  apf::writeVtkFiles("split",m);
  if(outMeshFN.size())
    m->writeNative(outMeshFN.c_str());
  m->destroyNative();
  apf::destroyMesh(m);
}

int main(int argc, char** argv) {

  int provided;
  MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided);
  assert(provided == MPI_THREAD_MULTIPLE);
  PCU_Comm_Init();

  if(argc < 3 || argc > 4) {
    fprintf(stderr, "Usage: %s <model> <mesh> [<outmesh>]\n", argv[0]);
    return 1;
  }

  if(argc==3) outMeshFN = "";
  if(argc==4) outMeshFN = argv[3];

  gmi_register_mesh();
  PCU_Protect();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1], argv[2]);
  m->verify();

  apf::MeshTag* w = m->createDoubleTag("testZtnWeight",1);
  applyElmWeights(m,w);

  Parma_PrintPtnStats(m, "before");
  apf::writeVtkFiles("before",m);
  const double tol = 1.05;
  apf::Balancer* b = apf::makeZoltanBalancer(m, apf::GRAPH, apf::REPARTITION);
  b->balance(w, tol);
  delete b;
  m->verify();
  Parma_PrintPtnStats(m, "after");
  apf::writeVtkFiles("balance",m);

  const int mult = 2;
  apf::Splitter* s = apf::makeZoltanSplitter(m, apf::RIB, apf::REPARTITION);
  double startTime = MPI_Wtime();
  apf::Migration* plan = s->split(w, tol, mult);
  double elapsedTime = MPI_Wtime() - startTime;
  PCU_Max_Doubles(&elapsedTime,1);
  if (PCU_Comm_Self()==0)
    printf("Partition time - %f seconds\n",elapsedTime);
  delete s;
  apf::removeTagFromDimension(m, w, m->getDimension());
  m->destroyTag(w);
  splitMdsMesh(m, plan, mult, getStats);

  PCU_Comm_Free();
  return 0;
}
