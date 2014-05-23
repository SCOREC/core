#include <gmi_mesh.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <PCU.h>
#include <parma.h>

namespace {

const char* modelFile = 0;
const char* meshFile = 0;
const char* outFile = 0;
int partitionFactor = 1;

apf::Mesh2* readMesh()
{
  double t0 = MPI_Wtime();
  apf::Mesh2* m = apf::loadMdsMesh(modelFile,meshFile);
  double t1 = MPI_Wtime();
  if ( ! PCU_Comm_Self())
    printf("time to load %s and %s: %f seconds\n",modelFile,meshFile,t1-t0);
  return m;
}

void freeMesh(apf::Mesh* m)
{
  m->destroyNative();
  apf::destroyMesh(m);
}

apf::Migration* globalMigr;

void preparePartition(apf::Mesh* m)
{
  apf::Splitter* splitter = Parma_MakeRibSplitter(m);
  double t0 = MPI_Wtime();
  globalMigr = splitter->split(0, 0, partitionFactor);
  double t1 = MPI_Wtime();
  if ( ! PCU_Comm_Self())
    printf("time to run RIB: %f seconds\n",t1-t0);
  delete splitter;
}

apf::Mesh2* globalMesh;

void* thread_main(void*)
{
  int threadSelf = PCU_Thrd_Self();
  apf::Mesh* m;
  apf::Migration* plan;
  if (!threadSelf) {
    m = globalMesh;
    plan = globalMigr;
  } else {
    m = apf::makeEmptyMdsMesh(getMdsModel(globalMesh),
        globalMesh->getDimension(), globalMesh->hasMatching());
    plan = new apf::Migration(m);
  }
  double t0 = MPI_Wtime();
  m->migrate(plan);
  double t1 = MPI_Wtime();
  if ( ! PCU_Comm_Self())
    printf("time to migrate: %f seconds\n",t1-t0);
  m->verify();
  double t2 = MPI_Wtime();
  m->writeNative(outFile);
  double t3 = MPI_Wtime();
  if ( ! PCU_Comm_Self())
    printf("time to write %s: %f seconds\n",outFile,t3-t2);
  freeMesh(m);
  return NULL;
}

void getConfig(int argc, char** argv)
{
  assert(argc==5);
  modelFile = argv[1];
  meshFile = argv[2];
  outFile = argv[3];
  partitionFactor = atoi(argv[4]);
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
  apf::Mesh2* m = readMesh();
  preparePartition(m);
  globalMesh = m;
  PCU_Thrd_Run(partitionFactor, thread_main, NULL);
  PCU_Comm_Free();
  MPI_Finalize();
}


