#include <PCU.h>
#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <gmi_mesh.h>
#include <ma.h>
#include <maShape.h>
#include <cassert>
#include <cstdlib>

namespace {

const char* modelFile = 0;
const char* meshFile = 0;
double qualityTol = 0.0;
double minQuality = 1.0;
double avgQuality = 0.0;
double avgQualityBelowTol = 0.0;
long numElemsBelowTol = 0;

static void fail(char** argv)
{
  if (!PCU_Comm_Self())
    printf("Usage: %s <model> <mesh> <quality tolerance>\n", argv[0]);
  MPI_Finalize();
  exit(EXIT_FAILURE);
}

static void getConfig(int argc, char** argv)
{
  if (argc != 4)
    fail(argv);
  modelFile = argv[1];
  meshFile = argv[2];
  qualityTol = atof(argv[3]);
}

static void parallelReduce(apf::Mesh2* m)
{
  long numElems = m->count(m->getDimension());
  numElems = PCU_Add_Long(numElems);
  avgQuality = PCU_Add_Double(avgQuality);
  avgQuality /= numElems;
  numElemsBelowTol = PCU_Add_Long(numElemsBelowTol);
  avgQualityBelowTol = PCU_Add_Double(avgQualityBelowTol);
  if (numElemsBelowTol > 0)
    avgQualityBelowTol /= numElemsBelowTol;
  minQuality = PCU_Min_Double(minQuality);
}

static void processMesh(apf::Mesh2* m)
{
  int d = m->getDimension();
  apf::MeshEntity* elem;
  apf::MeshIterator* elems = m->begin(d);
  ma::IdentitySizeField I(m);
  while ((elem = m->iterate(elems)))
  {
    double q = ma::measureElementQuality(m, &I, elem);
    q = pow(q, (1.0/d));
    avgQuality += q;
    minQuality = fmin(minQuality, q);
    if (q < qualityTol)
    {
      numElemsBelowTol++;
      avgQualityBelowTol += q;
    }
  }
  parallelReduce(m);
  m->end(elems);
}

void printDiagnostics()
{
  if (!PCU_Comm_Self())
  {
    printf("average element quality: %f\n", avgQuality);
    printf("quality tolerance: %f\n", qualityTol);
    printf("number of elements below tolerance: %ld\n", numElemsBelowTol);
    if (numElemsBelowTol > 0)
      printf("avgerage quality of elements below tolerance: %f\n",
          avgQualityBelowTol);
    printf("minumum quality in mesh: %f\n", minQuality);
  }
}

}

int main(int argc, char** argv)
{
  assert(argc==4);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_mesh();
  getConfig(argc,argv);
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  processMesh(m);
  printDiagnostics();
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}
