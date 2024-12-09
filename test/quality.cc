#include <lionPrint.h>
#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <gmi_mesh.h>
#include <ma.h>
#include <maShape.h>
#include <pcu_util.h>
#include <cstdlib>

namespace {

const char* modelFile = 0;
const char* meshFile = 0;
double qualityTol = 0.0;
double minQuality = 1.0;
double avgQuality = 0.0;
double avgQualityBelowTol = 0.0;
long numElemsBelowTol = 0;

static void fail(char** argv, pcu::PCU *PCUObj)
{
  if (!PCUObj->Self())
    printf("Usage: %s <model> <mesh> <quality tolerance>\n", argv[0]);
  MPI_Finalize();
  exit(EXIT_FAILURE);
}

static void getConfig(int argc, char** argv, pcu::PCU *PCUObj)
{
  if (argc != 4)
    fail(argv, PCUObj);
  modelFile = argv[1];
  meshFile = argv[2];
  qualityTol = atof(argv[3]);
}

static void parallelReduce(apf::Mesh2* m)
{
  long numElems = m->count(m->getDimension());
  numElems = m->getPCU()->Add<long>(numElems);
  avgQuality = m->getPCU()->Add<double>(avgQuality);
  avgQuality /= numElems;
  numElemsBelowTol = m->getPCU()->Add<long>(numElemsBelowTol);
  avgQualityBelowTol = m->getPCU()->Add<double>(avgQualityBelowTol);
  if (numElemsBelowTol > 0)
    avgQualityBelowTol /= numElemsBelowTol;
  minQuality = m->getPCU()->Min<double>(minQuality);
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

void printDiagnostics(pcu::PCU *PCUObj)
{
  if (!PCUObj->Self())
  {
    printf("average element quality: %f\n", avgQuality);
    printf("quality tolerance: %f\n", qualityTol);
    printf("number of elements below tolerance: %ld\n", numElemsBelowTol);
    if (numElemsBelowTol > 0)
      printf("average quality of elements below tolerance: %f\n",
          avgQualityBelowTol);
    printf("minimum quality in mesh: %f\n", minQuality);
  }
}

}

int main(int argc, char** argv)
{
  PCU_ALWAYS_ASSERT(argc==4);
  MPI_Init(&argc,&argv);
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
  gmi_register_mesh();
  getConfig(argc,argv,&PCUObj);
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2],&PCUObj);
  processMesh(m);
  printDiagnostics(&PCUObj);
  m->destroyNative();
  apf::destroyMesh(m);
  }
  MPI_Finalize();
}
