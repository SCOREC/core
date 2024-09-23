#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfBox.h>
#include <apfMesh2.h>
#include <apf.h>
#include <lionPrint.h>
#include <pcu_util.h>
#include <cstdlib>

namespace {

int nx = 0;
int ny = 0;
int nz = 0;
double wx = 0.0;
double wy = 0.0;
double wz = 0.0;
bool is = true;
const char* modelFile = 0;
const char* meshFile = 0;

void verifyArgs(int argc, char** argv, pcu::PCU *PCUObj)
{
  if (argc != 10)
  {
    if (!PCUObj->Self())
    {
      printf("Usage: %s <nx> <ny> <nz> <wx> "
             "<wy> <wz> <is> <model> <mesh>\n", argv[0]);
      printf(" <nx> num x elements\n");
      printf(" <ny> num y elements\n");
      printf(" <nz> num z elements\n");
      printf(" <wz> x length scale\n");
      printf(" <wy> y length scale\n");
      printf(" <wz> z length scale\n");
      printf(" <is> is simplical mesh\n");
      printf(" <model> .dmg model file name\n");
      printf(" <mesh> .smb mesh file name\n");
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
  }
}

void getArgs(char** argv)
{
  nx = atoi(argv[1]);
  ny = atoi(argv[2]);
  nz = atoi(argv[3]);
  wx = atof(argv[4]);
  wy = atof(argv[5]);
  wz = atof(argv[6]);
  is = atoi(argv[7]);
  modelFile = argv[8];
  meshFile = argv[9];
}

}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  {
  pcu::PCU pcu_obj = pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
  verifyArgs(argc, argv, &pcu_obj);
  if(pcu_obj.Peers()>1) {
    if(pcu_obj.Self())
      fprintf(stderr, "%s must be run on a single process!\n", argv[0]);
    exit(1);
  }
  getArgs(argv);
  gmi_register_mesh();
  apf::Mesh2* m = apf::makeMdsBox(nx,ny,nz,wx,wy,wz,is,&pcu_obj);
  gmi_model* g = m->getModel();
  m->verify();
  m->writeNative(meshFile);
  gmi_write_dmg(g, modelFile);
  m->destroyNative();
  apf::destroyMesh(m);
  }
  MPI_Finalize();
}
