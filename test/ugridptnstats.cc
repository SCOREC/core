#include <apf.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <cstdlib>
#include <pcu_util.h>

const double vtxw = 1.0;
const double edgew = 1.0;
const double triw = 1.0;
const double quadw = 1.0;
const double tetw = 1.0;
const double pyrw = 6.8;
const double przw = 7.5;
const double hexw = 13.8;
const double weights[8] = {vtxw, edgew, triw, quadw, tetw, hexw, przw, pyrw};

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_null();
  PCU_ALWAYS_ASSERT( 3 == argc );
  const char* ugridfile = argv[1];
  const char* ptnfile = argv[2];
  gmi_register_null();
  gmi_model* g = gmi_load(".null");
  apf::printUgridPtnStats(g,ugridfile,ptnfile,weights);
  PCU_Comm_Free();
  MPI_Finalize();
}
