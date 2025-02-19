#include <apf.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <lionPrint.h>
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
#ifndef SCOREC_NO_MPI
  MPI_Init(&argc,&argv);
#else
  (void) argc, (void) argv;
#endif
  {
  pcu::PCU PCUObj;
  lion_set_verbosity(1);
  gmi_register_null();
  PCU_ALWAYS_ASSERT( 3 == argc );
  const char* ugridfile = argv[1];
  const char* ptnfile = argv[2];
  gmi_register_null();
  gmi_model* g = gmi_load(".null");
  apf::printUgridPtnStats(g,ugridfile,ptnfile,weights,&PCUObj);
  }
#ifndef SCOREC_NO_MPI
  MPI_Finalize();
#endif
}
