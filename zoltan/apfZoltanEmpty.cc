#include "apfZoltan.h"
#include <apf.h>

namespace apf {

Splitter* makeZoltanSplitter(Mesh* mesh, int method, int approach, bool sync)
{
  fail("apf_zoltan compiled empty !");
  return 0;
}

Balancer* makeZoltanBalancer(Mesh* mesh, int method, int approach)
{
  fail("apf_zoltan compiled empty !");
  return 0;
}

}
