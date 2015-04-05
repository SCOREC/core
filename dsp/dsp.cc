#include "dsp.h"
#include <apf.h>

namespace dsp {

bool tryToDisplace(apf::Mesh2* m, apf::Field* df)
{
  double f = 1;
  apf::axpy(f, df, m->getCoordinateField());
  if (0 == verifyVolumes(m, false))
    return true;
  do {
    f /= 2;
    apf::axpy(-f, df, m->getCoordinateField());
  } while (0 != verifyVolumes(m, false));
  apf::axpy(-f, df, df);
  return false;
}

void displace(apf::Mesh2* m, apf::Field* df,
    Smoother* smoother, Adapter* adapter,
    Boundary& fixed, Boundary& moving)
{
  smoother->smooth(df, fixed, moving);
  while (!tryToDisplace(m, df))
    adapter->adapt(m);
}

}
