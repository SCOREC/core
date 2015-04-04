#include "dsp.h"

namespace dsp {

bool tryToDisplace(apf::Mesh2* m, apf::Field* df)
{
  double f = 1;
  displaceMesh(m, df, f);
  if (0 == verifyVolumes(m, false))
    return true;
  do {
    f /= 2;
    displaceMesh(m, df, -f);
  } while (0 != verifyVolumes(m, false));
  return false;
}

}
