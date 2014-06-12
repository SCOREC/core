#include "parma_rib_math.h"
#include <apfMatrix.h>

namespace parma {

void getPrincipalEigenvector(apf::Matrix3x3 const& A, apf::Vector3& v)
{
  apf::Vector3 vs[3];
  double ls[3];
  int n = apf::eigen(A, vs, ls);
  int best = 0;
  for (int i = 1; i < n; ++i)
    if (fabs(ls[i]) > fabs(ls[best]))
      best = i;
  v = vs[best];
}

}
