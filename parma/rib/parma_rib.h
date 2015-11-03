#ifndef PARMA_RIB_H
#define PARMA_RIB_H

#include <mthVector.h>

namespace parma {

struct Body
{
  mth::Vector3<double> point;
  double mass;
};

struct Bodies
{
  Bodies(Body* arr, int n_);
  ~Bodies();
  Bodies() {}
  void destroy();
  int n;
  Body** body;
};

void bisect(Bodies* all, Bodies* left, Bodies* right);

void recursivelyBisect(Bodies* all, int depth, Bodies out[]);

}

#endif
