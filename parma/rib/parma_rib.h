#ifndef PARMA_RIB_H
#define PARMA_RIB_H

#include <apfVector.h>

namespace parma {

struct Body
{
  apf::Vector3 point;
  double mass;
};

struct Bodies
{
  Bodies(Body* arr, int n_);
  Bodies() {}
  void destroy();
  int n;
  Body** body;
};

void testBisectionPlane(Body* b, int n);

void bisect(Bodies* all, Bodies* left, Bodies* right);

void recursivelyBisect(Bodies* all, int depth, Bodies out[]);

}

#endif
