#ifndef PH_ADAPT_H
#define PH_ADAPT_H

#include "phInput.h"

namespace apf {
class Mesh2;
}

namespace ph {

enum {
  PH_STRATEGY_UR = 7,
  PH_STRATEGIES = 9
};

void adapt(Input& in, apf::Mesh2* m);

void tetrahedronize(Input& in, apf::Mesh2* m);

}

#endif
