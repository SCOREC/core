#ifndef PARMA_STEP_H
#define PARMA_STEP_H
#include <apfMesh.h>
#include "parma_associative.h"
#include "parma_stop.h"

namespace parma {
  class Sides;
  class Weights;
  class Targets;
  class Selector;
  class Stepper {
    public:
      Stepper(apf::Mesh* mIn, double alphaIn,
        Sides* s, Weights* w, Targets* t, Selector* sel,
        const char* entType, Stop* stopper = new Less);
      virtual ~Stepper();
      bool step(double maxImb, int verbosity=0);
    private:
      Stepper();
      apf::Mesh* m;
      double alpha;
      int verbose;
      Sides* sides;
      Weights* weights;
      Targets* targets;
      Selector* selects;
      const char* name;
      Stop* stop;
  };
}
#endif
