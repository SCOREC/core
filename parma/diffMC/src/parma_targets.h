#ifndef PARMA_TARGETS_H
#define PARMA_TARGETS_H
#include <apfMesh.h>

namespace parma {
  class Sides;
  class Weights;
  class Ghosts;
  class Targets : public Associative<double> {
    public:
      Targets(Sides* s, Weights* w, double alpha);
      virtual double total()=0;
    private:
      Targets();
  };
  Targets* makeTargets(Sides* s, Weights* w, double alpha);
  Targets* makeGhostTargets(Sides* s, Weights* w, Ghost*g, double alpha);
}
#endif
