#ifndef PARMA_TARGETS_H
#define PARMA_TARGETS_H
#include <apfMesh.h>
#include "parma_associative.h"

namespace parma {
  class Sides;
  class Weights;
  class Ghosts;
  class Targets : public Associative<double> {
    public:
      virtual ~Targets() {}
      virtual double total()=0;
  };
  Targets* makeTargets(Sides* s, Weights* w, double alpha);
  Targets* makeVtxEdgeTargets(Sides* s, Weights* w[2], double alpha);
  Targets* makeVtxEdgeElmTargets(Sides* s, Weights* w[3], double alpha);
  Targets* makeShapeTargets(apf::Mesh* m, Sides* s, Weights* w, double alpha);
  Targets* makeGhostTargets(Sides* s, Weights* w, Ghosts* g, double alpha);
}
#endif
