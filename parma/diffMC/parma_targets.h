#ifndef PARMA_TARGETS_H
#define PARMA_TARGETS_H
#include <apfMesh.h>
#include "parma_associative.h"

namespace parma {
  class Sides;
  class SurfToVol;
  class Weights;
  class Ghosts;
  class Targets : public Associative<double> {
    public:
      virtual ~Targets() {}
      virtual double total()=0;
  };
  Targets* makeTargets(Sides* s, Weights* w, double alpha);
  Targets* makeVtxElmTargets(Sides* s, Weights* w[2], int aspectTol,
      double vtxTol, double alpha);
  Targets* makeWeightSideTargets(Sides* s, Weights* w, int sideTol,
      double alpha);
  Targets* makeVtxEdgeTargets(Sides* s, Weights* w[2], int sideTol,
      double vtxTol, double alpha);
  Targets* makeVtxEdgeElmTargets(Sides* s, Weights* w[3], double alpha);
  Targets* makeShapeTargets(apf::Mesh* m, Sides* s, Weights* w, double alpha,
                            double avgSideMult, double avgSide, 
                            double minSideMult, bool isInMIS);
  Targets* makeGhostTargets(Sides* s, Weights* w, Ghosts* g, double alpha);
}
#endif
