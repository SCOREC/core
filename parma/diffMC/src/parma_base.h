#ifndef PARMA_BASE_H
#define PARMA_BASE_H
#include <apfMesh.h>
#include "parma_associative.h"
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
#include "parma_selector.h"

namespace parma {
  class Balancer {
    public:
      Balancer(apf::Mesh* mIn, apf::MeshTag* wIn, double alphaIn);
      void setSelector(Selector* s);
      void setWeights(Weights* w);
      void setTargets(Targets* t);
      ~Balancer();
      bool run(double maxImb, int verbosity);
    private:
      Balancer();
      apf::Mesh* m;
      double maxImb;
      apf::MeshTag* w;
      double alpha;
      int verbose;
      virtual double imbalance();
      Sides* sides;
      Weights* weights;
      Targets* targets;
      Selector* selects;
  };
};
#endif
