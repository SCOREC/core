#ifndef PARMA_BASE_H
#define PARMA_BASE_H
#include <apfMesh.h>
#include "parma_associative.h"

namespace parma {
  class Sides;
  class Weights;
  class Targets;
  class Selector;
  class Balancer {
    public:
      Balancer(apf::Mesh* m, apf::MeshTag* w, double alpha);
      void setSides(Sides* s) {sides = s;}
      void setWeights(Weights* w) {weights = w;}
      void setTargets(Targets* t) {targets = t;}
      void setSelector(Selector* s) {selects = s;}
      ~Balancer();
      bool run(double maxImb, int verbosity=0);
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
