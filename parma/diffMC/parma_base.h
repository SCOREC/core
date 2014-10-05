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
      Balancer(apf::Mesh* mIn, apf::MeshTag* wIn, double alphaIn,
        Sides* s, Weights* w, Targets* t, Selector* sel, 
        bool (*stop)(double imb, double maxImb)=less);
      ~Balancer();
      bool run(double maxImb, int verbosity=0);
    private:
      Balancer();
      apf::Mesh* m;
      apf::MeshTag* w;
      double alpha;
      int verbose;
      virtual double imbalance();
      Sides* sides;
      Weights* weights;
      Targets* targets;
      Selector* selects;
      bool (*stop)(double imb, double maxImb);
      static bool less(double imb, double maxImb) {
        return imb < maxImb;
      }
  };
}
#endif
