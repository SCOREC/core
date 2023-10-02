#ifndef PARMA_ENTWEIGHTS_H
#define PARMA_ENTWEIGHTS_H
#include "parma_weights.h"
namespace parma {
  class EntWeights : public Weights {
    public:
      EntWeights(apf::Mesh* m, apf::MeshTag* w, Sides* s, int d);
      double self();
    private:
      EntWeights();
      int entDim;
      double weight;
      virtual double getEntWeight(apf::Mesh* m, apf::MeshEntity* e, apf::MeshTag* w);
      void init(apf::Mesh* m, apf::MeshTag* w, Sides* s);
  };
}
#endif
