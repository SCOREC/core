#ifndef PARMA_WEIGHTS_H
#define PARMA_WEIGHTS_H
#include <apfMesh.h>
#include "parma_associative.h"

namespace parma {
  class Sides;
  class Weights : public Associative<double> {
    public:
      Weights(apf::Mesh*, apf::MeshTag*, Sides*) {}
      virtual ~Weights() {}
      virtual double self()=0;
    private:
      Weights();
  };
  class GhostWeights;
  Weights* makeEntWeights(apf::Mesh* m, apf::MeshTag* w, Sides* s, int dim);
  Weights* makeGhostMPASWeights(apf::Mesh* m, apf::MeshTag* w, Sides* s,
      int layers, int bridge);
  GhostWeights* makeVtxGhostWeights(apf::Mesh* m, apf::MeshTag* w, Sides* s,
      int layers);
  GhostWeights* makeElmGhostWeights(apf::Mesh* m, apf::MeshTag* w, Sides* s);
  void destroyGhostWeights(GhostWeights* gw);
  Weights* convertGhostToEntWeight(GhostWeights* gw, int dim);
  double getEntWeight(apf::Mesh* m, apf::MeshEntity* e, apf::MeshTag* w);
  double getMaxWeight(apf::Mesh* m, apf::MeshTag* w, int entDim);
  double getAvgWeight(apf::Mesh* m, apf::MeshTag* w, int entDim);
  double getWeight(apf::Mesh* m, apf::MeshTag* w, int entDim);
  double getMaxWeight(Weights* w, pcu::PCU *PCUObj);
  void getImbalance(Weights* w, double& imb, double& avg, pcu::PCU *PCUObj);
}

#endif

