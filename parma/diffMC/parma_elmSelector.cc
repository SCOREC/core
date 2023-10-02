#include "parma_vtxSelector.h"

namespace {
  class ElmSelector : public parma::VtxSelector {
    public:
      ElmSelector(apf::Mesh* m, apf::MeshTag* w)
        : VtxSelector(m, w) {}
    protected:
      virtual double add(apf::MeshEntity*, apf::Up& cavity, const int destPid,
          apf::Migration* plan) {
        double w = 0;
        for(int i=0; i < cavity.n; i++) {
          plan->send(cavity.e[i], destPid);
          w += getWeight(cavity.e[i]);
        }
        return w;
      }
  };
}

namespace parma {
  Selector* makeElmSelector(apf::Mesh* m, apf::MeshTag* w) {
    return new ElmSelector(m, w);
  }
}
