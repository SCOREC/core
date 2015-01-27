#ifndef PARMA_VTXSELECTOR_H
#define PARMA_VTXSELECTOR_H

#include <map>
#include "parma_selector.h"

namespace parma {
  typedef std::map<int,double> Mid;

  class VtxSelector : public Selector {
    public:
      VtxSelector(apf::Mesh* m, apf::MeshTag* w);
      ~VtxSelector();
      apf::Migration* run(Targets* tgts);
    protected:
      virtual double getWeight(apf::MeshEntity* e);
      virtual double add(apf::MeshEntity* v, apf::Up& cavity, const int destPid,
          apf::Migration* plan);
      double select(Targets* tgts, apf::Migration* plan, double planW,
          int maxSize);
    private:
      apf::MeshTag* dist;
      VtxSelector();
      Mid sending;
  };
}

#endif
