#include "parma_vtxSelector.h"
#include <apf.h>

namespace {
  typedef std::set<apf::MeshEntity*> SetEnt;

  class EdgeSelector : public parma::VtxSelector {
    public:
      EdgeSelector(apf::Mesh* m, apf::MeshTag* w)
        : VtxSelector(m, w) {}
    protected:
      // if all the elements bounded by an edge are in the plan then it is being
      // sent and its weight is counted
      bool sent(apf::Migration* plan, apf::MeshEntity* e) {
        apf::Adjacent elms;
        mesh->getAdjacent(e, mesh->getDimension(), elms);
        APF_ITERATE(apf::Adjacent, elms, elm)
          if( ! plan->has(*elm) ) 
            return false;
        return true;
      }

      double cavityWeight(apf::Migration* plan, SetEnt& s) {
        double w = 0;
        APF_ITERATE(SetEnt, s, sItr)
          if( sent(plan,*sItr) )
            w += getWeight(*sItr);
        return w;
      }

      void addCavityEdge(apf::MeshEntity* e, SetEnt& s) {
        apf::Adjacent adjEdge;
        mesh->getAdjacent(e, 1, adjEdge);
        APF_ITERATE(apf::Adjacent, adjEdge, adjItr)
          s.insert(*adjItr);
      }

      virtual double add(apf::MeshEntity*, apf::Up& cavity, const int destPid,
          apf::Migration* plan) {
        SetEnt s;
        for(int i=0; i < cavity.n; i++) {
          addCavityEdge(cavity.e[i], s); 
          plan->send(cavity.e[i], destPid);
        }
        return cavityWeight(plan, s);
      }
  };
}

namespace parma {
  Selector* makeEdgeSelector(apf::Mesh* m, apf::MeshTag* w) {
    return new EdgeSelector(m, w);
  }
}
