#include "parma_selector.h"
#include "parma_targets.h"
#include "parma_weights.h"

namespace parma {

  class VtxSelector : public Selector {
    public:
      VtxSelector(apf::Mesh* m, apf::MeshTag* w)
        : Selector(m, w) {}
      apf::Migration* run(Targets* tgts) {
        apf::Migration* plan = new apf::Migration(mesh);
        vtag = mesh->createIntTag("selector_visited",1);
        const size_t maxBoundedElm = 6;
        double planW=0;
        for( size_t maxAdjElm=2; maxAdjElm<=maxBoundedElm; maxAdjElm+=2)
          planW += select(tgts, planW, maxAdjElm, plan);
        apf::removeTagFromDimension(mesh,vtag,0);
        mesh->destroyTag(vtag);
        return plan;
      }
    private:
      apf::MeshTag* vtag;
      VtxSelector();
      double add(apf::MeshEntity* vtx, const size_t maxAdjElm, 
          const int destPid, apf::Migration* plan) {
        apf::DynamicArray<apf::MeshEntity*> adjElms;
        mesh->getAdjacent(vtx, mesh->getDimension(), adjElms);
        if( adjElms.getSize() > maxAdjElm ) 
          return 0;
        for(size_t i=0; i<adjElms.getSize(); i++) {
          apf::MeshEntity* elm = adjElms[i];
          if ( mesh->hasTag(elm, vtag) ) continue;
          mesh->setIntTag(elm, vtag, &destPid); 
          plan->send(elm, destPid);
        }
        return getEntWeight(mesh,vtx,wtag);
      }
      double select(Targets* tgts, const double planW, 
          const size_t maxAdjElm, 
          apf::Migration* plan) {
        double planWeight = 0;
        apf::MeshEntity* vtx;
        apf::MeshIterator* itr = mesh->begin(0);
        while( (vtx = mesh->iterate(itr)) ) {
          if ( planW + planWeight > tgts->total() ) break;
          apf::Copies rmt;
          mesh->getRemotes(vtx, rmt);
          if( 1 == rmt.size() ) {
            int destPid = (rmt.begin())->first;
            if( tgts->has(destPid) )
              planWeight += add(vtx, maxAdjElm, destPid, plan);
          }
        }
        mesh->end(itr);
        return planWeight;
      }
  };
  Selector* makeVtxSelector(apf::Mesh* m, apf::MeshTag* w) {
    return new VtxSelector(m, w);
  }
} //end namespace parma
