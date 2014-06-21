#include <PCU.h>
#include <assert.h>
#include "parma_entWeights.h"

namespace parma {  
  class OwnedVtxWeights : public EntWeights {
    public:
      OwnedVtxWeights(apf::Mesh* m, apf::MeshTag* w, Sides* s)
        : EntWeights(m, w, s, 0) {}
    private:
      OwnedVtxWeights();
      double getEntWeight(apf::Mesh* m, apf::MeshEntity* e, apf::MeshTag* w) {
        assert(m->hasTag(e,w));
        double entW = 0;
        if (m->isOwned(e))
          m->getDoubleTag(e,w,&entW);
        return entW;
      }
  };
  Weights* makeOwnedVtxWeights(apf::Mesh* m, apf::MeshTag* w, Sides* s) {
    return new OwnedVtxWeights(m, w, s);
  }
} //end namespace
