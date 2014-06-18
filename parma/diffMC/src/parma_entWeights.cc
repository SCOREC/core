#include "parma_weights.h"
#include <PCU.h>
#include <assert.h>

namespace parma {  
  class EntWeights : public Weights {
    public:
      EntWeights(apf::Mesh* m, apf::MeshTag* w, Sides* s, int d) 
        : Weights(m, w, s), entDim(d) 
      {
        assert(entDim >= 0 && entDim <= 3);
        weight = selfWeight(m, w);
        init(m, w, s);
      }

      double self() {
        return weight;
      }
    private:
      Weights();
      virtual double getEntWeight(apf::Mesh* m, apf::MeshEntity* e, apf::MeshTag* w) {
        assert(m->hasTag(e,w));
        double w = 0;
        m->getDoubleTag(e,w,&w);
        return w;
      }
      double selfWeight(apf::Mesh* m, apf::MeshTag* w) {
        apf::MeshIterator* it = m->begin(entDim);
        apf::MeshEntity* e;
        double sum = 0;
        while ((e = m->iterate(it)))
          sum += getEntWeight(m, e, w);
        m->end(it);
        return sum;
      }
      void init(apf::Mesh* m, apf::MeshTag* w, Sides* s) {
        PCU_Comm_Begin();
        const Sides::Item* side;
        s->begin();
        while( (side = s->iterate()) ) 
          PCU_COMM_PACK(side->first, selfWeight);
        s->end();
        PCU_Comm_Send();
        while (PCU_Comm_Listen()) {
          double otherWeight;
          PCU_COMM_UNPACK(otherWeight);
          set(PCU_Comm_Sender(), otherWeight);
        }
      }
  };

  Weights* makeEntWeights(apf::Mesh* m, apf::MeshTag* w, Sides* s, int dim) {
    return new EntWeights(m, w, s, dim);
  }
} //end namespace
