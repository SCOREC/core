#include <pcu_util.h>
#include <PCU.h>
#include "parma_entWeights.h"
#include "parma_sides.h"

namespace parma {  
  double getMaxWeight(apf::Mesh* m, apf::MeshTag* w, int entDim) {
    double locW = getWeight(m,w,entDim);
    return PCU_Max_Double(locW);
  }

  double getAvgWeight(apf::Mesh* m, apf::MeshTag* w, int entDim) {
    double locW = getWeight(m,w,entDim);
    return PCU_Add_Double(locW) / PCU_Comm_Peers();
  }

  double getWeight(apf::Mesh* m, apf::MeshTag* w, int entDim) {
    PCU_ALWAYS_ASSERT(entDim >= 0 && entDim <= 3);
    apf::MeshIterator* it = m->begin(entDim);
    apf::MeshEntity* e;
    double sum = 0;
    while ((e = m->iterate(it)))
      sum += parma::getEntWeight(m, e, w);
    m->end(it);
    return sum;
  }

  void getImbalance(Weights* w, double& imb, double& avg) {
    double sum, max;
    sum = max = w->self();
    sum = PCU_Add_Double(sum);
    max = PCU_Max_Double(max);
    avg = sum/PCU_Comm_Peers();
    imb = max/avg;
  }

  double getMaxWeight(Weights* w) {
    return PCU_Max_Double(w->self());
  }

  double getEntWeight(apf::Mesh* m, apf::MeshEntity* e, apf::MeshTag* w) {
    PCU_ALWAYS_ASSERT(m->hasTag(e,w));
    double weight;
    m->getDoubleTag(e,w,&weight);
    return weight;
  }

  EntWeights::EntWeights(apf::Mesh* m, apf::MeshTag* w, Sides* s, int d) 
    : Weights(m, w, s), entDim(d) 
  {
    PCU_ALWAYS_ASSERT(entDim >= 0 && entDim <= 3);
    weight = getWeight(m, w, entDim);
    init(m, w, s);
  }
  double EntWeights::self() {
    return weight;
  }
  double EntWeights::getEntWeight(apf::Mesh* m, apf::MeshEntity* e, 
      apf::MeshTag* w) {
    PCU_ALWAYS_ASSERT(m->hasTag(e,w));
    double entW = 0;
    m->getDoubleTag(e,w,&entW);
    return entW;
  }

  void EntWeights::init(apf::Mesh*, apf::MeshTag*, Sides* s) {
    PCU_Comm_Begin();
    const Sides::Item* side;
    s->begin();
    while( (side = s->iterate()) ) 
      PCU_COMM_PACK(side->first, weight);
    s->end();
    PCU_Comm_Send();
    while (PCU_Comm_Listen()) {
      double otherWeight;
      PCU_COMM_UNPACK(otherWeight);
      set(PCU_Comm_Sender(), otherWeight);
    }
  }
  Weights* makeEntWeights(apf::Mesh* m, apf::MeshTag* w, Sides* s, int dim) {
    return new EntWeights(m, w, s, dim);
  }


} //end namespace
