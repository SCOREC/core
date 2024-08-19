#include <pcu_util.h>
#include "parma_entWeights.h"
#include "parma_sides.h"

namespace parma {  
  double getMaxWeight(apf::Mesh* m, apf::MeshTag* w, int entDim) {
    double locW = getWeight(m,w,entDim);
    return m->getPCU()->Max<double>(locW);
  }

  double getAvgWeight(apf::Mesh* m, apf::MeshTag* w, int entDim) {
    double locW = getWeight(m,w,entDim);
    return m->getPCU()->Add<double>(locW) / m->getPCU()->Peers();
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

  void getImbalance(Weights* w, double& imb, double& avg, pcu::PCU *PCUObj) {
    double sum, max;
    sum = max = w->self();
    sum = PCUObj->Add<double>(sum);
    max = PCUObj->Max<double>(max);
    avg = sum/PCUObj->Peers();
    imb = max/avg;
  }

  double getMaxWeight(Weights* w, pcu::PCU *PCUObj) {
    return PCUObj->Max<double>(w->self());
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

  void EntWeights::init(apf::Mesh* m, apf::MeshTag*, Sides* s) {
    m->getPCU()->Begin();
    const Sides::Item* side;
    s->begin();
    while( (side = s->iterate()) ) 
      m->getPCU()->Pack(side->first, weight);
    s->end();
    m->getPCU()->Send();
    while (m->getPCU()->Listen()) {
      double otherWeight;
      m->getPCU()->Unpack(otherWeight);
      set(m->getPCU()->Sender(), otherWeight);
    }
  }
  Weights* makeEntWeights(apf::Mesh* m, apf::MeshTag* w, Sides* s, int dim) {
    return new EntWeights(m, w, s, dim);
  }


} //end namespace
