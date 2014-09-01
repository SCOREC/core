#include <PCU.h>
#include "parma_base.h"
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
#include "parma_selector.h"

namespace parma {
  Balancer::Balancer(apf::Mesh* mIn, apf::MeshTag* wIn, double alphaIn,
     Sides* s, Weights* w, Targets* t, Selector* sel,
     bool (*fn)(double imb, double maxImb))
    : m(mIn), w(wIn), alpha(alphaIn), 
      sides(s), weights(w), targets(t), selects(sel), stop(fn) {}

  Balancer::~Balancer() {
    delete sides;
    delete weights;
    delete targets;
    delete selects;
  }

  bool Balancer::run(double maxImb, int verbosity) {
    const double imb = imbalance();
    if ( !PCU_Comm_Self() && verbosity )
      fprintf(stdout, "imbalance %.3f\n", imb);
    if ( stop(imb,maxImb) )
      return false;
    apf::Migration* plan = selects->run(targets);
    PCU_Debug_Print("plan size %d\n", plan->count());
    m->migrate(plan);
    return true;
  }

  double Balancer::imbalance() { 
    double maxWeight = 0, totalWeight = 0;
    maxWeight = totalWeight = weights->self();
    PCU_Add_Doubles(&totalWeight,1);
    PCU_Max_Doubles(&maxWeight, 1);
    double averageWeight = totalWeight / PCU_Comm_Peers();
    return maxWeight / averageWeight;
  }
}
