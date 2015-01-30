#include <PCU.h>
#include <parma.h>
#include "parma_step.h"
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
#include "parma_selector.h"
#include "parma_stop.h"

namespace parma {
  Stepper::Stepper(apf::Mesh* mIn, double alphaIn,
     Sides* s, Weights* w, Targets* t, Selector* sel,
     Stop* stopper) 
    : m(mIn), alpha(alphaIn), sides(s), weights(w), targets(t), 
    selects(sel), stop(stopper) {
  }

  Stepper::~Stepper() {
    delete sides;
    delete weights;
    delete targets;
    delete selects;
    delete stop;
  }

  bool Stepper::step(double maxImb, int verbosity) {
    const double imb = imbalance();
    if ( !PCU_Comm_Self() && verbosity )
      fprintf(stdout, "imbalance %.3f\n", imb);
    if ( stop->stop(imb,maxImb) )
      return false;
    apf::Migration* plan = selects->run(targets);
    const double t0 = PCU_Time();
    m->migrate(plan);
    if ( !PCU_Comm_Self() && verbosity )
      fprintf(stdout, "elements migrated in %f seconds\n", PCU_Time()-t0);
    if( verbosity > 1 ) 
      Parma_PrintPtnStats(m, "endStep", (verbosity>2));
    return true;
  }

  double Stepper::imbalance() { 
    double maxWeight = 0, totalWeight = 0;
    maxWeight = totalWeight = weights->self();
    PCU_Add_Doubles(&totalWeight,1);
    PCU_Max_Doubles(&maxWeight, 1);
    double averageWeight = totalWeight / PCU_Comm_Peers();
    return maxWeight / averageWeight;
  }
}
