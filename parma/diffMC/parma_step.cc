#include <PCU.h>
#include <parma.h>
#include "parma_step.h"
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
#include "parma_selector.h"
#include "parma_stop.h"
#include "parma_commons.h"

namespace parma {
  using parmaCommons::status;

  Stepper::Stepper(apf::Mesh* mIn, double alphaIn,
     Sides* s, Weights* w, Targets* t, Selector* sel,
     const char* entType, Stop* stopper)
    : m(mIn), alpha(alphaIn), sides(s), weights(w), targets(t),
    selects(sel), name(entType), stop(stopper) {
      verbose = 0;
  }

  Stepper::~Stepper() {
    delete sides;
    delete weights;
    delete targets;
    delete selects;
    delete stop;
  }

  bool Stepper::step(double maxImb, int verbosity) {
    double imb, avg;
    getImbalance(weights, imb, avg);
    if ( !PCU_Comm_Self() && verbosity )
      status("%s imbalance %.3f avg %.3f\n", name, imb, avg);
    if ( stop->stop(imb,maxImb) )
      return false;
    apf::Migration* plan = selects->run(targets);
    int planSz = PCU_Add_Int(plan->count());
    const double t0 = PCU_Time();
    m->migrate(plan);
    if ( !PCU_Comm_Self() && verbosity )
      status("%d elements migrated in %f seconds\n", planSz, PCU_Time()-t0);
    if( verbosity > 1 ) 
      Parma_PrintPtnStats(m, "endStep", (verbosity>2));
    return true;
  }
}
