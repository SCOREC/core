#include <PCU.h>
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
namespace parma {
  class WeightSideTargets : public Targets {
    public:
      WeightSideTargets(Sides* s, Weights* w, int sideTol, double alpha) {
        init(s, w, sideTol, alpha);
      }
      double total() {
        return totW;
      }
    private:
      WeightSideTargets();
      double totW;
      void init(Sides* s, Weights* w, int sideTol, double alpha) {
        const double selfW = w->self();
        totW = 0;
        const Sides::Item* side;
        s->begin();
        while( (side = s->iterate()) ) {
          const int peer = side->first;
          const double peerW = w->get(peer);
          const int peerSides = s->get(peer);
          if( selfW > peerW && 
              peerSides < sideTol ) {
            const double difference = selfW - peerW;
            double sideFraction = side->second;
            sideFraction /= s->total();
            double scaledW = difference * sideFraction * alpha;
            set(peer, scaledW);
            totW+=scaledW;
          }
        }
        s->end();
      }
  };
  Targets* makeWeightSideTargets(Sides* s, Weights* w, int sideTol, 
      double alpha) {
    return new WeightSideTargets(s, w, sideTol, alpha);
  }
}
