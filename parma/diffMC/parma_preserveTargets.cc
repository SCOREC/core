#include <PCU.h>
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
namespace parma {
  class PreserveTargets : public Targets {
    public:
      PreserveTargets(Sides* s, Weights* balance, Weights* preserve,
                int sideTol, double preserveTol, double alpha) {
        init(s, balance, preserve, sideTol, preserveTol, alpha);
      }
      double total() {
        return totW;
      }
    private:
      PreserveTargets();
      double totW;
      void init(Sides* s, Weights* balance, Weights* preserve,
          int sideTol, double preserveTol, double alpha) {
        totW = 0;
        const Sides::Item* side;
        s->begin();
        while( (side = s->iterate()) ) {
          const int peer = side->first;
          const double peerPresW = preserve->get(peer);
          const double selfBalW = balance->self();
          const double peerBalW = balance->get(peer);
          const int peerSides = s->get(peer);
          if( selfBalW > peerBalW  &&
              peerPresW < preserveTol &&
              peerSides < sideTol ) {
            const double difference = selfBalW - peerBalW;
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
  Targets* makePreservingTargets(Sides* s, Weights* balanceW, Weights* preserveW,
      int sideTol, double preserveTol, double alpha) {
    return new PreserveTargets(s, balanceW, preserveW, sideTol, preserveTol, alpha);
  }
}
