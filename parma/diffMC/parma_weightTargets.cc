#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
namespace parma {
  class WeightTargets : public Targets {
    public:
      WeightTargets(Sides* s, Weights* w, double alpha) {
        init(s, w, alpha);
      }
      double total() {
        return totW;
      }
    private:
      WeightTargets();
      double totW;
      void init(Sides* s, Weights* w, double alpha) {
        totW = 0;
        const Sides::Item* side;
        s->begin();
        while( (side = s->iterate()) ) {
          const int peer = side->first;
          const double selfW = w->self();
          const double peerW = w->get(peer);
          if ( selfW > peerW ) {
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
  Targets* makeTargets(Sides* s, Weights* w, double alpha) {
    return new WeightTargets(s,w,alpha);
  }
} //end namespace

