#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
namespace parma {
  class VtxElmTargets : public Targets {
    public:
      VtxElmTargets(Sides* s, Weights* w[2], double alpha) {
        init(s, w, alpha);
      }
      double total() {
        return totW;
      }
    private:
      VtxElmTargets();
      double totW;
      void init(Sides* s, Weights* w[2], double alpha) {
        totW = 0;
        const Sides::Item* side;
        s->begin();
        while( (side = s->iterate()) ) {
          const int peer = side->first;
          const double selfVtxW = w[0]->self();
          const double selfElmW = w[1]->self();
          const double peerVtxW = w[0]->get(peer);
          const double peerElmW = w[1]->get(peer);
          if( selfVtxW > peerVtxW && selfElmW > peerElmW ) {
            const double difference = selfElmW - peerElmW;
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
  Targets* makeVtxElmTargets(Sides* s, Weights* w[2], double alpha) {
    return new VtxElmTargets(s,w,alpha);
  }
}
