#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
namespace parma {
  class VtxEdgeTargets : public Targets {
    public:
      VtxEdgeTargets(Sides* s, Weights* w[2], double alpha) {
        init(s, w, alpha);
      }
      double total() {
        return totW;
      }
    private:
      VtxEdgeTargets();
      double totW;
      void init(Sides* s, Weights* w[2], double alpha) {
        totW = 0;
        const Sides::Item* side;
        s->begin();
        while( (side = s->iterate()) ) {
          /*
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
          */
        }
        s->end();
      }
  };
  Targets* makeVtxEdgeTargets(Sides* s, Weights* w[2], double alpha) {
    return new VtxEdgeTargets(s,w,alpha);
  }

  class VtxEdgeElmTargets : public Targets {
    public:
      VtxEdgeElmTargets(Sides* s, Weights* w[3], double alpha) {
        init(s, w, alpha);
      }
      double total() {
        return totW;
      }
    private:
      VtxEdgeElmTargets();
      double totW;
      void init(Sides* s, Weights* w[3], double alpha) {
        totW = 0;
        const Sides::Item* side;
        s->begin();
        while( (side = s->iterate()) ) {
          /*
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
          */
        }
        s->end();
      }
  };
  Targets* makeVtxEdgeElmTargets(Sides* s, Weights* w[3], double alpha) {
    return new VtxEdgeElmTargets(s,w,alpha);
  }

} //end namespace

