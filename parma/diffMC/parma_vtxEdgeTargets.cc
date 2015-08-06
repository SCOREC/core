#include <PCU.h>
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
namespace parma {
  class VtxEdgeTargets : public Targets {
    public:
      VtxEdgeTargets(Sides* s, Weights* w[2], int sideTol, double vtxTol,
          double alpha) {
        init(s, w, sideTol, vtxTol, alpha);
      }
      double total() {
        return totW;
      }
    private:
      VtxEdgeTargets();
      double totW;
      void init(Sides* s, Weights* w[2], int sideTol, double vtxTol,
          double alpha) {
        totW = 0;
        const double selfEdgeW = w[1]->self();
        const Sides::Item* side;
        s->begin();
        while( (side = s->iterate()) ) {
          const int peer = side->first;
          const double peerVtxW = w[0]->get(peer);
          const double peerEdgeW = w[1]->get(peer);
          const int peerSides = s->get(peer);
          if( selfEdgeW > peerEdgeW &&
              peerVtxW < vtxTol &&
              peerSides < sideTol ) {
            const double difference = selfEdgeW - peerEdgeW;
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
  Targets* makeVtxEdgeTargets(Sides* s, Weights* w[2], int sideTol,
      double vtxTol, double alpha) {
    return new VtxEdgeTargets(s, w, sideTol, vtxTol, alpha);
  }
} //end namespace
