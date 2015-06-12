#include <PCU.h>
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
namespace parma {
  class ElmLtVtxEdge : public Targets {
    public:
      ElmLtVtxEdge(Sides* s, Weights* w[3], int sideTol,
          double vtxTol, double edgeTol, double alpha) {
        init(s, w, sideTol, vtxTol, edgeTol, alpha);
      }
      double total() {
        return totW;
      }
    private:
      ElmLtVtxEdge();
      double totW;
      void init(Sides* s, Weights* w[3], int sideTol,
          double vtxTol, double edgeTol, double alpha) {
        totW = 0;
        const Sides::Item* side;
        s->begin();
        while( (side = s->iterate()) ) {
          const int peer = side->first;
          const double peerVtxW = w[0]->get(peer);
          const double peerEdgeW = w[1]->get(peer);
          const double selfElmW = w[2]->self();
          const double peerElmW = w[2]->get(peer);
          const int peerSides = s->get(peer);
          if( selfElmW > peerElmW  && 
              peerVtxW < vtxTol && 
              peerEdgeW < edgeTol && 
              peerSides < sideTol ) {
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
  Targets* makeElmLtVtxEdgeTargets(Sides* s, Weights* w[3], int sideTol,
      double vtxTol, double edgeTol, double alpha) {
    return new ElmLtVtxEdge(s, w, sideTol, vtxTol, edgeTol, alpha);
  }
}
