#include "parma_targets.h"
#include "parma_sides.h"
#include "parma_targets.h"
#include "parma_ghosts.h"
namespace parma {
  class GhostTargets : public Associative<double> {
    public:
      GhostTargets(Sides* s, Weights* w, Ghosts* g, double alpha) 
        : Targets(s, w, alpha)
      {
        init(s, w, g, alpha);
      }
      double total() {
        return totW;
      }
    private:
      GhostTargets();
      double totW;
      void init(Sides* s, Weights* w, Ghosts* g, double alpha) {
        totW = 0;
        const Sides::Item* side;
        s->begin();
        while( (side = s->iterate()) ) {
          const int peer = side->first;
          const double selfW = w->self() + g->self();
          const double peerW = g->get(peer) + w->get(peer);
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
  Targets* makeGhostTargets(Sides* s, Weights* w, Ghosts* g, double alpha) {
    return new GhostTargets(s,w,g,alpha);
  }
} //end namespace

