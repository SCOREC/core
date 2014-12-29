#include <PCU.h>
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
#include "parma_surfToVol.h"
namespace parma {
  class WeightSvTargets : public Targets {
    public:
      WeightSvTargets(Sides* s, Weights* w, SurfToVol* sv, 
          double aspectTol, double alpha) {
        init(s, w, sv, aspectTol, alpha);
      }
      double total() {
        return totW;
      }
    private:
      WeightSvTargets();
      double totW;
      void init(Sides* s, Weights* w, SurfToVol* sv, 
          double aspectTol, double alpha) {
        const double selfElmW = w->self();
        totW = 0;
        const Sides::Item* side;
        s->begin();
        while( (side = s->iterate()) ) {
          const int peer = side->first;
          const double peerElmW = w->get(peer);
          const double peerAspect = sv->get(peer);
          if( selfElmW > peerElmW && 
              peerAspect < aspectTol ) {
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
  Targets* makeWeightSvTargets(Sides* s, Weights* w, SurfToVol* sv, 
      double aspectTol, double alpha) {
    return new WeightSvTargets(s, w, sv, aspectTol, alpha);
  }
}
