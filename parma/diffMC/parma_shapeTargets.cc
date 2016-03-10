#include <PCU.h>
#include <parma.h>
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
#include "parma_convert.h"
#include "parma_commons.h"
#include <apf.h>
#include <string>

namespace parma {
  using parmaCommons::status;

  class ShapeTargets : public Targets {
    public:
      ShapeTargets(Sides* s, int small) {
        init(s,small);
        totW = 0;
      }
      double total() {
        return totW;
      }
    private:
      ShapeTargets();
      double totW;
      void init(Sides* s, int small) {
        PCU_Debug_Print("small %d\n", small);
        std::string sstr = s->print("sides");
        PCU_Debug_Print("%s\n", sstr.c_str());
        s->begin();
        const Sides::Item* side;
        while( (side = s->iterate()) ) {
          PCU_Debug_Print("side %d size %d\n", side->first, side->second);
          if( side->second <= small ) {
            PCU_Debug_Print("adding %d to targets\n", side->first);
            set(side->first, small);
          }
        }
        s->end();
        std::string tgtstr = print("targets");
        PCU_Debug_Print("%s\n", tgtstr.c_str());
      }
  };
  Targets* makeShapeTargets(Sides* s, int small) {
    return new ShapeTargets(s,small);
  }
} //end namespace
