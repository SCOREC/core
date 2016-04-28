#include <PCU.h>
#include <parma.h>
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
#include "parma_convert.h"
#include "parma_commons.h"
#include <apf.h>
#include <string>
#include <limits.h>

namespace parma {
  using parmaCommons::status;

  class ShapeTargets : public Targets {
    public:
      ShapeTargets(Sides* s) {
        smallLimit = 10;
        init(s);
        totW = 0;
      }
      double total() {
        return totW;
      }
    private:
      ShapeTargets();
      int smallLimit;
      double totW;
      void init(Sides* s) {
        const unsigned maxNb = TO_UINT(PCU_Max_Int(s->size()));
        if( s->size() != maxNb ) return;
        PCU_Debug_Print("maxNb %d\n", maxNb);
        std::string sstr = s->print("sides");
        PCU_Debug_Print("%s\n", sstr.c_str());
        int small = INT_MAX;
        s->begin();
        const Sides::Item* side;
        while( (side = s->iterate()) )
          if( side->second < small )
            small = side->second;
        s->end();
        PCU_Debug_Print("small %d\n", small);
        if( small > smallLimit ) return;
        s->begin();
        while( (side = s->iterate()) )
          if( side->second <= small )
            set(side->first, small);
        s->end();
        std::string tgtstr = print("targets");
        PCU_Debug_Print("%s\n", tgtstr.c_str());
      }
  };
  Targets* makeShapeTargets(Sides* s) {
    return new ShapeTargets(s);
  }
} //end namespace
