#include <assert.h>
#include <PCU.h>
#include "parma_surfToVol.h"
#include "parma_sides.h"

namespace parma {  
  class SidesToElements : public SurfToVol {
    public:
      SidesToElements(apf::Mesh* m, Sides* s) {
        const int surf = s->total();
        const int vol = m->count(m->getDimension());
        avgSurfToVol = surfToVolImb = surf/(double)vol;
        PCU_Add_Doubles(&avgSurfToVol, 1);
        avgSurfToVol /= PCU_Comm_Peers();
        surfToVolImb /= avgSurfToVol;
        init(m,s);
      }
    private:
      double avgSurfToVol;
      void init(apf::Mesh*, Sides* s) {
        PCU_Comm_Begin();
        const Sides::Item* side;
        s->begin();
        while( (side = s->iterate()) ) 
          PCU_COMM_PACK(side->first, surfToVolImb);
        s->end();
        PCU_Comm_Send();
        while (PCU_Comm_Listen()) {
          double otherSurfToVolImb;
          PCU_COMM_UNPACK(otherSurfToVolImb);
          set(PCU_Comm_Sender(), otherSurfToVolImb);
        }
      }
  };
  SurfToVol* makeSidesToElements(apf::Mesh* m, Sides* s) {
    return new SidesToElements(m,s);
  }
} //end namespace
