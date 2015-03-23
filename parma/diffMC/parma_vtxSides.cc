#include "parma_sides.h"
#include "parma_dcpart.h"
#include <apf.h>


namespace parma {  
  class VtxSides : public Sides {
    public:
      VtxSides(apf::Mesh* m) : Sides(m) {
        init(m);
      }
    private:
      void init(apf::Mesh* m) {
        apf::MeshEntity* s;
        apf::MeshIterator* it = m->begin(0);
        totalSides = 0;
	dcPart dc(m);
        while ((s = m->iterate(it))) {
	  apf::Adjacent adj;
          m->getAdjacent(s,m->getDimension(),adj);
          if (m->isShared(s)&&!dc.isIsolated(s)) {
            apf::Copies rmts;
            m->getRemotes(s, rmts);
            APF_ITERATE(apf::Copies, rmts, r)
              set(r->first, get(r->first)+1);
            ++totalSides;
          }
	}
        m->end(it);
      }
  };

  Sides* makeVtxSides(apf::Mesh* m) {
    return new VtxSides(m);
  }
}
