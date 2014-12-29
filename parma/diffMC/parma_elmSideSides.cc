#include "parma_sides.h"
#include <apf.h>

namespace parma {  
  class ElmSideSides : public Sides {
    public:
      ElmSideSides(apf::Mesh* m) : Sides(m) {
        init(m);
      }
    private:
      void init(apf::Mesh* m) {
        apf::MeshEntity* s;
        apf::MeshIterator* it = m->begin(m->getDimension()-2);
        totalSides = 0;
        while ((s = m->iterate(it)))
          if (m->isShared(s)) {
            apf::Copies rmts;
            m->getRemotes(s, rmts);
            APF_ITERATE(apf::Copies, rmts, r)
              set(r->first, get(r->first)+1);
            ++totalSides;
          }
        m->end(it);
      }
  };

  Sides* makeElmSideSides(apf::Mesh* m) {
    return new ElmSideSides(m);
  }
}


