#include "parma_bdryVtx.h"
#include "parma_distQ.h"

namespace {
  class DistItr : public parma::BdryVtxItr {
    public:
      DistItr(apf::Mesh* m, apf::MeshTag* d) 
      {
        dq = new parma::DistanceQueue<parma::Greater>(m);

        apf::MeshEntity* v;
        apf::MeshIterator* it = m->begin(0);
        while( (v = m->iterate(it)) )
          if( m->isShared(v) ) {
            int dist; m->getIntTag(v,d,&dist);
            dq->push(v, dist);
          }
        m->end(it);
      }

      ~DistItr() 
      {
        delete dq;
      }

      apf::MeshEntity* next() 
      {
        if( dq->empty() )
          return NULL;
        else
          return dq->pop();
      }
    private:
      parma::DistanceQueue<parma::Greater>* dq;
  };
}

namespace parma {
  BdryVtxItr* makeBdryVtxDistItr(apf::Mesh* m, apf::MeshTag* d) {
    return new DistItr(m, d);
  }
}
