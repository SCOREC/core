#include <limits.h>
#include <apf.h>
#include <PCU.h>
#include "parma_dijkstra.h" 
#include "parma_meshaux.h" 
#include "parma_distQ.h"

namespace parma {
  void dijkstra(apf::Mesh* m, DijkstraContains* c,
      apf::MeshEntity* src, apf::MeshTag* d) {
    parma::DistanceQueue<parma::Less> pq(m);
    int zero = 0;
    m->setIntTag(src, d, &zero);
    pq.push(src,0);

    unsigned count=0;

    while( !pq.empty() ) {
      apf::MeshEntity* v = pq.pop();
      if( ! c->has(v) ) continue;
      int vd; m->getIntTag(v, d, &vd);
      assert( vd >= 0 && vd != INT_MAX );
      apf::Adjacent adjVtx;
      getEdgeAdjVtx(m,v,adjVtx);
      APF_ITERATE(apf::Adjacent, adjVtx, eItr) {
        apf::MeshEntity* u = *eItr;
        int ud; m->getIntTag(u,d,&ud);
        if( vd+1 < ud && c->has(u) ) {
          int l = vd+1;
          m->setIntTag(u,d,&l);
          pq.push(u,l);
          count++;
        }
      }
    }
    PCU_Debug_Print("count %u\n", count);
  }
}
