#include <limits.h>
#include <apf.h>
#include <PCU.h>
#include "parma_dijkstra.h" 
#include "parma_meshaux.h" 
#include "parma_distQ.h"

/*
 * Components can have a seperator that may already be distanced.  If the
 * seperator vertex distances are lower than the distance the current component
 * would set than the vertices on the other side of the seperator belonging
 * exclusively to the current component will not be distanced.
 */
void resetDist(apf::Mesh* m, parma::DijkstraContains* c, apf::MeshTag* d) {
  const int far = INT_MAX;
  apf::MeshEntity* v;
  apf::MeshIterator* it = m->begin(0);
  while( (v = m->iterate(it)) ) {
    if( !c->has(v) ) continue;
    m->setIntTag(v,d,&far);
  }
  m->end(it);
}

namespace parma {
  void dijkstra(apf::Mesh* m, DijkstraContains* c,
      apf::MeshEntity* src, apf::MeshTag* d) {
    resetDist(m,c,d);
    parma::DistanceQueue<parma::Less> pq(m);
    int zero = 0;
    m->setIntTag(src, d, &zero);
    pq.push(src,0);

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
        }
      }
    }
  }
}
