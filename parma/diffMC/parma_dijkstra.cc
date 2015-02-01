#include <limits.h>
#include <apf.h>
#include <PCU.h>
#include "parma_dijkstra.h" 
#include "parma_meshaux.h" 
#include "parma_distQ.h"

namespace {
  bool disconnected(apf::Mesh*m, apf::MeshTag* conn, apf::MeshEntity* e) {
    int c;
    apf::Adjacent adjElms;
    m->getAdjacent(e, m->getDimension(), adjElms);
    APF_ITERATE(apf::Adjacent, adjElms, adjItr) {
      m->getIntTag(*adjItr,conn,&c);
      if( c ) return false;
    }
    return true;
  }
}

namespace parma {
  void dijkstra(apf::Mesh* m, apf::MeshTag* c, apf::MeshTag* d,
      apf::MeshEntity* src) {
    parma::DistanceQueue<parma::Less> pq(m);
    int zero = 0;
    m->setIntTag(src, d, &zero);
    pq.push(src,0);

    while( !pq.empty() ) {
      apf::MeshEntity* v = pq.pop();
      int vd; m->getIntTag(v, d, &vd);
      if( vd == INT_MAX) continue;
      apf::Adjacent adjVtx;
      getEdgeAdjVtx(m,v,adjVtx);
      APF_ITERATE(apf::Adjacent, adjVtx, eItr) {
        apf::MeshEntity* u = *eItr;
        int ud; m->getIntTag(u,d,&ud);
        if( vd+1 < ud ) {
          int l = disconnected(m,c,u) ? INT_MAX : vd+1;
          m->setIntTag(u,d,&l);
          pq.push(u,l);
        }
      }
    }
  }

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
