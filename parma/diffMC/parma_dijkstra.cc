#include <limits.h>
#include <apf.h>
#include <PCU.h>
#include "parma_dijkstra.h" 
#include "parma_meshaux.h" 
#include "parma_distQ.h"
#include "parma_convert.h"

/*
 * Components can have a seperator that may already be distanced.  If the
 * seperator vertex distances are lower than the distance the current component
 * would set than the vertices on the other side of the seperator belonging
 * exclusively to the current component will not be distanced.
 */
namespace {
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

  void getCavity(apf::Mesh* m, apf::MeshEntity* v, apf::Up& cavity) {
    cavity.n = 0;
    apf::Adjacent elms;
    m->getAdjacent(v, m->getDimension(), elms);
    APF_ITERATE(apf::Adjacent, elms, adjItr)
      cavity.e[(cavity.n)++] = *adjItr;
  }

  typedef std::set<apf::MeshEntity*> CavEnts;
  CavEnts* getCavityEdges(apf::Mesh* m, parma::DijkstraContains* c, 
      apf::Up& cavity) {
    CavEnts* ce = new CavEnts;
    for(int i=0; i<cavity.n; i++) {
      apf::Downward edges;
      int ne = m->getDownward(cavity.e[i], 1, edges);
      for(int j=0; j<ne; j++) {
        apf::Downward verts;
        int nv = m->getDownward(edges[j], 0, verts);
        assert(nv==2);
        if( c->has(verts[0]) && c->has(verts[1]) )
          ce->insert(edges[j]);
      }
    }
    return ce;
  }

  apf::MeshEntity* parentVtx(apf::Mesh* m, parma::DijkstraContains* c,
      apf::MeshTag* d, apf::MeshEntity* u) {
    apf::MeshEntity* parent = NULL;
    int ud; m->getIntTag(u,d,&ud);
    apf::Adjacent adjVtx;
    getEdgeAdjVtx(m,u,adjVtx);
    APF_ITERATE(apf::Adjacent, adjVtx, v) {
      int vd; m->getIntTag(*v,d,&vd);
      if( vd == ud-1 && c->has(*v) ) {
        parent = *v;
        break;
      }
    }
    return parent;
  }

  typedef std::set<apf::MeshEntity*> Level;
  void walkCavEdges(apf::Mesh* m, CavEnts* ce, apf::MeshEntity* p,
      apf::MeshEntity* s, apf::Adjacent& adjVtx) {
    typedef apf::DynamicArray<apf::MeshEntity*> entArr;
    //should be an upper bound
    const unsigned maxVisited = TO_UINT(ce->size()*2);
    entArr visited(maxVisited);
    adjVtx.setSize(maxVisited);
    unsigned visCnt=0;
    size_t numAdj=0;
    apf::MeshTag* lvlT = m->createIntTag("parmaCavWalk",1);
    Level cur;
    Level next; 
    next.insert(p);
    int treeDepth = 0;
    while( ! next.empty() ) {
      cur = next;
      next.clear();
      treeDepth++;
      APF_ITERATE(Level, cur, vtxItr) {
        apf::MeshEntity* u = *vtxItr;
        if( m->hasTag(u,lvlT) ) continue;
        m->setIntTag(u,lvlT,&treeDepth);
        visited[visCnt++] = u;
        assert(visCnt < maxVisited);
        apf::Up edges;
        m->getUp(u,edges);
        for(int i=0; i<edges.n; i++) {
          if( ! ce->count(edges.e[i]) ) 
            continue; //not a cavity edge
          apf::MeshEntity* v = 
            apf::getEdgeVertOppositeVert(m,edges.e[i],u);
          if( v != s && ! m->hasTag(v,lvlT) && ! cur.count(v) ) {
            next.insert(v);
            adjVtx[numAdj++] = v;
          }
        }
      }
    }
    adjVtx.setSize(numAdj);
    for(std::size_t i = 0; i<visCnt; i++) {
      assert(m->hasTag(visited[i],lvlT));
      m->removeTag(visited[i],lvlT);
    }
    m->destroyTag(lvlT);
  }

  void getConnectedVtx(apf::Mesh* m, parma::DijkstraContains* c, 
      apf::MeshTag* d, apf::MeshEntity* u, apf::Adjacent& adjVtx) {
    if( c->bdryHas(u) ) {
      apf::MeshEntity* p = parentVtx(m,c,d,u);
      if( p ) {
        apf::Up cav;
        getCavity(m,u,cav);
        CavEnts* ce = getCavityEdges(m,c,cav);
        walkCavEdges(m,ce,p,u,adjVtx);
        delete ce;
      }
    }
    if( !adjVtx.getSize() )
      getEdgeAdjVtx(m,u,adjVtx); 
  }

  void dijkstra_(apf::Mesh* m, parma::DijkstraContains* c,
      parma::DistanceQueue<parma::Less>& pq, apf::MeshTag* d) {
    while( !pq.empty() ) {
      apf::MeshEntity* v = pq.pop();
      if( ! c->has(v) ) continue;
      int vd; m->getIntTag(v, d, &vd);
      assert( vd >= 0 && vd != INT_MAX );
      apf::Adjacent adjVtx;
      getConnectedVtx(m,c,d,v,adjVtx);
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

namespace parma {
  void dijkstra(apf::Mesh* m, DijkstraContains* c,
      apf::MeshEntity* src, apf::MeshTag* d) {
    resetDist(m,c,d);
    parma::DistanceQueue<parma::Less> pq(m);
    int zero = 0;
    m->setIntTag(src, d, &zero);
    pq.push(src,0);
    dijkstra_(m,c,pq,d);
  }
  void dijkstra(apf::Mesh* m, DijkstraContains* c,
      DistanceQueue<Less>& pq, apf::MeshTag* d) {
    dijkstra_(m,c,pq,d);
  }
}
