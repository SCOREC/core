#include <apf.h>
#include <PCU.h>
#include <list>
#include <set>
#include <limits.h>
#include "parma_graphDist.h"
#include "parma_dijkstra.h"
#include "parma_meshaux.h"
#include "parma_commons.h"

namespace {
  apf::MeshTag* initTag(apf::Mesh* m, const char* name,
      int initVal=0, int dim=0) {
    apf::MeshTag* t = m->createIntTag(name,1);
    apf::MeshEntity* e;
    apf::MeshIterator* it = m->begin(dim);
    while( (e = m->iterate(it)) )
      m->setIntTag(e,t,&initVal);
    m->end(it);
    return t;
  }

  typedef std::set<apf::MeshEntity*> Level;

  bool onBoundary(apf::Mesh* m, apf::MeshEntity* e) {
    int gd = m->getModelType(m->toModel(e));
    int md = m->getDimension();
    bool shared = m->isShared(e);
    return shared || (!shared && gd < md);
  }

  int walkInward(apf::MeshTag* n, apf::Mesh* m) {
    Level cur;
    Level next;
    apf::MeshEntity* v;
    apf::MeshIterator* it = m->begin(0);
    while( (v = m->iterate(it)) )
      if( onBoundary(m,v) )
        cur.insert(v);
    m->end(it);

    size_t count = 0;
    int depth = 1;
    while( count != m->count(0) ) {
      APF_ITERATE(Level, cur, vtxItr) {
        v = *vtxItr;
        int lvl; m->getIntTag(v,n,&lvl);
        if( lvl ) continue;
        m->setIntTag(v,n,&depth);
        count++;
        assert(count <= m->count(0));
        apf::Adjacent adjVtx;
        getEdgeAdjVtx(m,v,adjVtx);
        APF_ITERATE(apf::Adjacent, adjVtx, vItr) {
          m->getIntTag(*vItr,n,&lvl);
          if( !lvl )
            next.insert(*vItr);
        }
      }
      cur = next;
      assert(cur.size() == cur.size());
      next.clear();
      depth++;
    }
    return depth-1;
  }

  Level* reduce(apf::MeshTag* n, apf::Mesh* m, Level* l, int depth) {
    Level* g = new Level;
    while( !l->empty() ) {
      apf::MeshEntity* v = *(l->begin()); l->erase(v);
      g->insert(v);
      Level q;
      q.insert(v);
      while( !q.empty() ) {
        v = *(q.begin()); q.erase(v);
        apf::Adjacent adjVtx;
        getEdgeAdjVtx(m,v,adjVtx);
        APF_ITERATE(apf::Adjacent, adjVtx, u) {
          int d; m->getIntTag(*u,n,&d);
          if( d == depth && l->count(*u) ) {
            l->erase(*u);
            q.insert(*u);
          }
        }
      }
    }
    delete l;
    return g;
  }

  Level* getVerts(apf::MeshTag* n, apf::Mesh* m, int depth) {
    Level* l = new Level;
    apf::MeshEntity* v;
    apf::MeshIterator* it = m->begin(0);
    while( (v = m->iterate(it)) ) {
      int d; m->getIntTag(v,n,&d);
      if( d == depth )
        l->insert(v);
    }
    m->end(it);
    l = reduce(n,m,l,depth);
    return l;
  }

  Level* getCentralVerts(apf::Mesh* m) {
    double t0 = PCU_Time();
    apf::MeshTag* lvlsT = initTag(m, "parmaSelectorLevels");
    int depth = walkInward(lvlsT, m);
    Level* deepest = getVerts(lvlsT,m,depth);
    apf::removeTagFromDimension(m,lvlsT,0);
    m->destroyTag(lvlsT);
    parmaCommons::printElapsedTime("getCentralVerts", PCU_Time()-t0);
    return deepest;
  }

  Level* getCentralElms(apf::Mesh* m, Level& verts) {
    Level* elms = new Level;
    APF_ITERATE(Level, verts, itr) {
      apf::Adjacent adjElms;
      m->getAdjacent(*itr, m->getDimension(), adjElms);
      APF_ITERATE(apf::Adjacent, adjElms, adjItr)
        elms->insert(*adjItr);
    }
    return elms;
  }

  inline void getFaceAdjElms(apf::Mesh* m, apf::MeshEntity* e,
      apf::Adjacent& adj) {
    int bridge = m->getDimension()-1; int tgt = m->getDimension();
    getBridgeAdjacent(m, e, bridge, tgt, adj);
  }

  void walkElms(apf::Mesh* m, apf::MeshTag* conn, apf::MeshEntity* src) {
    int one = 1;
    int count = 0;
    std::list<apf::MeshEntity*> elms;
    elms.push_back(src);
    while( !elms.empty() ) {
      apf::MeshEntity* e = elms.front();
      elms.pop_front();
      int c; m->getIntTag(e,conn,&c);
      if( c ) continue;
      m->setIntTag(e,conn,&one);
      count++;
      apf::Adjacent adjElms;
      getFaceAdjElms(m,e,adjElms);
      APF_ITERATE(apf::Adjacent, adjElms, eItr) {
        m->getIntTag(*eItr,conn,&c);
        if( !c )
          elms.push_back(*eItr);
      }
    }
  }

  apf::MeshTag* computeDistance(apf::Mesh* m, Level& verts) {
    double t0 = PCU_Time();
    Level* centralElms = getCentralElms(m, verts);
    int initVal = 0;
    apf::MeshTag* connT =
      initTag(m, "parmaElmConnectivity", initVal, m->getDimension());
    APF_ITERATE(Level, *centralElms, itr)
      walkElms(m,connT,*itr);
    delete centralElms;
    parmaCommons::printElapsedTime("computeConnectivity", PCU_Time()-t0);

    t0 = PCU_Time();
    initVal = INT_MAX;
    apf::MeshTag* distT = initTag(m, "parmaDistance", initVal);

    APF_ITERATE(Level, verts, itr)
      parma::dijkstra(m,connT,distT,*itr);

    apf::removeTagFromDimension(m,connT,m->getDimension());
    m->destroyTag(connT);

    parmaCommons::printElapsedTime("computeDistance", PCU_Time()-t0);
    return distT;
  }

} //end namespace

namespace parma {
  apf::MeshTag* measureGraphDist(apf::Mesh* m) {
    Level* centralVerts = getCentralVerts(m);
    apf::MeshTag* t = computeDistance(m, *centralVerts);
    delete centralVerts;
    return t;
  }

} //end namespace parma
