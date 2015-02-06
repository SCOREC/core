#include <apf.h>
#include <PCU.h>
#include <list>
#include <set>
#include <limits.h>
#include "parma_graphDist.h"
#include "parma_dijkstra.h"
#include "parma_dcpart.h" 

#define TO_UINT(a) static_cast<unsigned>(a)
#define TO_INT(a) static_cast<int>(a)

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

  unsigned* getMaxDist(apf::Mesh* m, parma::dcComponents& c, apf::MeshTag* dt) {
    unsigned* rmax = new unsigned[c.size()];
    for(unsigned i=0; i<c.size(); i++) {
      rmax[i] = 0;
      apf::MeshEntity* v;
      c.beginBdry(i);
      while( (v = c.iterateBdry()) ) {
        int d; m->getIntTag(v,dt,&d);
        unsigned du = TO_UINT(d);
        if( du > rmax[i] )
          rmax[i] = du;
      }
      c.endBdry();
    }
    return rmax;
  }

  void offset(apf::Mesh* m, parma::dcComponents& c, 
      apf::MeshTag* dt, unsigned* rmax) {
    if (!c.size())
      return;
    unsigned* rsum = new unsigned[c.size()];
    rsum[0] = 0;
    for(unsigned i=1; i<c.size(); i++)
      rsum[i] = rsum[i-1] + rmax[i-1] + 1;

    apf::MeshEntity* v;
    apf::MeshIterator* it = m->begin(0);
    while( (v = m->iterate(it)) )
      if( c.has(v) ) {
        unsigned id = c.getId(v);
        int d; m->getIntTag(v,dt,&d);
        d+=TO_INT(rsum[id]);
        m->setIntTag(v,dt,&d);
      }
    m->end(it);
    delete [] rsum;
  }


  class CompContains : public parma::DijkstraContains {
    public:
      CompContains(parma::dcComponents& comps, unsigned compid) 
        : c(comps), id(compid) {}
      ~CompContains() {}
      // Without the onBdry check some boundary vertices would be excluded from
      // distancing.  When there are multiple boundary and interior vertices
      // this is generally not a problem, but for cases where the component is a
      // single element the core vertex could be a boundary vertex that is not
      // assiged to the component as queried by getId(e).
      bool has(apf::MeshEntity* e) {
        bool onBdry = c.bdryHas(id,e);
        bool inComp = (c.has(e) && c.getId(e) == id);
        return ( inComp || onBdry );
      }
    private:
      parma::dcComponents& c;
      unsigned id;
  };

  apf::MeshTag* computeDistance(apf::Mesh* m, parma::dcComponents& c) {
    int initVal = INT_MAX;
    apf::MeshTag* distT = initTag(m, "parmaDistance", initVal);
    for(unsigned i=0; i<c.size(); i++) {
      CompContains* contains = new CompContains(c,i);
      apf::MeshEntity* src = c.getCore(i);
      parma::dijkstra(m, contains, src, distT);
      delete contains;
    }
    return distT;
  }

  bool hasDistance(apf::Mesh* m, apf::MeshTag* dist) {
    apf::MeshEntity* e;
    apf::MeshIterator* it = m->begin(0);
    while( (e = m->iterate(it)) ) {
      int d;
      m->getIntTag(e,dist,&d);
      if( d == INT_MAX )
        return false;
    }
    m->end(it);
    return true;
  }
} //end namespace

namespace parma {
  apf::MeshTag* measureGraphDist(apf::Mesh* m) {
    dcComponents c = dcComponents(m);
    apf::MeshTag* t = computeDistance(m,c);
    if( PCU_Comm_Peers() > 1 && !c.numIso() )
      if( !hasDistance(m,t) ) {
        fprintf(stderr, "CAKE rank %d comp %u iso %u ... "
            "some vertices don't have distance computed\n", 
            PCU_Comm_Self(), c.size(), c.numIso());
        assert(false);
      }
    unsigned* rmax = getMaxDist(m,c,t);
    offset(m,c,t,rmax);
    delete [] rmax;
    return t;
  }
}
