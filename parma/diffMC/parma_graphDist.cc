#include <apf.h>
#include <PCU.h>
#include <list>
#include <set>
#include <limits.h>
#include "parma_graphDist.h"
#include "parma_dijkstra.h"
#include "parma_dcpart.h"
#include "parma_meshaux.h"

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
    //If maxDistanceIncrease number of diffusion steps were ran there could be
    //at most an increase in the distance of maxDistanceIncrease.  This can be
    //seen with the worst case of a mesh of a one element thick rectangle that
    //has N elements along it's length, where N = 2*maxDistanceIncrease.
    //Suppose that one element is assigned to part zero and N-1 elements
    //assigned to part one.  Diffusion will only be able to migrate one element
    //per step from part one to zero.  The max distance of part zero increases
    //by one each step.  Thus in maxDistanceIncrease steps the distance can at
    //most increase by maxDistanceIncrease for a given component.
    const int maxDistanceIncrease = 1000;
    if (!c.size())
      return;
    unsigned* rsum = new unsigned[c.size()];
    rsum[0] = 0;
    for(unsigned i=1; i<c.size(); i++)
      rsum[i] = rsum[i-1] + rmax[i-1] + 1 + maxDistanceIncrease;

    for(unsigned i=1; i<c.size(); i++)
      PCU_Debug_Print("offset %u is %u\n", i, rsum[i]);

    // Go backwards so that the largest bdry vtx changes are made first
    //  and won't be augmented in subsequent bdry traversals.
    for(unsigned i=c.size()-1; i>0; i--) {
      apf::MeshEntity* v;
      c.beginBdry(i);
      while( (v = c.iterateBdry()) ) {
        int d; m->getIntTag(v,dt,&d);
        int rsi = TO_INT(rsum[i]);
        if(d < rsi) { //not visited
          d+=rsi;
          m->setIntTag(v,dt,&d);
        }
      }
      c.endBdry();
    }

    // Offset the interior vertices
    apf::MeshEntity* v;
    apf::MeshIterator* it = m->begin(0);
    while( (v = m->iterate(it)) ) {
      //skip if
      if( !c.has(v) ) continue; // not part of a component
      unsigned id = c.getId(v);
      //also skip if on a component boundary
      if( c.bdryHas(id,v) ) continue;
      int d; m->getIntTag(v,dt,&d);
      d += TO_INT(rsum[id]);
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
      bool bdryHas(apf::MeshEntity* e) {
        return c.bdryHas(id,e);
      }
    private:
      parma::dcComponents& c;
      unsigned id;
  };

  const char* distanceTagName() {
    return "parmaDistance";
  }

  apf::MeshTag* createDistTag(apf::Mesh* m) {
    int initVal = INT_MAX;
    return initTag(m, distanceTagName(), initVal);
  }

  apf::MeshTag* computeDistance(apf::Mesh* m, parma::dcComponents& c) {
    apf::MeshTag* distT = createDistTag(m);
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

  bool hasDistance(apf::Mesh* m) {
    apf::MeshTag* dist = parma::getDistTag(m);
    if( dist )
      return true;
    else
      return false;
  }

  inline bool onMdlBdry(apf::Mesh* m, apf::MeshEntity* v) {
    const int dim = m->getDimension();
    apf::ModelEntity* g = m->toModel(v);
    const int gdim = m->getModelType(g);
    return gdim < dim;
  }

  /**
   * @brief construct list of part boundary and geometric boundary
   *   vertices to have their distance updated
   * @remark There are two side effects. Each vertex in the list:
   *   (1) has the distance set to INT_MAX.  Any vertices that are not visited
   *   during the walks will have a INT_MAX distance which will bias diffusion
   *   to migrate the bounded cavities before other cavities.
   */
  void getBdryVtx(apf::Mesh* m, apf::MeshTag* dist, 
      parma::DistanceQueue<parma::Less>& pq) {
    int dmax = INT_MAX;
    apf::MeshEntity* u;
    apf::MeshIterator* it = m->begin(0);
    while( (u = m->iterate(it)) ) {
      if( !m->isShared(u) && !onMdlBdry(m,u) ) continue;
      m->setIntTag(u,dist,&dmax); // (1)
      apf::Adjacent verts;
      getEdgeAdjVtx(m,u,verts);
      APF_ITERATE(apf::Adjacent, verts, v) {
        if( !m->isShared(*v) && !onMdlBdry(m,*v) ) {
          int vd; m->getIntTag(*v,dist,&vd);
          if( vd == INT_MAX ) continue;
          pq.push(*v,vd);
        }
      }
    }
    m->end(it);
  }

  class CompUpdateContains : public parma::DijkstraContains {
    public:
      CompUpdateContains() {}
      ~CompUpdateContains() {}
      bool has(apf::MeshEntity*) {
        return true;
      }
      //disable the non-manifold boundary detection mechanism
      bool bdryHas(apf::MeshEntity*) {
        return false;
      }
  };

  apf::MeshTag* updateDistance(apf::Mesh* m) {
    PCU_Debug_Print("updateDistance\n");
    apf::MeshTag* dist = parma::getDistTag(m);
    parma::DistanceQueue<parma::Less> pq(m);
    getBdryVtx(m,dist,pq);
    CompUpdateContains c;
    parma::dijkstra(m,&c,pq,dist);
    return dist;
  }
} //end namespace

namespace parma {
  apf::MeshTag* getDistTag(apf::Mesh* m) {
    return m->findTag(distanceTagName());
  }

  apf::MeshTag* measureGraphDist(apf::Mesh* m) {
    apf::MeshTag* t = NULL;
    if( hasDistance(m) ) {
      t = updateDistance(m);
    } else {
      PCU_Debug_Print("computeDistance\n");
      dcComponents c = dcComponents(m);
      t = computeDistance(m,c);
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
    }
    return t;
  }
}
