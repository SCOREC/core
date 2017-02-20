#include <apf.h>
#include <PCU.h>
#include "parma_graphDist.h"
#include "parma_dijkstra.h"
#include "parma_dcpart.h"
#include "parma_meshaux.h"
#include "parma_convert.h"
#include "parma_commons.h"
#include <list>
#include <set>
#include <limits.h>
#include <stdlib.h>

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
    const unsigned check = m->getTagChecksum(dt,apf::Mesh::VERTEX);
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
    assert(check == m->getTagChecksum(dt,apf::Mesh::VERTEX));
    return rmax;
  }

  void warnAboutDistanceProblem(const unsigned csA, const unsigned csB,
      const unsigned dtchanges) {
    if( !dtchanges && csA != csB ) {
      parmaCommons::error(
          "rank %d there were no distance changes but the checksum "
          "on the distance array does not match\n",
          PCU_Comm_Self());
    }
    if( dtchanges && csA == csB ) {
      parmaCommons::error(
          "rank %d there were distance changes but the checksum "
          "on the distance array did not change\n",
          PCU_Comm_Self());
    }
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
    const unsigned csStart = m->getTagChecksum(dt,apf::Mesh::VERTEX);
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
    unsigned dtchanges = 0;
    for(unsigned i=c.size()-1; i>0; i--) {
      apf::MeshEntity* v;
      c.beginBdry(i);
      while( (v = c.iterateBdry()) ) {
        int d; m->getIntTag(v,dt,&d);
        int rsi = TO_INT(rsum[i]);
        if(d < rsi) { //not visited
          d+=rsi;
          m->setIntTag(v,dt,&d);
          dtchanges++;
        }
      }
      c.endBdry();
    }
    const unsigned csMid = m->getTagChecksum(dt,apf::Mesh::VERTEX);
    warnAboutDistanceProblem(csStart,csMid,dtchanges);
    assert((!dtchanges && csStart == csMid) || (dtchanges && csStart != csMid));

    // Offset the interior vertices
    dtchanges = 0;
    apf::MeshEntity* v;
    apf::MeshIterator* it = m->begin(0);
    while( (v = m->iterate(it)) ) {
      //skip if not part of a component
      if( !c.has(v) ) continue;
      unsigned id = c.getId(v);
      //also skip if on a component boundary
      if( c.bdryHas(id,v) ) continue;
      //also skip if the offset is zero
      if( !rsum[id] ) continue;
      int d; m->getIntTag(v,dt,&d);
      int dist = d + TO_INT(rsum[id]);
      m->setIntTag(v,dt,&dist);
      if(dist != d) dtchanges++;
    }
    m->end(it);
    const unsigned csEnd = m->getTagChecksum(dt,apf::Mesh::VERTEX);
    warnAboutDistanceProblem(csMid,csEnd,dtchanges);
    assert((!dtchanges && csMid == csEnd) || (dtchanges && csMid != csEnd));
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
      getElmAdjVtx(m,u,verts);
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

namespace parma_ordering {
  typedef std::list<apf::MeshEntity*> queue;

  inline apf::MeshEntity* pop(queue& q) {
    apf::MeshEntity* e = q.front();
    assert(e);
    q.pop_front();
    return e;
  }

  int bfs(apf::Mesh* m, parma::DijkstraContains* c,
       apf::MeshEntity* src, apf::MeshTag* order, int num) {
    if( !src )
      return num;
    queue q;
    q.push_back(src);
    while( !q.empty() ) {
      apf::MeshEntity* v = pop(q);
      if( m->hasTag(v,order) ) continue;
      m->setIntTag(v,order,&num); num++;
      apf::Adjacent adjVtx;
      getEdgeAdjVtx(m,v,adjVtx);
      APF_ITERATE(apf::Adjacent, adjVtx, eItr) {
        apf::MeshEntity* u = *eItr;
        if( c->has(u) && ! m->hasTag(u,order) )
          q.push_back(u);
      }
    }
    return num;
  }

  apf::MeshEntity* getMaxDistSeed(apf::Mesh* m, CompContains* c,
      apf::MeshTag* dt, apf::MeshTag* order) {
    int rmax = -1;
    apf::MeshEntity* emax = NULL;
    int cnt=0;
    apf::MeshIterator* it = m->begin(0);
    apf::MeshEntity* e;
    while( (e = m->iterate(it)) ) {
      if( ! c->has(e) ) continue;
      cnt++;
      int d; m->getIntTag(e,dt,&d);
      PCU_Debug_Print("cnt %d d %d hasTag %d\n", cnt, d, m->hasTag(e,order));
      if( !m->hasTag(e,order) && d > rmax ) {
        rmax = d;
        emax = e;
      }
    }
    m->end(it);
    return emax;
  }

  apf::MeshTag* reorder(apf::Mesh* m, parma::dcComponents& c, apf::MeshTag* dist) {
    const unsigned check = c.getIdChecksum();
    apf::MeshTag* order = m->createIntTag("parma_ordering",1);
    int start = 0;
    for(int i=TO_INT(c.size())-1; i>=0; i--) {
      CompContains* contains = new CompContains(c,i);
      apf::MeshEntity* src = getMaxDistSeed(m,contains,dist,order);
      PCU_Debug_Print("comp %d starting vertex found? %d\n", i, (src != NULL));
      start = bfs(m, contains, src, order, start);
      assert(check == c.getIdChecksum());
      delete contains;
      if(start == TO_INT(m->count(0))) {
        if( i != 0 ) //if not the last component to order
          parmaCommons::status("%d all vertices visited comp %u of %u\n",
              PCU_Comm_Self(), i, c.size());
        break;
      }
    }
    assert(start == TO_INT(m->count(0)));

    int* sorted = new int[m->count(0)];
    for(unsigned i=0; i<m->count(0); i++)
      sorted[i] = 0;
    apf::MeshIterator* it = m->begin(0);
    apf::MeshEntity* e;
    while( (e = m->iterate(it)) ) {
      assert( m->hasTag(e,order) );
      int id; m->getIntTag(e,order,&id);
      assert(id < TO_INT(m->count(0)));
      sorted[id] = 1;
    }
    m->end(it);
    for(unsigned i=0; i<m->count(0); i++)
      assert(sorted[i]);
    delete [] sorted;
    assert(check == c.getIdChecksum());
    return order;
  }

  // linear arrangement - estimate effectiveness of reordering
  void la(apf::Mesh* m, apf::MeshTag* order=NULL) {
    int setOrder = (order == NULL);
    if( setOrder ) {
      order = m->createIntTag("parma_default_ordering",1);
      apf::MeshIterator* it = m->begin(0);
      apf::MeshEntity* e;
      int i = 0;
      while( (e = m->iterate(it)) ) {
        m->setIntTag(e,order,&i);
        i++;
      }
      m->end(it);
    }
    const unsigned check = m->getTagChecksum(order,apf::Mesh::VERTEX);
    int la = 0;
    apf::Downward verts;
    apf::MeshIterator* it = m->begin(1);
    apf::MeshEntity* e;
    while( (e = m->iterate(it)) ) {
      m->getDownward(e,0,verts);
      int vid; m->getIntTag(verts[0],order,&vid);
      int uid; m->getIntTag(verts[1],order,&uid);
      la += abs(vid-uid);
    }
    m->end(it);
    PCU_Debug_Print("la %d\n", la);
    long tot=PCU_Add_Long(TO_LONG(la));
    int max=PCU_Max_Int(la);
    int min=PCU_Min_Int(la);
    double avg = TO_DOUBLE(tot)/PCU_Comm_Peers();
    if( !PCU_Comm_Self() )
      parmaCommons::status("la min %d max %d avg %.3f\n", min, max, avg);
    assert(check == m->getTagChecksum(order,apf::Mesh::VERTEX));
    if( setOrder )
      m->destroyTag(order);
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
          parmaCommons::error("rank %d comp %u iso %u ... "
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

apf::MeshTag* Parma_BfsReorder(apf::Mesh* m, int) {
  double t0 = PCU_Time();
  assert( !hasDistance(m) );
  parma::dcComponents c = parma::dcComponents(m);
  const unsigned checkIds = c.getIdChecksum();
  apf::MeshTag* dist = computeDistance(m,c);
  const unsigned check = m->getTagChecksum(dist,apf::Mesh::VERTEX);
  if( PCU_Comm_Peers() > 1 && !c.numIso() )
    if( !hasDistance(m,dist) ) {
      parmaCommons::error("rank %d comp %u iso %u ... "
          "some vertices don't have distance computed\n",
          PCU_Comm_Self(), c.size(), c.numIso());
      assert(false);
    }
  parma_ordering::la(m);
  apf::MeshTag* order = parma_ordering::reorder(m,c,dist);
  parma_ordering::la(m,order);
  assert(checkIds == c.getIdChecksum());
  assert(check == m->getTagChecksum(dist,apf::Mesh::VERTEX));
  m->destroyTag(dist);
  parmaCommons::printElapsedTime(__func__,PCU_Time()-t0);
  return order;
}
