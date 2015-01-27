#include <PCU.h>
#include "parma_vtxSelector.h"
#include "parma_targets.h"
#include "parma_weights.h"
#include "parma_graphDist.h"
#include "parma_bdryVtx.h"
#include "parma_commons.h"
#include <apf.h>
#include <limits.h>

typedef unsigned int uint;

namespace {
  struct UintArr {
    uint s; //size of d array
    uint l; //used entries in d
    uint d[1];
  };

  UintArr* makeUintArr(uint n) {
    UintArr* a = (UintArr*) malloc(sizeof(UintArr) + sizeof(uint)*(n-1));
    a->s = n;
    a->l = 0;
    return a;
  }
  void destroy(UintArr* u) {
    free(u);
  }

  UintArr* getCavityPeers(apf::Mesh* m, apf::MeshEntity* v) {
    typedef std::map<uint,uint> MUiUi;
    MUiUi pc;
    apf::Adjacent sideSides;
    m->getAdjacent(v, m->getDimension()-2, sideSides);
    APF_ITERATE(apf::Adjacent, sideSides, ss) {
      apf::Copies rmts;
      m->getRemotes(*ss,rmts);
      APF_ITERATE(apf::Copies, rmts, r)
         pc[static_cast<uint>(r->first)]++;
    }
    uint max = 0;
    APF_ITERATE(MUiUi, pc, p)
      if( p->second > max )
         max = p->second;
    UintArr* peers = makeUintArr(pc.size());
    APF_ITERATE(MUiUi, pc, p)
      if( p->second == max )
        peers->d[peers->l++] = p->first;
    return peers;
  }

  void getCavity(apf::Mesh* m, apf::MeshEntity* v, apf::Migration* plan, 
      apf::Up& cavity) {
    cavity.n = 0;
    apf::Adjacent elms;
    m->getAdjacent(v, m->getDimension(), elms);
    APF_ITERATE(apf::Adjacent, elms, adjItr)
      if( !plan->has(*adjItr) )
        cavity.e[(cavity.n)++] = *adjItr; 
  }
}

namespace parma {
  VtxSelector::VtxSelector(apf::Mesh* m, apf::MeshTag* w)
    : Selector(m, w)
  {
    dist = measureGraphDist(m);
  }

  VtxSelector::~VtxSelector() {
    apf::removeTagFromDimension(mesh,dist,0);
    mesh->destroyTag(dist);
  }

  apf::Migration* VtxSelector::run(Targets* tgts) {
    apf::Migration* plan = new apf::Migration(mesh);
    double planW = 0;
    for(int max=2; max <= 12; max+=2)
      planW += select(tgts, plan, planW, max);
    return plan;
  }

  double VtxSelector::getWeight(apf::MeshEntity* e) {
    return getEntWeight(mesh,e,wtag);
  }

  double VtxSelector::add(apf::MeshEntity* v, apf::Up& cavity, const int destPid,
      apf::Migration* plan) {
    for(int i=0; i < cavity.n; i++)
      plan->send(cavity.e[i], destPid);
    return getWeight(v);
  }

  double VtxSelector::select(Targets* tgts, apf::Migration* plan, double planW,
      int maxSize) {
    double t0 = PCU_Time();
    BdryVtxItr* bdryVerts = makeBdryVtxDistItr(mesh, dist);
    apf::Up cavity;
    apf::MeshEntity* e;
    while( (e = bdryVerts->next()) ) {
      if( planW > tgts->total() ) break;
      getCavity(mesh, e, plan, cavity);
      UintArr* peers = getCavityPeers(mesh,e);
      int d; mesh->getIntTag(e,dist,&d);
      for( uint i=0; i<peers->l; i++ ) {
        uint destPid = peers->d[i];
        if( (tgts->has(destPid) &&
              sending[destPid] < tgts->get(destPid) &&
              cavity.n <= maxSize ) ||
            INT_MAX == d ) {
          double ew = add(e, cavity, destPid, plan);
          sending[destPid] += ew;
          planW += ew;
          break;
        }
      }
      destroy(peers);
    }
    parmaCommons::printElapsedTime("select", PCU_Time()-t0);
    delete bdryVerts;
    return planW;
  }

  Selector* makeVtxSelector(apf::Mesh* m, apf::MeshTag* w) {
    return new VtxSelector(m, w);
  }
} //end namespace parma
