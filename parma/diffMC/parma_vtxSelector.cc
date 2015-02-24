#include <PCU.h>
#include "parma.h"
#include "parma_vtxSelector.h"
#include "parma_targets.h"
#include "parma_weights.h"
#include "parma_graphDist.h"
#include "parma_bdryVtx.h"
#include "parma_commons.h"
#include <apf.h>

typedef unsigned int uint;
#define TO_UINT(a) static_cast<unsigned>(a)
#define TO_INT(a) static_cast<int>(a)
typedef std::map<uint,uint> muu;

namespace {
  struct UintArr {
    uint s; //size of d array
    uint l; //used entries in d
    uint d[1];
  };

  UintArr* makeUintArr(uint n) {
    UintArr* a = 
      static_cast<UintArr*>(malloc(sizeof(UintArr) + sizeof(uint)*(n-1)));
    a->s = n;
    a->l = 0;
    return a;
  }
  void destroy(UintArr* u) {
    free(u);
  }

  UintArr* getCavityPeers(apf::Mesh* m, apf::MeshEntity* v) {
    muu pc;
    apf::Adjacent sideSides;
    m->getAdjacent(v, m->getDimension()-2, sideSides);
    APF_ITERATE(apf::Adjacent, sideSides, ss) {
      apf::Copies rmts;
      m->getRemotes(*ss,rmts);
      APF_ITERATE(apf::Copies, rmts, r)
         pc[TO_UINT(r->first)]++;
    }
    uint max = 0;
    APF_ITERATE(muu, pc, p)
      if( p->second > max )
         max = p->second;
    UintArr* peers = makeUintArr(TO_UINT(pc.size()));
    APF_ITERATE(muu, pc, p)
      if( p->second == max )
        peers->d[peers->l++] = p->first;
    return peers;
  }

  typedef std::map<apf::MeshEntity*,unsigned> meu;
  // construct a map of <faces, occurances in cavity>
  meu* getCavityFaces(apf::Mesh* m, apf::Up& cavity) {
    const int dim = m->getDimension();
    meu* cf = new meu;
    for(int i=0; i<cavity.n; i++) {
      apf::Downward faces;
      int nf = m->getDownward(cavity.e[i], dim-1, faces);
      for(int j=0; j<nf; j++)
        (*cf)[faces[j]]++;
    }
    return cf;
  }

  // return true if face-connected to the part
  bool disconnected(apf::Mesh* m, apf::Migration* plan, apf::Up& cavity) {
    meu* faces = getCavityFaces(m, cavity);
    bool dc = true;
    APF_ITERATE(meu, *faces, f) {
      if( 2 == f->second ) continue; //cavity interior
      apf::Up elms;
      m->getUp(f->first, elms);
      if( elms.n == 2 &&
          !plan->has(elms.e[0]) &&
          !plan->has(elms.e[1]) ) {
        dc = false;
        break;
      }
    }
    delete faces;
    return dc;
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

  VtxSelector::~VtxSelector() { }

  apf::Migration* VtxSelector::run(Targets* tgts) {
    apf::Migration* plan = new apf::Migration(mesh);
    double t0 = PCU_Time();
    double planW = 0;
    for(int max=2; max <= 12; max+=2)
      planW += select(tgts, plan, planW, max);
    parmaCommons::printElapsedTime("select", PCU_Time()-t0);
    clearDistTag(mesh, dist, plan);
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
    BdryVtxItr* bdryVerts = makeBdryVtxDistItr(mesh, dist);
    apf::Up cavity;
    apf::MeshEntity* e;
    unsigned dcCnt = 0;
    while( (e = bdryVerts->next()) ) {
      if( planW > tgts->total() ) break;
      getCavity(mesh, e, plan, cavity);
      UintArr* peers = getCavityPeers(mesh,e);
      bool sent = false;
      for( uint i=0; i<peers->l; i++ ) {
        int destPid = TO_INT(peers->d[i]);
        if( tgts->has(destPid) &&
            sending[destPid] < tgts->get(destPid) &&
            cavity.n <= maxSize ) {
          double ew = add(e, cavity, destPid, plan);
          sending[destPid] += ew;
          planW += ew;
          sent = true;
          break;
        }
      }
      if( !sent && disconnected(mesh, plan, cavity) ) {
        assert(peers->l);
        int destPid = TO_INT(peers->d[0]);
        dcCnt++;
        double ew = add(e, cavity, destPid, plan);
        sending[destPid] += ew;
        planW += ew;
      }
      destroy(peers);
    }
    PCU_Debug_Print("sent %u disconnected cavities\n", dcCnt);
    delete bdryVerts;
    return planW;
  }

  Selector* makeVtxSelector(apf::Mesh* m, apf::MeshTag* w) {
    return new VtxSelector(m, w);
  }
} //end namespace parma
