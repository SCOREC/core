#include <PCU.h>
#include <apf.h>
#include <apfMesh.h>
#include "parma_selector.h"
#include "parma_targets.h"
#include "parma_weights.h"
#include "parma_convert.h"
#include <stdlib.h>

typedef unsigned uint;
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

  void getCavityPeers(apf::Mesh* m, apf::MeshEntity* v, UintArr* peers) {
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
    assert( peers->s >= pc.size() );
    peers->l=0;
    APF_ITERATE(muu, pc, p)
      if( p->second == max )
        peers->d[peers->l++] = p->first;
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

  bool sharedWithPeer(apf::Mesh* m, apf::MeshEntity* e, int peer) {
    apf::Parts res;
    m->getResidence(e,res);
    return res.count(peer);
  }

  bool sharedWithTarget(apf::Mesh* m, apf::MeshEntity* e, parma::Targets* t) {
    const parma::Targets::Item* tgt;
    t->begin();
    int shared = 0;
    while( (tgt = t->iterate()) )
      shared += sharedWithPeer(m,e,tgt->first);
    t->end();
    return shared;
  }

  void add(apf::Up& cavity, const int destPid,
      apf::Migration* plan) {
    for(int i=0; i < cavity.n; i++)
      plan->send(cavity.e[i], destPid);
  }

  class ShapeSelector : public parma::Selector {
    public:
      ShapeSelector(apf::Mesh* m, apf::MeshTag* w)
        : Selector(m, w) { }

      /**
       * @brief eliminate sides with neighbors in Targets
       */
      apf::Migration* run(parma::Targets* tgts) {
        apf::Migration* plan = new apf::Migration(mesh);
        apf::Up cavity;
        UintArr* peers = makeUintArr(1024);
        apf::MeshEntity* e;
        apf::MeshIterator* it = mesh->begin(0);
        while( (e = mesh->iterate(it)) ) {
          if( !mesh->isShared(e) ) continue;
          if( !sharedWithTarget(mesh,e,tgts) ) continue;
          getCavity(mesh, e, plan, cavity);
          getCavityPeers(mesh,e,peers);
          for( uint i=0; i<peers->l; i++ ) {
            int destPid = TO_INT(peers->d[i]);
            if( !tgts->has(destPid) ) {
              add(cavity, destPid, plan);
              break;
            }
          }
        }
        destroy(peers);
        mesh->end(it);
        return plan;
      }
  };
}

namespace parma {
  Selector* makeShapeSelector(apf::Mesh* m, apf::MeshTag* wtag) {
    return new ShapeSelector(m, wtag);
  }
}

