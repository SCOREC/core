#include "parma_selector.h"
#include "parma_targets.h"
#include "parma_weights.h"
#include "parma_graphDist.h"
#include "parma_bdryVtx.h"
#include "parma_commons.h"
#include <apf.h>
#include <PCU.h>
#include <set>
#include <stdio.h>
#include <limits.h>
#include <vector>

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

  typedef std::set<apf::MeshEntity*> SetEnt;
  typedef std::map<int,double> Mid;

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
  class VtxSelector : public Selector {
    public:
      VtxSelector(apf::Mesh* m, apf::MeshTag* w)
        : Selector(m, w)
      {
        dist = measureGraphDist(m);
      }
      ~VtxSelector() {
        apf::removeTagFromDimension(mesh,dist,0);
        mesh->destroyTag(dist);
      }
      apf::Migration* run(Targets* tgts) {
        apf::Migration* plan = new apf::Migration(mesh);
        double planW = 0;
        for(int max=2; max <= 12; max+=2)
          planW += select(tgts, plan, planW, max);
        return plan;
      }
    protected:
      virtual double getWeight(apf::MeshEntity* e) {
        return getEntWeight(mesh,e,wtag);
      }
      virtual double add(apf::MeshEntity* v, apf::Up& cavity, const int destPid,
          apf::Migration* plan) {
        for(int i=0; i < cavity.n; i++)
           plan->send(cavity.e[i], destPid);
        return getWeight(v);
      }
      double select(Targets* tgts, apf::Migration* plan, double planW,
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
    private:
      apf::MeshTag* dist;
      VtxSelector();
      Mid sending;
  };
  Selector* makeVtxSelector(apf::Mesh* m, apf::MeshTag* w) {
    return new VtxSelector(m, w);
  }

  class EdgeSelector : public VtxSelector {
    public:
      EdgeSelector(apf::Mesh* m, apf::MeshTag* w)
        : VtxSelector(m, w) {}
    protected:
      // if all the elements bounded by an edge are in the plan then it is being
      // sent and its weight is counted
      bool sent(apf::Migration* plan, apf::MeshEntity* e) {
        apf::Adjacent elms;
        mesh->getAdjacent(e, mesh->getDimension(), elms);
        APF_ITERATE(apf::Adjacent, elms, elm)
          if( ! plan->has(*elm) ) 
            return false;
        return true;
      }

      double cavityWeight(apf::Migration* plan, SetEnt& s) {
        double w = 0;
        APF_ITERATE(SetEnt, s, sItr)
          if( sent(plan,*sItr) )
            w += getWeight(*sItr);
        return w;
      }

      void addCavityEdge(apf::MeshEntity* e, SetEnt& s) {
        apf::Adjacent adjEdge;
        mesh->getAdjacent(e, 1, adjEdge);
        APF_ITERATE(apf::Adjacent, adjEdge, adjItr)
          s.insert(*adjItr);
      }

      virtual double add(apf::MeshEntity*, apf::Up& cavity, const int destPid,
          apf::Migration* plan) {
        SetEnt s;
        for(int i=0; i < cavity.n; i++) {
           addCavityEdge(cavity.e[i], s); 
           plan->send(cavity.e[i], destPid);
        }
        return cavityWeight(plan, s);
      }
  };
  Selector* makeEdgeSelector(apf::Mesh* m, apf::MeshTag* w) {
    return new EdgeSelector(m, w);
  }

  class ElmSelector : public VtxSelector {
    public:
      ElmSelector(apf::Mesh* m, apf::MeshTag* w)
        : VtxSelector(m, w) {}
    protected:
      virtual double add(apf::MeshEntity*, apf::Up& cavity, const int destPid,
          apf::Migration* plan) {
        double w = 0;
        for(int i=0; i < cavity.n; i++) {
           plan->send(cavity.e[i], destPid);
           w += getWeight(cavity.e[i]);
        }
        return w;
      }
  };
  Selector* makeElmSelector(apf::Mesh* m, apf::MeshTag* w) {
    return new ElmSelector(m, w);
  }

  class ElmLtVtxSelector : public VtxSelector {
    public:
      ElmLtVtxSelector(apf::Mesh* m, apf::MeshTag* w, double maxV)
        : VtxSelector(m, w), maxVtx(maxV) {}
      apf::Migration* run(Targets* tgts) {
        apf::Migration* plan = new apf::Migration(mesh);
        double planW = 0;
        for(int max=2; max <= 12; max+=2)
          planW += select(tgts, plan, planW, max);
        cancel(&plan, trim(tgts));
        return plan;
      }
    protected:
      void addCavityVtx(apf::MeshEntity* e, SetEnt& s) {
        apf::Adjacent adjVtx;
        mesh->getAdjacent(e, 0, adjVtx);
        APF_ITERATE(apf::Adjacent, adjVtx, adjItr)
          s.insert(*adjItr);
      }

      double cavityWeight(SetEnt& s) {
        double w = 0;
        APF_ITERATE(SetEnt, s, sItr)
          w += getWeight(*sItr);
        return w;
      }

      double cavityWeight(SetEnt& s, int dest) {
        double w = 0;
        APF_ITERATE(SetEnt, s, sItr) {
          apf::Copies rmts;
          mesh->getRemotes(*sItr,rmts);
          bool shared = false;
          APF_ITERATE(apf::Copies, rmts, r)
            if( r->first == dest ) {
              shared = true;
              break;
            }
          if( !shared )
            w += getWeight(*sItr);
        }
        return w;
      }

      void cancel(apf::Migration** plan, Mid* order) {

        apf::Migration* planA = *plan;
        typedef std::pair<apf::MeshEntity*, int> PairEntInt;
        std::vector<PairEntInt > keep;
        const size_t planSz = static_cast<size_t>(planA->count());
        keep.reserve(planSz);

        std::map<int,SetEnt > peerSelections;
        for(int i=0; i < planA->count(); i++) {
           apf::MeshEntity* e = planA->get(i);
           int dest = planA->sending(e);
           SetEnt vset = peerSelections[dest];
           addCavityVtx(e, vset);
           if( cavityWeight(vset) <= (*order)[dest] ) {
             keep.push_back(PairEntInt(e,dest));
             peerSelections[dest] = vset;
           }
        }
        delete order;
        delete planA;
        *plan = new apf::Migration(mesh);
        for(size_t i=0; i < keep.size(); i++)
          (*plan)->send(keep[i].first, keep[i].second);
      }

      typedef std::pair<int,double> Migr;
      struct CompareMigr {
        bool operator()(const Migr& a, const Migr& b) const {
          if( a.second < b.second )
            return true;
          else
            return false;
        }
      };

      Mid* trim(Targets*) {
        typedef std::set<Migr,CompareMigr> MigrComm;

        PCU_Comm_Begin();
        APF_ITERATE(Mid, sendingVtx, s) {
          PCU_COMM_PACK(s->first, s->second);
        }
        PCU_Comm_Send();

        MigrComm incoming;
        double w;
        while (PCU_Comm_Listen()) {
          PCU_COMM_UNPACK(w);
          incoming.insert(Migr(PCU_Comm_Sender(),w));
        }

        double selfW = parma::getWeight(mesh,wtag,0);
        Mid accept;
        double totW = selfW;
        APF_ITERATE(MigrComm, incoming, in) {
          double avail = maxVtx - totW;
          if( avail > 0 )
            if( in->second <= avail )
              accept[in->first] = in->second;
            else
              accept[in->first] = avail;
          else
            accept[in->first] = 0;
          totW += accept[in->first];
        }

        PCU_Comm_Begin();
        APF_ITERATE(Mid, accept, a)
          PCU_COMM_PACK(a->first, a->second);
        PCU_Comm_Send();
        Mid* order = new Mid;
        double outw;
        while (PCU_Comm_Listen()) {
          PCU_COMM_UNPACK(outw);
          (*order)[PCU_Comm_Sender()] = outw;
        }
        return order;
      }

      virtual double add(apf::MeshEntity*, apf::Up& cavity,
          const int destPid, apf::Migration* plan) {
        SetEnt cav;
        double w = 0;
        for(int i=0; i < cavity.n; i++) {
          addCavityVtx(cavity.e[i], cav);
          plan->send(cavity.e[i], destPid);
          w += getWeight(cavity.e[i]);
        }
        sendingVtx[destPid] += cavityWeight(cav);
        return w;
      }
    private:
      Mid sendingVtx;
      double maxVtx;
  };
  Selector* makeElmLtVtxSelector(apf::Mesh* m, apf::MeshTag* w, double maxVtx) {
    return new ElmLtVtxSelector(m, w, maxVtx);
  }
} //end namespace parma
