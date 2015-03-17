#include <set>
#include <apf.h>
#include <PCU.h>
#include "parma_vtxSelector.h"
#include "parma_targets.h"
#include "parma_weights.h"

#define TO_SIZET(a) static_cast<size_t>(a)

namespace {
  typedef std::set<apf::MeshEntity*> SetEnt;

  class ElmLtVtxSelector : public parma::VtxSelector {
    private:
      double maxVtx;

    public:
      ElmLtVtxSelector(apf::Mesh* m, apf::MeshTag* w, double maxV)
        : VtxSelector(m, w), maxVtx(maxV) { }

      apf::Migration* run(parma::Targets* tgts) {
        apf::Migration* plan = new apf::Migration(mesh);
        double planW = 0;
        for(int max=2; max <= 12; max+=2)
          planW += select(tgts, plan, planW, max);
        parma::Mid* capacity = trim(tgts,plan);
        cancel(&plan, capacity);
        return plan;
      }

    protected:
      void insertInteriorVerts(apf::MeshEntity* e, int dest, SetEnt& s) {
        apf::Adjacent adjVtx;
        mesh->getAdjacent(e, 0, adjVtx);
        APF_ITERATE(apf::Adjacent, adjVtx, v) {
          apf::Parts res;
         mesh->getResidence(*v,res);
          if( !res.count(dest) ) //not on the boundary with dest
            s.insert(*v);
        }
      }

      double weight(SetEnt& s) {
        double w = 0;
        APF_ITERATE(SetEnt, s, sItr)
          w += getWeight(*sItr);
        return w;
      }

      void cancel(apf::Migration** plan, parma::Mid* capacity) {
        typedef std::pair<apf::MeshEntity*, int> PairEntInt;
        apf::Migration* planA = *plan;
        std::vector<PairEntInt > keep; //temporary plan container
        const size_t planSz = TO_SIZET(planA->count());
        keep.reserve(planSz);

        //The plan is a vector so this loop will visit the plan's elements in
        //the same order that they were selected in, which, with graph distance
        //sorting, is from far to near.  An element is kept in the plan if
        //adding its vertices to the weight does not exceed the peer's vtx
        //weight capacity, capacity[dest].
        std::map<int,SetEnt > peerVerts;
        for(int i=0; i < planA->count(); i++) {
           apf::MeshEntity* e = planA->get(i);
           int dest = planA->sending(e);
           SetEnt vset = peerVerts[dest];
           insertInteriorVerts(e, dest, vset);
           if( weight(vset) <= (*capacity)[dest] ) {
             keep.push_back(PairEntInt(e,dest));
             peerVerts[dest] = vset;
           }
        }
        delete capacity;
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

      parma::Mid* trim(parma::Targets*, apf::Migration* plan) {
        //compute the weight of the vertices being sent to each peer
        typedef std::map<int,SetEnt > PeerEntSet;
        PeerEntSet peerVerts;
        for(int i=0; i < plan->count(); i++) {
          apf::MeshEntity* elm = plan->get(i);
          const int dest = plan->sending(elm);
          insertInteriorVerts(elm, dest, peerVerts[dest]);
        }

        parma::Mid sendingVtx;
        APF_ITERATE(PeerEntSet, peerVerts, pv)
          sendingVtx[pv->first] = weight(pv->second);

        typedef std::set<Migr,CompareMigr> MigrComm;

        PCU_Comm_Begin();
        APF_ITERATE(parma::Mid, sendingVtx, s)
          PCU_COMM_PACK(s->first, s->second);
        PCU_Comm_Send();

        MigrComm incoming; //map<sender,weight sending>
        double w;
        while (PCU_Comm_Listen()) {
          PCU_COMM_UNPACK(w);
          incoming.insert(Migr(PCU_Comm_Sender(),w));
        }

        double selfW = parma::getWeight(mesh,wtag,0);
        parma::Mid accept;
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
        APF_ITERATE(parma::Mid, accept, a)
          PCU_COMM_PACK(a->first, a->second);
        PCU_Comm_Send();
        parma::Mid* capacity = new parma::Mid;
        double outw;
        while (PCU_Comm_Listen()) {
          PCU_COMM_UNPACK(outw);
          (*capacity)[PCU_Comm_Sender()] = outw;
        }
        return capacity;
      }

      virtual double add(apf::MeshEntity*, apf::Up& cavity,
          const int destPid, apf::Migration* plan) {
        SetEnt cav;
        double w = 0;
        for(int i=0; i < cavity.n; i++) {
          plan->send(cavity.e[i], destPid);
          w += getWeight(cavity.e[i]);
        }
        return w;
      }
  };
}//end namespace

namespace parma {
  Selector* makeElmLtVtxSelector(apf::Mesh* m, apf::MeshTag* w, double maxVtx) {
    return new ElmLtVtxSelector(m, w, maxVtx);
  }
}
