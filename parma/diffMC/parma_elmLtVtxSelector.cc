#include <set>
#include <apf.h>
#include <PCU.h>
#include "parma_vtxSelector.h"
#include "parma_targets.h"
#include "parma_weights.h"

namespace {
  typedef std::set<apf::MeshEntity*> SetEnt;

  class ElmLtVtxSelector : public parma::VtxSelector {
    private:
      parma::Mid sendingVtx;
      double maxVtx;

    public:
      ElmLtVtxSelector(apf::Mesh* m, apf::MeshTag* w, double maxV)
        : VtxSelector(m, w), maxVtx(maxV) {}

      apf::Migration* run(parma::Targets* tgts) {
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

      void cancel(apf::Migration** plan, parma::Mid* order) {
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

      parma::Mid* trim(parma::Targets*) {
        typedef std::set<Migr,CompareMigr> MigrComm;

        PCU_Comm_Begin();
        APF_ITERATE(parma::Mid, sendingVtx, s) {
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
        parma::Mid* order = new parma::Mid;
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
  };
}//end namespace

namespace parma {
  Selector* makeElmLtVtxSelector(apf::Mesh* m, apf::MeshTag* w, double maxVtx) {
    return new ElmLtVtxSelector(m, w, maxVtx);
  }
}
