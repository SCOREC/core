#include <set>
#include <apf.h>
#include "parma_vtxSelector.h"
#include "parma_targets.h"
#include "parma_weights.h"
#include "parma_convert.h"

namespace {
  typedef std::set<apf::MeshEntity*> SetEnt;

  class EdgeEqVtx : public parma::VtxSelector {
    private:
      double maxVtx;

    public:
      EdgeEqVtx(apf::Mesh* m, apf::MeshTag* w, double maxV)
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

        mesh->getPCU()->Begin();
        APF_ITERATE(parma::Mid, sendingVtx, s)
          mesh->getPCU()->Pack(s->first, s->second);
        mesh->getPCU()->Send();

        MigrComm incoming; //map<sender,weight sending>
        double w;
        while (mesh->getPCU()->Listen()) {
          mesh->getPCU()->Unpack(w);
          incoming.insert(Migr(mesh->getPCU()->Sender(),w));
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

        mesh->getPCU()->Begin();
        APF_ITERATE(parma::Mid, accept, a)
          mesh->getPCU()->Pack(a->first, a->second);
        mesh->getPCU()->Send();
        parma::Mid* capacity = new parma::Mid;
        double outw;
        while (mesh->getPCU()->Listen()) {
          mesh->getPCU()->Unpack(outw);
          (*capacity)[mesh->getPCU()->Sender()] = outw;
        }
        return capacity;
      }
  };
}//end namespace

namespace parma {
  Selector* makeEdgeEqVtxSelector(apf::Mesh* m, apf::MeshTag* w, double maxVtx) {
    return new EdgeEqVtx(m, w, maxVtx);
  }
}
