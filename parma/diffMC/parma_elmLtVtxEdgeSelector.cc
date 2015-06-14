#include <set>
#include <apf.h>
#include <PCU.h>
#include "parma_vtxSelector.h"
#include "parma_targets.h"
#include "parma_weights.h"
#include "parma_convert.h"

namespace {
  struct dtwo {
    dtwo(double x, double y)
      : a(x), b(y) {}
    dtwo()
      : a(0), b(0) {}
    double a;
    double b;
  };
  typedef std::map<int, dtwo> Midd;
  typedef std::set<apf::MeshEntity*> SetEnt;

  class ElmLtVtxEdgeSelector : public parma::VtxSelector {
    private:
      double maxVtx;
      double maxEdge;

    public:
      ElmLtVtxEdgeSelector(apf::Mesh* m, apf::MeshTag* w, double maxV, double maxE)
        : VtxSelector(m, w), maxVtx(maxV), maxEdge(maxE) { }

      apf::Migration* run(parma::Targets* tgts) {
        apf::Migration* plan = new apf::Migration(mesh);
        double planW = 0;
        for(int max=2; max <= 12; max+=2)
          planW += select(tgts, plan, planW, max);
        Midd* capacity = trim(tgts,plan);
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

      void insertInteriorEdges(apf::MeshEntity* e, int dest, SetEnt& s) {
        apf::Adjacent adjEdges;
        mesh->getAdjacent(e, 1, adjEdges);
        APF_ITERATE(apf::Adjacent, adjEdges, edge) {
          apf::Parts res;
          mesh->getResidence(*edge,res);
          if( !res.count(dest) ) //not on the boundary with dest
            s.insert(*edge);
        }
      }

      double weight(SetEnt& s) {
        double w = 0;
        APF_ITERATE(SetEnt, s, sItr)
          w += getWeight(*sItr);
        return w;
      }

      void cancel(apf::Migration** plan, Midd* capacity) {
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
        std::map<int,SetEnt > peerEdges;
        for(int i=0; i < planA->count(); i++) {
           apf::MeshEntity* e = planA->get(i);
           int dest = planA->sending(e);
           SetEnt vset = peerVerts[dest];
           SetEnt eset = peerEdges[dest];
           insertInteriorVerts(e, dest, vset);
           insertInteriorEdges(e, dest, eset);
           if( weight(vset) <= (*capacity)[dest].a &&
               weight(eset) <= (*capacity)[dest].b ) {
             keep.push_back(PairEntInt(e,dest));
             peerVerts[dest] = vset;
             peerEdges[dest] = eset;
           }
        }
        delete capacity;
        delete planA;
        *plan = new apf::Migration(mesh);
        for(size_t i=0; i < keep.size(); i++)
          (*plan)->send(keep[i].first, keep[i].second);
      }

      struct Migr {
        Migr(int i, double v, double e) :
          id(i), vw(v), ew(e) {}
        int id;
        double vw;
        double ew;
      };
      struct CompareMigr {
        //sort by vertex weight... because
        bool operator()(const Migr& a, const Migr& b) const {
          if( a.vw < b.vw )
            return true;
          else
            return false;
        }
      };
      typedef std::set<Migr,CompareMigr> MigrComm;

      //return  map<int neighbor, pair<double vtxW, double edgeW> > where
      //  neighbor is a neighbor's part id
      //  vtxW is the vtx weight capacity of neighbor
      //  edgeW is the edge weight capacity of neighbor
      Midd* trim(parma::Targets*, apf::Migration* plan) {
        //compute the weight of the vertices and edges being sent to each peer
        typedef std::map<int,SetEnt > PeerEntSet;
        PeerEntSet peerVerts;
        PeerEntSet peerEdges;
        for(int i=0; i < plan->count(); i++) {
          apf::MeshEntity* elm = plan->get(i);
          const int dest = plan->sending(elm);
          insertInteriorVerts(elm, dest, peerVerts[dest]);
          insertInteriorEdges(elm, dest, peerEdges[dest]);
        }

        //send vtx and edge weight
        PCU_Comm_Begin();
        APF_ITERATE(PeerEntSet, peerVerts, pv) {
          const int dest = pv->first;
          double vw = weight(peerVerts[dest]);
          PCU_COMM_PACK(dest, vw);
          double ew = weight(peerEdges[dest]);
          PCU_COMM_PACK(dest, ew);
        }
        PCU_Comm_Send();
        MigrComm incoming;
        double vw, ew;
        while (PCU_Comm_Listen()) {
          PCU_COMM_UNPACK(vw);
          PCU_COMM_UNPACK(ew);
          incoming.insert(Migr(PCU_Comm_Sender(),vw,ew)); //Migr ctor does not exist
        }

        Midd accept;
        double totW[2] = {
          parma::getWeight(mesh,wtag,0), //vtx
          parma::getWeight(mesh,wtag,1)  //edge
        };
        APF_ITERATE(MigrComm, incoming, in) {
          const double inW[2] = {
            (*in).vw,  //vtx
            (*in).ew   //edge
          };
          double avail[2] = {
            maxVtx - totW[0], //vtx
            maxEdge - totW[1] //edge
          };
          const int nbor = (*in).id;
          if( avail[0] > 0 && avail[1] > 0 ) {
            if( inW[0] <= avail[0] && inW[1] <= avail[1] ) {
              accept[nbor] = dtwo(inW[0],inW[1]);
            } else {
              accept[nbor] = dtwo(avail[0],avail[1]);
            }
          } else {
            accept[nbor] = dtwo(0,0);
          }
          totW[0] += accept[nbor].a;
          totW[1] += accept[nbor].b;
        }

        PCU_Comm_Begin();
        APF_ITERATE(Midd, accept, acc) {
          PCU_COMM_PACK(acc->first, acc->second.a);
          PCU_COMM_PACK(acc->first, acc->second.b);
        }
        PCU_Comm_Send();
        Midd* capacity = new Midd;
        dtwo outw;
        while (PCU_Comm_Listen()) {
          PCU_COMM_UNPACK(outw.a);
          PCU_COMM_UNPACK(outw.b);
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
  Selector* makeElmLtVtxEdgeSelector(apf::Mesh* m, apf::MeshTag* w, double maxVtx, double maxEdge) {
    return new ElmLtVtxEdgeSelector(m, w, maxVtx, maxEdge);
  }
}
