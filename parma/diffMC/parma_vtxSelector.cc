#include "parma_selector.h"
#include "parma_targets.h"
#include "parma_weights.h"
#include "parma_commons.h"
#include <apf.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <PCU.h>
#include <set>
#include <list>
#include <stdio.h>
#include <limits.h>
#include <sstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <vector>

namespace {

  typedef std::set<apf::MeshEntity*> SetEnt;
  typedef std::vector<int> VecInt;
  typedef std::map<int,double> Mid;
  typedef std::map<int,int> Mii;
  typedef std::set<apf::MeshEntity*> Level;

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

  inline void getEdgeAdjVtx(apf::Mesh* m, apf::MeshEntity* v,
      apf::Adjacent& adj) {
    int bridge = 1; int tgt = 0;
    getBridgeAdjacent(m, v, bridge, tgt, adj);
  }

  inline void getFaceAdjElms(apf::Mesh* m, apf::MeshEntity* e,
      apf::Adjacent& adj) {
    int bridge = m->getDimension()-1; int tgt = m->getDimension();
    getBridgeAdjacent(m, e, bridge, tgt, adj);
  }

  bool onBoundary(apf::Mesh* m, apf::MeshEntity* e) {
     int gd = m->getModelType(m->toModel(e));
     int md = m->getDimension();
     bool shared = m->isShared(e);
     return shared || (!shared && gd < md);
  }

  int walkInward(apf::MeshTag* n, apf::Mesh* m) {
    Level cur;
    Level next;
    apf::MeshEntity* v;
    apf::MeshIterator* it = m->begin(0);
    while( (v = m->iterate(it)) )
      if( onBoundary(m,v) )
        cur.insert(v);
    m->end(it);

    size_t count = 0;
    int depth = 1;
    while( count != m->count(0) ) {
      APF_ITERATE(Level, cur, vtxItr) {
        v = *vtxItr;
        int lvl; m->getIntTag(v,n,&lvl);
        if( lvl ) continue;
        m->setIntTag(v,n,&depth);
        count++;
        assert(count <= m->count(0));
        apf::Adjacent adjVtx;
        getEdgeAdjVtx(m,v,adjVtx);
        APF_ITERATE(apf::Adjacent, adjVtx, vItr) {
          m->getIntTag(*vItr,n,&lvl);
          if( !lvl )
            next.insert(*vItr);
        }
      }
      cur = next;
      assert(cur.size() == cur.size());
      next.clear();
      depth++;
    }
    return depth-1;
  }

  Level* reduce(apf::MeshTag* n, apf::Mesh* m, Level* l, int depth) {
    Level* g = new Level;
    while( !l->empty() ) {
      apf::MeshEntity* v = *(l->begin()); l->erase(v);
      g->insert(v);
      Level q;
      q.insert(v);
      while( !q.empty() ) {
        v = *(q.begin()); q.erase(v);
        apf::Adjacent adjVtx;
        getEdgeAdjVtx(m,v,adjVtx);
        APF_ITERATE(apf::Adjacent, adjVtx, u) {
          int d; m->getIntTag(*u,n,&d);
          if( d == depth && l->count(*u) ) {
            l->erase(*u);
            q.insert(*u);
          }
        }
      }
    }
    delete l;
    return g;
  }

  Level* getVerts(apf::MeshTag* n, apf::Mesh* m, int depth) {
    Level* l = new Level;
    apf::MeshEntity* v;
    apf::MeshIterator* it = m->begin(0);
    while( (v = m->iterate(it)) ) {
      int d; m->getIntTag(v,n,&d);
      if( d == depth )
        l->insert(v);
    }
    m->end(it);
    l = reduce(n,m,l,depth);
    return l;
  }

  Level* getCentralVerts(apf::Mesh* m) {
    double t0 = PCU_Time();
    apf::MeshTag* lvlsT = initTag(m, "parmaSelectorLevels");
    int depth = walkInward(lvlsT, m);
    Level* deepest = getVerts(lvlsT,m,depth);
    //apf::destroyNumbering(lvls);
    apf::removeTagFromDimension(m,lvlsT,0);
    m->destroyTag(lvlsT);
    parmaCommons::printElapsedTime("getCentralVerts", PCU_Time()-t0);
    return deepest;
  }

  Level* getCentralElms(apf::Mesh* m, Level& verts) {
    Level* elms = new Level;
    APF_ITERATE(Level, verts, itr) {
      apf::Adjacent adjElms;
      m->getAdjacent(*itr, m->getDimension(), adjElms);
      APF_ITERATE(apf::Adjacent, adjElms, adjItr)
        elms->insert(*adjItr);
    }
    return elms;
  }

  struct Greater {
    bool operator() (const int& l, const int& r) const {
      return (l > r);
    }
  };
  struct Less {
    bool operator() (const int& l, const int& r) const {
      return (l < r);
    }
  };

  template <class Compare> class DistanceQueue {
    typedef typename std::multimap<int, apf::MeshEntity*, Compare> DistanceQ;
    typedef typename DistanceQ::iterator DistanceQIter;
    public:
      DistanceQueue(apf::Mesh* mesh) : m(mesh) {
        t = m->createIntTag("parmaDistanceQueue",1);
      }
      ~DistanceQueue() {
        apf::removeTagFromDimension(m,t,0);
        m->destroyTag(t);
      }
      void push(apf::MeshEntity* e, int dist) {
        DistanceQIter it = q.begin();
        if ( m->hasTag(e, t) ) {
          it = erase(dist, e);
        }
        int one = 1;
        m->setIntTag(e, t, &one);
        q.insert(it, std::make_pair(dist, e));
      }
      apf::MeshEntity* pop() {
        DistanceQIter it;
        it = q.begin();
        apf::MeshEntity* e = it->second;
        q.erase(it);
        return e;
      }
      bool empty() {
        return q.empty();
      }
    private:
      DistanceQIter erase(int dist, apf::MeshEntity* e) {
        assert( m->hasTag(e, t) );
        DistanceQIter it = q.find(dist);
        DistanceQIter rit = ( it != q.begin() ) ? it-- : q.end();
        while( it != q.end() ) {
          if( it->second == e ) {
            q.erase(it);
            break;
          }
          rit = it;
          it++;
        }
        return rit;
      }
      apf::Mesh* m;
      apf::MeshTag* t;
      DistanceQ q;
  };

  void walkElms(apf::Mesh* m, apf::MeshTag* conn, apf::MeshEntity* src) {
    int one = 1;
    int count = 0;
    std::list<apf::MeshEntity*> elms;
    elms.push_back(src);
    while( !elms.empty() ) {
      apf::MeshEntity* e = elms.front();
      elms.pop_front();
      int c; m->getIntTag(e,conn,&c);
      if( c ) continue;
      m->setIntTag(e,conn,&one);
      count++;
      apf::Adjacent adjElms;
      getFaceAdjElms(m,e,adjElms);
      APF_ITERATE(apf::Adjacent, adjElms, eItr) {
        m->getIntTag(*eItr,conn,&c);
        if( !c )
          elms.push_back(*eItr);
      }
    }
  }

  bool disconnected(apf::Mesh*m, apf::MeshTag* conn, apf::MeshEntity* e) {
    int c;
    apf::Adjacent adjElms;
    m->getAdjacent(e, m->getDimension(), adjElms);
    APF_ITERATE(apf::Adjacent, adjElms, adjItr) {
      m->getIntTag(*adjItr,conn,&c);
      if( c ) return false;
    }
    return true;
  }

  void dijkstra(apf::Mesh* m, apf::MeshTag* c, apf::MeshTag* d,
      apf::MeshEntity* src) {
    DistanceQueue<Less> pq(m);
    int zero = 0;
    m->setIntTag(src, d, &zero);
    pq.push(src,0);

    while( !pq.empty() ) {
      apf::MeshEntity* v = pq.pop();
      int vd; m->getIntTag(v, d, &vd);
      if( vd == INT_MAX) continue;
      apf::Adjacent adjVtx;
      getEdgeAdjVtx(m,v,adjVtx);
      APF_ITERATE(apf::Adjacent, adjVtx, eItr) {
        apf::MeshEntity* u = *eItr;
        int ud; m->getIntTag(u,d,&ud);
        if( vd+1 < ud ) {
          int l = disconnected(m,c,u) ? INT_MAX : vd+1;
          m->setIntTag(u,d,&l);
          pq.push(u,l);
        }
      }
    }
  }

  apf::MeshTag* computeDistance(apf::Mesh* m, Level& verts) {
    double t0 = PCU_Time();
    Level* centralElms = getCentralElms(m, verts);
    int initVal = 0;
    apf::MeshTag* connT =
      initTag(m, "parmaElmConnectivity", initVal, m->getDimension());
    APF_ITERATE(Level, *centralElms, itr)
      walkElms(m,connT,*itr);
    delete centralElms;
    parmaCommons::printElapsedTime("computeConnectivity", PCU_Time()-t0);

    t0 = PCU_Time();
    initVal = INT_MAX;
    apf::MeshTag* distT = initTag(m, "parmaDistance", initVal);

    APF_ITERATE(Level, verts, itr)
      dijkstra(m,connT,distT,*itr);

    apf::removeTagFromDimension(m,connT,m->getDimension());
    m->destroyTag(connT);

    parmaCommons::printElapsedTime("computeDistance", PCU_Time()-t0);
    return distT;
  }

  DistanceQueue<Greater> * BoundaryVertices(apf::Mesh* m, apf::MeshTag* d) {
    DistanceQueue<Greater> * dq = new DistanceQueue<Greater>(m);

    apf::MeshEntity* v;
    apf::MeshIterator* it = m->begin(0);
    while( (v = m->iterate(it)) )
      if( m->isShared(v) ) {
        int dist; m->getIntTag(v,d,&dist);
        dq->push(v, dist);
      }
    m->end(it);
    return dq;
  }

  int getCavityPeer(apf::Mesh* m, apf::MeshEntity* v) {
    Mii pc;
    apf::Adjacent sideSides;
    m->getAdjacent(v, m->getDimension()-2, sideSides);
    APF_ITERATE(apf::Adjacent, sideSides, ss) {
      apf::Copies rmts;
      m->getRemotes(*ss,rmts);
      APF_ITERATE(apf::Copies, rmts, r)
         pc[r->first]++;
    }
    int max = -1;
    int peer = -1;
    APF_ITERATE(Mii, pc, p)
      if( p->second > max ) {
         peer = p->first;
         max = p->second;
      }
    assert(peer>=0);
    return peer;
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
        Level* centralVerts = getCentralVerts(m);
        dist = computeDistance(m, *centralVerts);
        delete centralVerts;
      }
      ~VtxSelector() {
        apf::removeTagFromDimension(mesh,dist,0);
        mesh->destroyTag(dist);
      }
      apf::Migration* run(Targets* tgts) {
        apf::Migration* plan = new apf::Migration(mesh);
        double planW = 0;
        for(size_t max=2; max <= 12; max+=2)
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
        DistanceQueue<Greater>* bdryVerts = BoundaryVertices(mesh, dist);
        apf::Parts peers;
        apf::Up cavity;
        while( !bdryVerts->empty() ) {
          if( planW > tgts->total() ) break;
          apf::MeshEntity* e = bdryVerts->pop();
          getCavity(mesh, e, plan, cavity);
          int destPid = getCavityPeer(mesh,e);
          int d; mesh->getIntTag(e,dist,&d);
          if( (tgts->has(destPid) &&
               sending[destPid] < tgts->get(destPid) &&
               cavity.n <= maxSize ) ||
              INT_MAX == d ) {
            double ew = add(e, cavity, destPid, plan);
            sending[destPid] += ew;
            planW += ew;
          }
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
      double getWeight(apf::MeshEntity* vtx) {
        apf::Up edges;
        mesh->getUp(vtx, edges);
        double w = 0;
        for(int i=0; i<edges.n; i++)
          if( mesh->isShared(edges.e[i]) )
            w += getWeight(edges.e[i]);
        return w;
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
        apf::DynamicArray<apf::MeshEntity*> adjVtx;
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
        PCU_Debug_Print("plan count %d\n", planA->count());
        typedef std::pair<apf::MeshEntity*, int> PairEntInt;
        std::vector<PairEntInt > keep;
        keep.reserve(planA->count());

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
        PCU_Debug_Print("plan count %d\n", (*plan)->count());
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
          PCU_Debug_Print("trim sending to %d weight %.3f\n",
              s->first, s->second);
        }
        PCU_Comm_Send();

        MigrComm incoming;
        double w;
        while (PCU_Comm_Listen()) {
          PCU_COMM_UNPACK(w);
          PCU_Debug_Print("trim recv from %d weight %.3f\n", PCU_Comm_Sender(), w);
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
        PCU_Debug_Print("trim selfW %.3f\n", selfW);

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
      int maxVtx;
  };
  Selector* makeElmLtVtxSelector(apf::Mesh* m, apf::MeshTag* w, double maxVtx) {
    return new ElmLtVtxSelector(m, w, maxVtx);
  }
} //end namespace parma
