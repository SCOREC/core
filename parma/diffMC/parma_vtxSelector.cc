#include "parma_selector.h"
#include "parma_targets.h"
#include "parma_weights.h"
#include "parma_commons.h"
#include <apf.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <PCU.h>
#include "../viz/viz.h"
#include "maximalIndependentSet/mis.h"

#include <set>
#include <list>
#include <stdio.h>
#include <limits.h>
#include <sstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <vector>

static int vtxSelectorCalls = 0;
Visualization* viz;

namespace {
  typedef std::pair<int,double> Migr;
  struct CompareMigr {
    bool operator()(const Migr& a, const Migr& b) {
      if( a.second < b.second )
        return true;
      else
        return false;
    }
  };
  typedef std::set<Migr,CompareMigr> MigrComm;
  void print(const char* key, MigrComm& d) {
    std::stringstream ss;
    ss << key << " ";
    APF_ITERATE(MigrComm, d, sItr)
      ss << "(" << sItr->first << "," << sItr->second << ") ";
    ss << '\n';
    std::string s = ss.str();
    PCU_Debug_Print(s.c_str());
  }

  typedef std::set<apf::MeshEntity*> SetEnt;
  typedef std::map<int,SetEnt > MiSetEnt;
  void print(const char* key, MiSetEnt& d) {
    std::stringstream ss;
    ss << key << " ";
    APF_ITERATE(MiSetEnt, d, sItr)
      ss << "(" << sItr->first << "," << sItr->second.size() << ") ";
    ss << '\n';
    std::string s = ss.str();
    PCU_Debug_Print(s.c_str());
  }

  typedef std::vector<int> VecInt;
  void print(const char* key, VecInt& d) {
    std::stringstream ss;
    ss << key << " ";
    APF_ITERATE(VecInt, d, sItr)
      ss << *sItr << " ";
    ss << '\n';
    std::string s = ss.str();
    PCU_Debug_Print(s.c_str());
  }

  typedef std::map<int,double> Mid;
  void print(const char* key, Mid& d) {
    std::stringstream ss;
    ss << key << " ";
    APF_ITERATE(Mid, d, sItr)
      ss << "(" << sItr->first << "," << sItr->second << ") ";
    ss << '\n';
    std::string s = ss.str();
    PCU_Debug_Print(s.c_str());
  }

  typedef std::map<int,int> Mii;
  void printMii(const char* key, Mii& d) {
    std::stringstream ss;
    ss << key << " ";
    APF_ITERATE(Mii, d, sItr)
      ss << "(" << sItr->first << "," << sItr->second << ") ";
    ss << '\n';
    std::string s = ss.str();
    PCU_Debug_Print(s.c_str());
  }

  typedef std::set<apf::MeshEntity*> Level;

  void writeVtk(apf::Mesh* m, const char* pre, int n=0) {
    std::stringstream ss;
    ss << pre << std::setfill('0') << std::setw(4) << vtxSelectorCalls << "_";
    if(n)
      ss << n << "_";
    std::string s = ss.str();
    apf::writeVtkFiles(s.c_str(), m);
  }

  inline void setNumber(apf::Numbering* n, apf::MeshEntity* e, int val) {
    int node = 0; int comp = 0;
    apf::number(n,e,node,comp,val);
  }
  inline int getNumber(apf::Numbering* n, apf::MeshEntity* e) {
    int node = 0; int comp = 0;
    return apf::getNumber(n,e,node,comp);
  }
  inline bool numbered(apf::Numbering* n, apf::MeshEntity* e) {
    int node = 0; int comp = 0;
    if ( 0 == apf::getNumber(n,e,node,comp) )
      return false;
    else
      return true;
  }

  apf::Numbering* initNumbering(apf::Mesh* m, const char* name, int initVal=0) {
    apf::FieldShape* s = m->getShape();
    apf::Numbering* n = apf::createNumbering(m,name,s,1);
    apf::MeshEntity* e;
    apf::MeshIterator* it = m->begin(0);
    while( (e = m->iterate(it)) )
      setNumber(n,e,initVal);
    m->end(it);
    return n;
  }

  apf::Numbering* initElmNumbering(apf::Mesh* m, const char* name, 
      int initVal=0) {
    const int dim = m->getDimension();
    apf::FieldShape* s = apf::getConstant(dim);
    apf::Numbering* n = apf::createNumbering(m,name,s,1);
    apf::MeshEntity* e;
    apf::MeshIterator* it = m->begin(dim);
    while( (e = m->iterate(it)) )
      setNumber(n,e,initVal);
    m->end(it);
    return n;
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

  int walkInward(apf::Numbering* n, apf::Mesh* m) {
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
        if( numbered(n, v) ) continue;
        setNumber(n,v,depth);
        count++;
        assert(count <= m->count(0));
        apf::Adjacent adjVtx;
        getEdgeAdjVtx(m,v,adjVtx);
        APF_ITERATE(apf::Adjacent, adjVtx, vItr) {
          if( !numbered(n, *vItr) )
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
  
  Level* reduce(apf::Numbering* n, apf::Mesh* m, Level* l, int depth) {
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
        APF_ITERATE(apf::Adjacent, adjVtx, u)
          if( getNumber(n,*u) == depth && l->count(*u) ) {
            l->erase(*u);
            q.insert(*u);
          }
      }
    }
    delete l;
    return g;
  }

  Level* getVerts(apf::Numbering* n, apf::Mesh* m, int depth) {
    Level* l = new Level;
    apf::MeshEntity* v;
    apf::MeshIterator* it = m->begin(0);
    while( (v = m->iterate(it)) ) {
      int d = getNumber(n,v);
      if( d == depth )
        l->insert(v);
    }
    m->end(it);
    l = reduce(n,m,l,depth);
    return l;
  }

  Level* getCentralVerts(apf::Mesh* m) {
    double t0 = PCU_Time();
    apf::Numbering* lvls = initNumbering(m, "parmaSelectorLevels");
    int depth = walkInward(lvls, m);
    Level* deepest = getVerts(lvls,m,depth);
    apf::destroyNumbering(lvls);
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

  struct CompareByDist {
    apf::Numbering* dist;
    bool isMore;
    bool operator()(apf::MeshEntity* u, apf::MeshEntity* v) {
      const int ud = getNumber(dist,u);
      const int vd = getNumber(dist,v);
      if( (isMore && ud > vd) || (!isMore && ud < vd) ) 
        return true;
      else
        return false;
    }
  };

  class Heap {
    apf::Mesh* m;
    CompareByDist cmp;
    apf::MeshTag* t;
    std::vector<apf::MeshEntity*> pq;
    public:
    Heap(apf::Mesh* mesh, CompareByDist c) : m(mesh), cmp(c) {
      t = m->createIntTag("parmaHeap", 1);
    }
    ~Heap() {
      apf::removeTagFromDimension(m, t, 0);
      m->destroyTag(t);
    }
    inline apf::MeshEntity* pop() {
      pop_heap(pq.begin(), pq.end(), cmp);
      apf::MeshEntity* v = pq.back();
      pq.pop_back();
      m->removeTag(v,t);
      return v;
    }
    inline void push(apf::MeshEntity* v) {
      pq.push_back(v);
      push_heap(pq.begin(), pq.end(), cmp);
      int one = 1;
      m->setIntTag(v,t,&one);
    }
    inline void heapify() {
      make_heap(pq.begin(), pq.end(), cmp);
    }
    inline bool has(apf::MeshEntity* v) {
      return m->hasTag(v,t);
    }
    inline bool empty() {
      return 0 == pq.size();
    }
  };

  void walkElms(apf::Mesh* m, apf::Numbering* d, apf::MeshEntity* src) {
    int count = 0;
    std::list<apf::MeshEntity*> elms;
    elms.push_back(src);
    while( !elms.empty() ) {
      apf::MeshEntity* e = elms.front();
      elms.pop_front();
      if( numbered(d,e) ) continue;
      setNumber(d,e,1);
      count++;
      apf::Adjacent adjElms;
      getFaceAdjElms(m,e,adjElms);
      APF_ITERATE(apf::Adjacent, adjElms, eItr)
        if( !numbered(d,*eItr) )
          elms.push_back(*eItr);
    }
  }

  bool disconnected(apf::Mesh*m, apf::Numbering* c, apf::MeshEntity* e) {
    apf::Adjacent adjElms;
    m->getAdjacent(e, m->getDimension(), adjElms);
    APF_ITERATE(apf::Adjacent, adjElms, adjItr) 
      if( getNumber(c,*adjItr) )
        return false;
    return true;
  }

  void dijkstra(apf::Mesh* m, apf::Numbering* c, apf::Numbering* d, 
      apf::MeshEntity* src) {
    CompareByDist compare = {d,true};
    Heap pq(m,compare);
    setNumber(d, src, 0);
    pq.push(src);

    while( !pq.empty() ) {
      apf::MeshEntity* v = pq.pop();
      int vd = getNumber(d, v);
      if( vd == INT_MAX) continue;
      apf::Adjacent adjVtx;
      getEdgeAdjVtx(m,v,adjVtx);
      APF_ITERATE(apf::Adjacent, adjVtx, eItr) {
        apf::MeshEntity* u = *eItr;
        int ud = getNumber(d,u);
        if( vd+1 < ud ) {
          int l = disconnected(m,c,u) ? INT_MAX : vd+1;
          setNumber(d,u,l);
          if( pq.has(u) )
            pq.heapify();
          else
            pq.push(u);
        }
      }
    }
  }

  apf::Numbering* computeDistance(apf::Mesh* m, Level& verts) {
    double t0 = PCU_Time();
    Level* centralElms = getCentralElms(m, verts);
    apf::Numbering* conn = initElmNumbering(m, "parmaElmConnectivity");
    APF_ITERATE(Level, *centralElms, itr)
      walkElms(m,conn,*itr);
    delete centralElms;
    parmaCommons::printElapsedTime("computeConnectivity", PCU_Time()-t0);

    t0 = PCU_Time();
    apf::Numbering* dist = initNumbering(m, "parmaDistance", INT_MAX);
    APF_ITERATE(Level, verts, itr)
      dijkstra(m,conn,dist,*itr);
    apf::destroyNumbering(conn);
    parmaCommons::printElapsedTime("computeDistance", PCU_Time()-t0);
    return dist;
  }

  Heap* BoundaryVertices(apf::Mesh* m, apf::Numbering* d) {
    CompareByDist compare = {d,false};
    Heap* dh = new Heap(m, compare);

    apf::MeshEntity* v;
    apf::MeshIterator* it = m->begin(0);
    while( (v = m->iterate(it)) )
      if( m->isShared(v) )
        dh->push(v);
    return dh;
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
}

namespace parma {
  class VtxSelector : public Selector {
    public:
      VtxSelector(apf::Mesh* m, apf::MeshTag* w)
        : Selector(m, w)
      {
        Level* centralVerts = getCentralVerts(m);
        dist = computeDistance(m, *centralVerts);
/*
        if(!vtxSelectorCalls)
          viz = new Visualization(4242); //FIXME no dtor call
*/
        vtxSelectorCalls++;
        delete centralVerts;
      }
      ~VtxSelector() {
        apf::destroyNumbering(dist);
      }
      apf::Migration* run(Targets* tgts) {
        apf::Migration* plan = new apf::Migration(mesh);
        select(tgts, plan);
        return plan;
      }
    protected:
      virtual double getWeight(apf::MeshEntity* e) {
        return getEntWeight(mesh,e,wtag);
      }
      virtual double add(apf::MeshEntity* v, const int destPid,
          apf::Migration* plan) {
        apf::Adjacent adjElms;
        mesh->getAdjacent(v, mesh->getDimension(), adjElms);
        APF_ITERATE(apf::Adjacent, adjElms, adjItr)
          if( !plan->has(*adjItr) )
            plan->send(*adjItr, destPid);
        return getWeight(v);
      }
      void select(Targets* tgts, apf::Migration* plan) {
        double t0 = PCU_Time();
        double planW = 0;
        Heap* bdryVerts = BoundaryVertices(mesh, dist);
        apf::Parts peers;
        while( !bdryVerts->empty() ) {
          if( planW > tgts->total() ) break;
          apf::MeshEntity* e = bdryVerts->pop();
          int destPid = getCavityPeer(mesh,e);
          if( (tgts->has(destPid) && sending[destPid] < tgts->get(destPid)) ||
              INT_MAX == getNumber(dist,e) ) {
            double ew = add(e, destPid, plan);
            sending[destPid] += ew;
            planW += ew;
          }
        }
        print("Sending", sending);
        parmaCommons::printElapsedTime("select", PCU_Time()-t0);
        delete bdryVerts;
      }
    private:
      apf::Numbering* dist;
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
      virtual double add(apf::MeshEntity* vtx, const int destPid, 
          apf::Migration* plan) {
        double w = 0;
        apf::DynamicArray<apf::MeshEntity*> adjElms;
        mesh->getAdjacent(vtx, mesh->getDimension(), adjElms);
        APF_ITERATE(apf::Adjacent, adjElms, adjItr)
          if( !plan->has(*adjItr) ) {
            plan->send(*adjItr, destPid);
            w += getWeight(*adjItr);
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
        select(tgts, plan);
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
        print("order", *order);
        PCU_Debug_Print("plan count %d\n", planA->count());
        typedef std::pair<apf::MeshEntity*, int> PairEntInt;
        std::vector<PairEntInt > keep;
        keep.reserve(planA->count());
        MiSetEnt sendingVtx;
        for(int i=0; i < planA->count(); i++) {
           apf::MeshEntity* e = planA->get(i);
           int dest = planA->sending(e);
           SetEnt vset = sendingVtx[dest];
           addCavityVtx(e, vset);
           if( cavityWeight(vset) <= (*order)[dest] ) {
             keep.push_back(PairEntInt(e,dest));
             sendingVtx[dest] = vset;
           }
        }
        print("sendingVtx", sendingVtx);
        delete order;
        delete planA;
        *plan = new apf::Migration(mesh);
        for(size_t i=0; i < keep.size(); i++)
          (*plan)->send(keep[i].first, keep[i].second);
        PCU_Debug_Print("plan count %d\n", (*plan)->count());
      }


      Mid* trim(Targets* t) {
        print("sendingVtx", sendingVtx);
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
        print("incoming", incoming);

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
        print("accept", accept);

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

      virtual double add(apf::MeshEntity* vtx,
          const int destPid, apf::Migration* plan) {
        SetEnt cav;
        double w = 0;
        apf::DynamicArray<apf::MeshEntity*> adjElms;
        mesh->getAdjacent(vtx, mesh->getDimension(), adjElms);
        APF_ITERATE(apf::Adjacent, adjElms, adjItr)
          if( !plan->has(*adjItr) ) {
            addCavityVtx(*adjItr, cav);
            plan->send(*adjItr, destPid);
            w += getWeight(*adjItr);
          }
        sendingVtx[destPid] += cavityWeight(cav);
        return w;
      }
    private:
      Mid sendingVtx;
      int maxVtx;
      apf::MeshTag* migrTag;
  };
  Selector* makeElmLtVtxSelector(apf::Mesh* m, apf::MeshTag* w, double maxVtx) {
    return new ElmLtVtxSelector(m, w, maxVtx);
  }
} //end namespace parma
