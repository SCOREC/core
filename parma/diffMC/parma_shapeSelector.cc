#include <PCU.h>
#include "parma_selector.h"
#include "parma_targets.h"
#include "parma_weights.h"
#include "parma_centroids.h"
#include <apf.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <float.h>

#include <sstream>
#include <string.h>

namespace {
  apf::MeshEntity* getOtherElement(apf::Mesh* m, 
      apf::MeshEntity* e, apf::MeshEntity* s) {
    apf::Up up;
    m->getUp(s,up);
    for (int i=0; i < up.n; ++i) {
      apf::MeshEntity* o = up.e[i];
      if (o != e)
        return o;
    }
    return 0;
  }

  int getAdjacentElements(apf::Mesh* m, apf::MeshEntity* e, apf::Downward oe) {
    int n = 0;
    apf::Downward s;
    int ns = m->getDownward(e,m->getDimension()-1,s);
    for (int i=0; i < ns; ++i) {
      oe[n] = getOtherElement(m, e, s[i]);
      if (oe[n]) ++n;
    }
    return n;
  }

  class Distance {
    public:
      struct Closest {
        int peer;
        double dist;
      };
      Distance(apf::Mesh* mesh, parma::Centroids* centroids, 
          parma::Targets* tgts) : m(mesh), c(centroids), t(tgts) {}
      ~Distance() {}
      double get(apf::MeshEntity* e, int to) {
        return (apf::getLinearCentroid(m, e) - c->get(to)).getLength();
      }
      double get(apf::MeshEntity* e) {
        return (apf::getLinearCentroid(m, e) - c->self()).getLength();
      }
      Closest closest(apf::MeshEntity* elm) {
        if( !t->size() )
          PCU_Debug_Print("t->size() %lu\n", t->size());
        assert(t->size());
        Closest cl = {-1, DBL_MAX};
        const parma::Targets::Item* tgt;
        t->begin();
        while( (tgt = t->iterate()) ) {
          double d = get(elm, tgt->first);
          if( d < cl.dist ) {
            cl.peer = tgt->first;
            cl.dist = d;
          }
        }
        t->end();
        assert(cl.peer != -1);
        return cl;
      }
    private:
      apf::Mesh* m;
      parma::Centroids* c;
      parma::Targets* t;
  };

  void setNumber(apf::Numbering* n, apf::MeshEntity* e, int val) {
    int node = 0;
    apf::number(n,e,node,0,val);
  }

  apf::Numbering* initNumbering(apf::Mesh* m, const char* name) {
    apf::FieldShape* s = apf::getConstant(m->getDimension());
    apf::Numbering* n = apf::createNumbering(m,name,s,1);

    int negOne = 42;
    apf::MeshEntity* e;
    apf::MeshIterator* it = m->begin(m->getDimension());
    while( (e = m->iterate(it)) ) 
      setNumber(n,e,negOne);
    m->end(it);
    return n;
  }

  class DistanceQueue {
    typedef std::multimap< double, 
      apf::MeshEntity*, 
      std::greater<double> > DistanceQ;
    public:
      DistanceQueue(apf::Mesh* mesh, apf::MeshTag* tag, Distance* dist, apf::Numbering* popN) 
        : m(mesh), t(tag), d(dist), pt(popN) {
          popCnt = 0;
      }
      ~DistanceQueue() {}
      void push(apf::MeshEntity* e) {
        if ( m->hasTag(e, t) )
          return;
        Distance::Closest c = d->closest(e);
        double dist = d->get(e);
        m->setIntTag(e, t, &(c.peer));
        q.insert(std::make_pair(dist, e));
      }
      apf::MeshEntity* pop() {
        DistanceQ::iterator it = q.begin();
        apf::MeshEntity* e = it->second;
        q.erase(it);
        setNumber(pt,e,popCnt++);
        return e;
      }
      bool empty() {
        return q.empty();
      }
    private:
      apf::Mesh* m;
      apf::MeshTag* t;
      Distance* d;
      DistanceQ q;
      apf::Numbering* pt;
      int popCnt;
  };
}

namespace parma {
  class ShapeSelector : public Selector {
    public:
      ShapeSelector(apf::Mesh* m, apf::MeshTag* w, Centroids* c)
        : Selector(m, w), centroids(c), wTag(w) {
      }
      apf::Migration* run(Targets* tgts) {
        apf::Migration* plan = new apf::Migration(mesh);
        sendTag = mesh->createIntTag("centroid_send",1);
        Distance dist(mesh, centroids, tgts);
        apf::Numbering* popNum = initNumbering(mesh, "pop_order");
        DistanceQueue distQ(mesh, sendTag, &dist, popNum);
        init(tgts, &distQ, &dist);
        while ( !distQ.empty() ) 
          trySending(tgts, &distQ, plan);
        clearSendTag();
	//        writeVtk(mesh, "pop");
        apf::destroyNumbering(popNum);
        writeSending();
        return plan;
      }
    private:
      void writeSending() {
        std::stringstream ss;
        typedef std::map<int,double> mid;
        ss << "sending: ";
        APF_ITERATE(mid, sending, send) 
          ss << send->first << ',' << send->second << " ";
        ss << '\n';
        std::string s = ss.str();
        PCU_Debug_Print("%s",s.c_str());
      }
      Centroids* centroids;
      std::map<int,double> sending;
      apf::MeshTag* sendTag;
      apf::MeshTag* wTag;
      ShapeSelector();
      ~ShapeSelector() {}
      void clearSendTag() {
        apf::removeTagFromDimension(mesh,sendTag,mesh->getDimension());
        mesh->destroyTag(sendTag);
      }
      int getTaggedSend(apf::MeshEntity* e) {
        int to;
        mesh->getIntTag(e,sendTag,&to);
        return to;
      }
      void trySending(Targets* tgts, DistanceQueue* distQ, 
          apf::Migration* plan) {
        apf::MeshEntity* e = distQ->pop();
        assert( !plan->has(e) );
        int peer = getTaggedSend(e);
        if (sending[peer] >= tgts->get(peer))
          return;
        plan->send(e,peer);
        sending[peer] += getEntWeight(mesh, e, wTag);
        apf::Downward oe;
        int noe = getAdjacentElements(mesh, e, oe);
        for (int i=0; i < noe; ++i)
          distQ->push(oe[i]);
      }
      void init(Targets* tgts, DistanceQueue* distQ, Distance* dist) {
        if( !tgts->size() ) 
          return;
        apf::MeshEntity* v;
        apf::MeshIterator* it = mesh->begin(0);
        double maxDist = DBL_MIN;
        apf::MeshEntity* maxElm = NULL;
        while( (v = mesh->iterate(it)) ) {
          if ( hasTargets(v, tgts) ) {
            apf::Adjacent elms;
            mesh->getAdjacent(v, mesh->getDimension(), elms);
            double d = dist->get(elms[0]);
            if( d > maxDist ) {
              maxElm = elms[0];
              maxDist = d; 
            }
          }
        }
        mesh->end(it);
        if( maxElm )
          distQ->push(maxElm);
      }
      bool hasTargets(apf::MeshEntity* e, parma::Targets* tgts) {
        apf::Parts res;
        mesh->getResidence(e, res);
        size_t found = 0;
        tgts->begin();
        const parma::Targets::Item* tgt;
        while( (tgt = tgts->iterate()) )
          found += res.count(tgt->first);
        tgts->end();
        return (found == tgts->size());
      }
  };
  Selector* makeShapeSelector(apf::Mesh* m, apf::MeshTag* w, Centroids* c) {
    return new ShapeSelector(m, w, c);
  }
} //end namespace parma
