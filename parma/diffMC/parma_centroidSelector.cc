#include "PCU.h"
#include "parma_selector.h"
#include "parma_targets.h"
#include "parma_weights.h"
#include "parma_centroids.h"
#include "apf.h"

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
      Distance(apf::Mesh* mesh, parma::Centroids* centroids) 
        : m(mesh), c(centroids) {
        }
      ~Distance() {}
      double get(apf::MeshEntity* e, int to) {
        return (apf::getLinearCentroid(m, e) - c->get(to)).getLength();
      }
    private:
      apf::Mesh* m;
      parma::Centroids* c;
  };

  class DistanceQueue {
    typedef std::multimap<double,apf::MeshEntity*> DistanceQ;
    public:
      DistanceQueue(apf::Mesh* mesh, apf::MeshTag* tag, Distance* dist) 
        : m(mesh), t(tag), d(dist) {
      }
      ~DistanceQueue() {}
      void push(apf::MeshEntity* e, int to) {
        if ( m->hasTag(e, t))
          return;
        double dist = d->get(e, to);
        m->setIntTag(e, t, &to);
        q.insert(std::make_pair(dist, e));
      }
      apf::MeshEntity* pop() {
        DistanceQ::iterator it = q.begin();
        apf::MeshEntity* e = it->second;
        q.erase(it);
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
  };
}

namespace parma {
  class CentroidSelector : public Selector {
    public:
      CentroidSelector(apf::Mesh* m, apf::MeshTag* w, Centroids* c)
        : Selector(m, w), centroids(c), wTag(w) {
      }
      apf::Migration* run(Targets* tgts) {
        apf::Migration* plan = new apf::Migration(mesh);
        sendTag = mesh->createIntTag("centroid_send",1);
        Distance dist(mesh, centroids);
        DistanceQueue distQ(mesh, sendTag, &dist);
        init(tgts, &distQ, &dist);
        while ( ! distQ.empty()) 
          trySending(tgts, &distQ, plan);
        clearSendTag();
        return plan;
      }
    private:
      Centroids* centroids;
      std::map<int,double> sending;
      apf::MeshTag* sendTag;
      apf::MeshTag* wTag;
      CentroidSelector();
      ~CentroidSelector() {}
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
        assert( ! plan->has(e));
        int peer = getTaggedSend(e);
        if (sending[peer] >= tgts->get(peer))
          return;
        plan->send(e,peer);
        sending[peer] += getEntWeight(mesh, e, wTag);
        apf::Downward oe;
        int noe = getAdjacentElements(mesh, e, oe);
        for (int i=0; i < noe; ++i)
          distQ->push(oe[i],peer);
      }
      void init(Targets* tgts, DistanceQueue* distQ, Distance* dist) {
        apf::MeshEntity* s;
        apf::MeshIterator* it = mesh->begin(mesh->getDimension()-1);
        typedef std::map<int,std::pair<double,apf::MeshEntity*> > ClosestMap;
        ClosestMap closest;
        while ((s = mesh->iterate(it))) {
          if (isSide(s)) {
            int peer = getPeer(s);
            if ( ! tgts->has(peer))
              continue;
            apf::MeshEntity* elm = mesh->getUpward(s,0);
            double d = dist->get(elm,peer);
            if (( ! closest.count(peer)) ||
                (d < closest[peer].first))
              closest[peer] = std::make_pair(d,elm);
          }
        }
        mesh->end(it);
        APF_ITERATE(ClosestMap,closest,itr)
          distQ->push(itr->second.second,itr->first);
      }
      int getPeer(apf::MeshEntity* e) {
        return apf::getOtherCopy(mesh,e).peer;
      }
      bool isSide(apf::MeshEntity* e) {
        apf::Parts p;
        mesh->getResidence(e, p);
        return (p.size()==2);
      }
  };
  Selector* makeCentroidSelector(apf::Mesh* m, apf::MeshTag* w, Centroids* c) {
    return new CentroidSelector(m, w, c);
  }
} //end namespace parma
