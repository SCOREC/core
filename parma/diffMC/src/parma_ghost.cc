#include <apfPartition.h>
#include <apf.h>
#include <PCU.h>
#include <stdio.h>
#include <map>
#include <assert.h>
#include <sstream>
#include <string>
#include <apfNumbering.h>
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
#include "parma_selector.h"

using parma::Sides;
using parma::Weights;
using parma::Targets;
using parma::Selector;
namespace {
  double runBFS(apf::Mesh* m, int layers, std::vector<apf::MeshEntity*> current,
      std::vector<apf::MeshEntity*> next, apf::MeshTag* visited,
      apf::MeshTag* wtag)
  {
    int yes=1;
    double weight = 0;
    apf::MeshEntity* checkVertex=NULL;
    for (unsigned int i=0;i<next.size();i++) {
      checkVertex=next[0];
      weight += getEntWeight(m,next[i],wtag);
      m->setIntTag(next[i],visited,&yes);
    }
    for (int i=1;i<=layers;i++) {
      for (unsigned int j=0;j<current.size();j++) {
        apf::MeshEntity* vertex = current[j];
        apf::Up edges;
        m->getUp(vertex,edges);
        for (int k=0;k<edges.n;k++) {
          apf::MeshEntity* v = getOtherVtx(m,edges.e[k],vertex);
          if (!m->isOwned(v))
            continue;
          if (m->hasTag(v,visited))
            continue;
          assert(v!=checkVertex);
          next.push_back(v);
          m->setIntTag(v,visited,&i);
          weight += getEntWeight(m,v,wtag);
        } 
      }
      current=next;
      next.clear();
    }
    return weight;
  }
 

  class GhostFinder {
    public:
      GhostFinder(apf::Mesh* m, apf::MeshTag* w, int l, int b) 
        : mesh(m), wtag(w), layers(l), bridge(b) {
        
        
      }
      /**
       * @brief get the weight of vertices ghosted to peer
       */
      double weight(int peer) {
        depth = mesh->createIntTag("depths",1);
        apf::MeshIterator* itr = mesh->begin(0);
        apf::MeshEntity* v;
        std::vector<apf::MeshEntity*> current;
        std::vector<apf::MeshEntity*> next;
        while ((v=mesh->iterate(itr))) {
          if (isSharedWithTarget(mesh,v,peer)) {
            if (mesh->isOwned(v))
              next.push_back(v);
            else if (mesh->getOwner(v)==peer)
              current.push_back(v);
          }
        }
        double weight = runBFS(mesh,layers,current,next,depth,wtag);
        renderIntTag(mesh,depth,"depth",peer);
        destroyTag(mesh,depth,0);
        return weight;

      }
    private:
      GhostFinder();
      apf::Mesh* mesh;
      apf::MeshTag* wtag;
      int layers;
      int bridge;
      apf::MeshTag* depth;
  };


  class Ghosts : public Associative<double> {
    public:
      Ghosts(GhostFinder* finder, Sides* sides) {
        weight = 0;
        init(finder, sides);
        exchangeGhostsFrom();
        exchangeGhostsTo();
      }
      double self() {
        return weight;
      }
    private:
      double weight;
      void init(GhostFinder* finder, Sides* sides) {
        const Sides::Item* side;
        sides->begin();
        while( (side = sides->iterate()) )
          set(side->first, finder->weight(side->first));
        sides->end(); 
      }
    void exchangeGhostsFrom() {
      //ghosts from this part
      PCU_Comm_Begin();
      const Ghosts::Item* ghost;
      begin();
      while( (ghost = iterate()) ) 
        PCU_COMM_PACK(ghost->first, ghost->second);
      end();
      PCU_Comm_Send();
      while (PCU_Comm_Listen()) {
        double otherWeight;
        PCU_COMM_UNPACK(otherWeight);
        weight += otherWeight;
      }
    }
    void exchangeGhostsTo() {
      //all elements ghosted to this part
      PCU_Comm_Begin();
      const Ghosts::Item* ghost;
      begin();
      while( (ghost = iterate()) ) 
        PCU_COMM_PACK(ghost->first, weight);
      end();
      PCU_Comm_Send();
      while (PCU_Comm_Listen()) {
        double peerGhostWeight = 0;
        PCU_COMM_UNPACK(peerGhostWeight);
        set(PCU_Comm_Sender(), peerGhostWeight);
      }
    }

  class GhostTargets : public Targets {
    public:
      GhostTargets(Sides* s, Weights* w, Ghosts* g, double alpha) 
        : Targets(s, w, alpha)
      {
        init(s, w, g, alpha);
      }
      double total() {
        return totW;
      }
    private:
      GhostTargets();
      double totW;
      void init(Sides* s, Weights* w, Ghosts* g, double alpha) {
        totW = 0;
        const Sides::Item* side;
        s->begin();
        while( (side = s->iterate()) ) {
          const int peer = side->first;
          const double selfW = w->self() + g->self();
          const double peerW = g->get(peer) + w->get(peer);
          if ( selfW > peerW ) {
            const double difference = selfW - peerW;
            double sideFraction = side->second;
            sideFraction /= s->total();
            double scaledW = difference * sideFraction * alpha;
            set(peer, scaledW);
            totW+=scaledW;
          }
        }
        s->end();
      }
  };

  class ParmaGhost {
    public:
      ParmaGhost(apf::Mesh* mIn, apf::MeshTag* wIn, 
          int layersIn, int bridgeIn, double alphaIn) 
        : m(mIn), w(wIn), layers(layersIn), bridge(bridgeIn), alpha(alphaIn)
      {
        sides = parma::makeElmBdrySides(m);
        weights = parma::makeGhostWeights(m,w,sides);
        targets = new GhostTargets(sides, weights, ghosts, alpha);
        selects = parma::makeVtxSelector(m, w);
	iters=0;
      }

      ~ParmaGhost();
      bool run(double maxImb);
    private:
      ParmaGhost();
      apf::Mesh* m;
      apf::MeshTag* w;
      int layers;
      int bridge;
      double alpha;
      int verbose;
      double imbalance();
      Sides* sides;
      Weights* weights;
      GhostFinder* ghostFinder;
      Ghosts* ghosts;
      Targets* targets;
      Selector* selects;
      int iters;
  };

  ParmaGhost::~ParmaGhost() {
    delete sides;
    delete weights;
    delete ghostFinder;
    delete ghosts;
    delete targets;
    delete selects;
  }

  bool ParmaGhost::run(double maxImb) {
    const double imb = imbalance();
    if ( 0 == PCU_Comm_Self() )
      fprintf(stdout, "imbalance %.3f\n", imb);
    if ( imb < maxImb) 
      return false;
    apf::Migration* plan = selects->run(targets);
    PCU_Debug_Print("plan size %d\n",plan->count());
    m->migrate(plan);
    return true;
  }

  double ParmaGhost::imbalance() { 
    double maxWeight = 0, totalWeight = 0;
    sides->print("sides");
    targets->print("tgts");
    ghosts->print("ghosts");
    const double w = weights->self();
    const double gw = ghosts->self();
    PCU_Debug_Print("w %.3f gw %.3f\n", w, gw);
    maxWeight = totalWeight = w + gw;
    PCU_Add_Doubles(&totalWeight,1);
    PCU_Max_Doubles(&maxWeight, 1);
    double averageWeight = totalWeight / PCU_Comm_Peers();
    return maxWeight / averageWeight;
  }
}

class GhostBalancer : public apf::Balancer {
  public:
    GhostBalancer(apf::Mesh* m, int l, int b, double f, int v)
      : mesh(m), factor(f), layers(l), bridge(b), verbose(v) {
        (void) verbose; // silence!
    }
    bool runStep(apf::MeshTag* weights, double tolerance) {
        parma::Balancer ghost(mesh, wtag, factor);
        Sides* s = makeElmBdrySides(mesh);
        ghost.setSides(s);
        Weights* w = makeEntWeights(mesh, wtag, s, mesh->getDimension());
        ghost.setWeights(w);
        Targets* t = makeTargets(s, w, factor);
        ghost.setTargets(t);
        ghost.setSelector(makeVtxSelector(mesh, wtag));
        return ghost.run(tolerance);

      ParmaGhost ghost(mesh, weights, layers, bridge, factor);
      return ghost.run(tolerance);
    }
    virtual void balance(apf::MeshTag* weights, double tolerance) {
      double t0 = MPI_Wtime(); int iters=0;
      while (runStep(weights,tolerance)&&iters<20) {iters++;}
      double t1 = MPI_Wtime();
      if (!PCU_Comm_Self())
        printf("ghost balanced to %f in %f seconds\n", tolerance, t1-t0);
      if (!PCU_Comm_Self())
        printf("number of iterations %d\n", iters);
    }
  private:
    apf::Mesh* mesh;
    double factor;
    int layers;
    int bridge;
    int verbose;
};

apf::Balancer* Parma_MakeGhostDiffuser(apf::Mesh* m, 
    int layers, int bridge, double stepFactor, int verbosity) {
  return new GhostBalancer(m, layers, bridge, stepFactor, verbosity);
}
