#include <apfPartition.h>
#include <PCU.h>
#include <stdio.h>
#include <map>
#include <assert.h>
#include <sstream>
#include <string>

namespace parma {
  template <class T> class Associative {
    typedef std::map<int, T> Container;
    public:
      Associative() {
        iteratorActive = false;
      }
      void begin() {
        assert(!iteratorActive);
        iteratorActive = true;
        cItr = c.begin();
      }
      typedef std::pair<const int, T> Item;
      const Item* iterate() {
        assert(iteratorActive);
        if( cItr == c.end() ) 
          return NULL;
        else
          return &(*cItr++); // there is no spoon ... and this is likely crap
      }
      void end() {
        assert(iteratorActive);
        iteratorActive = false;
      }
      T get(int key) {
        return c[key];
      }
      void set(int key, T value) {
        c[key] = value;
      }
      void print(const char* key) {
        std::stringstream s;
        s << key << " ";
        const Item* i;
        begin();
        while( (i = iterate()) ) 
          s << i->first << " " << i->second << " ";
        end();
        std::string str = s.str();
        PCU_Debug_Print(str.c_str());
      }
    protected:
      Container c;
    private:
      typename Container::iterator cItr;
      bool iteratorActive;
  };


  class Sides : public Associative<int> {
    public:
      Sides(apf::Mesh* m) {
        totalSides = 0;
        init(m);
      }
      int total() {
        return totalSides;
      }
    private:
      int totalSides;
      void init(apf::Mesh* m) {
        apf::MeshEntity* s;
        apf::MeshIterator* it = m->begin(m->getDimension()-1);
        totalSides = 0;
        while ((s = m->iterate(it)))
          if (m->countUpward(s)==1 && m->isShared(s)) {
            const int peerId = apf::getOtherCopy(m,s).first;
            set(peerId, get(peerId)+1);
            ++totalSides;
          }
        m->end(it);
      }
  };

  double getEntWeight(apf::Mesh* m, apf::MeshEntity* e, apf::MeshTag* w) {
    assert(m->hasTag(e,w));
    double weight;
    m->getDoubleTag(e,w,&weight);
    return weight;
  }

  /**
   * @brief compute the weight of the mesh vertices
   * @remark since the mesh being partitioned is the delauney triangularization 
   *         of the voroni mesh so the vertex weights will correspond to element 
   *         weights of the voroni mesh
   */
  double getWeight(apf::Mesh* m, apf::MeshTag* w) {
    apf::MeshIterator* it = m->begin(0);
    apf::MeshEntity* e;
    double sum = 0;
    while ((e = m->iterate(it)))
      sum += getEntWeight(m, e, w);
    m->end(it);
    return sum;
  }

  class Weights : public Associative<double> {
    public:
      Weights(apf::Mesh* m, apf::MeshTag* w, Sides* s) {
        selfWeight = getWeight(m, w);
        init(m, w, s);
      }
      double self() {
        return selfWeight;
      }
    private:
      Weights();
      double selfWeight;
      void init(apf::Mesh* m, apf::MeshTag* w, Sides* s) {
        PCU_Comm_Begin();
        const Sides::Item* side;
        s->begin();
        while( (side = s->iterate()) ) 
          PCU_COMM_PACK(side->first, selfWeight);
        s->end();
        PCU_Comm_Send();
        while (PCU_Comm_Listen()) {
          double otherWeight;
          PCU_COMM_UNPACK(otherWeight);
          set(PCU_Comm_Sender(), otherWeight);
        }
      }
  };

  class GhostFinder {
    public:
      GhostFinder(apf::Mesh* m, apf::MeshTag* weight, int layer, int bridge) {
      }
      /**
       * @brief get the weight of vertices ghosted to peer
       */
      double weight(int peer) {
        return 0;
      }
    private:
      GhostFinder();
      int layers;
      int bridge;
      apf::Mesh* m;
      apf::MeshTag* w;
      apf::MeshTag* depth;
  };

  class Ghosts : public Associative<double> {
    public:
      Ghosts(apf::Mesh* m, GhostFinder* finder, Sides* sides) {
        weight = 0;
        init(m, finder, sides);
        exchange();
      }
      double self() {
        return weight;
      }
    private:
      double weight;
      void init(apf::Mesh* m, GhostFinder* finder, Sides* sides) {
        const Sides::Item* side;
        sides->begin();
        while( (side = sides->iterate()) )
          set(side->first, finder->weight(side->first));
        sides->end();
      }
      void exchange() {
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
  };

  class Targets : public Associative<double> {
    public:
      Targets(Sides* s, Weights* w, Ghosts* g, double alpha) {
        init(s, w, g, alpha);
      }
    private:
      void init(Sides* s, Weights* w, Ghosts* g, double alpha) {
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
            set(peer, difference * sideFraction * alpha);
          }
        }
        s->end();
      }
  };

  class Selector {
    public:
      Selector(apf::Mesh* m, Targets* tgts) {
      }
      void run() {
      }
    private:
      Selector();
  };

  class ParmaGhost {
    public:
      ParmaGhost(apf::Mesh* mIn, apf::MeshTag* wIn, 
          int layersIn, int bridgeIn, double alphaIn) 
        : m(mIn), w(wIn), layers(layersIn), bridge(bridgeIn), alpha(alphaIn)
      {
        sides = new Sides(m);
        weights = new Weights(m, w, sides);
        ghostFinder = new GhostFinder(m, w, layers, bridge);
        ghosts = new Ghosts(m, ghostFinder, sides);
        targets = new Targets(sides, weights, ghosts, alpha);
        targets->print("tgts");
        selects = new Selector(m, targets); //FIXME take targets and run selection
      }

      ~ParmaGhost();
      bool run(double maxImb, int verbosity);
    private:
      ParmaGhost();
      apf::Mesh* m;
      double maxImb;
      apf::MeshTag* w;
      int layers;
      int bridge;
      double alpha;
      int verbose;
      double imbalance();
      void exchangeGhosts();
      Sides* sides;
      Weights* weights;
      GhostFinder* ghostFinder;
      Ghosts* ghosts;
      Targets* targets;
      Selector* selects;
  };

  ParmaGhost::~ParmaGhost() {
    delete sides;
    delete weights;
    delete ghostFinder;
    delete ghosts;
  }

  bool ParmaGhost::run(double maxImb, int verbosityIn) {
    if ( imbalance() < maxImb ) 
      return false;
    selects->run();
    return true;
  }

  double ParmaGhost::imbalance() { 
    double maxWeight = 0, totalWeight = 0;
    maxWeight = totalWeight = weights->self() + ghosts->self();
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
        (void) factor; //suppress compiler warning
    }
    bool runStep(apf::MeshTag* weights, double tolerance) {
      const double alpha = 0.1;
      parma::ParmaGhost ghost(mesh, weights, layers, bridge, alpha);
      return ghost.run(tolerance, verbose);
    }
    virtual void balance(apf::MeshTag* weights, double tolerance) {
      double t0 = MPI_Wtime();
      while (runStep(weights,tolerance));
      double t1 = MPI_Wtime();
      if (!PCU_Comm_Self())
        printf("ghost balanced to %f in %f seconds\n", tolerance, t1-t0);
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
