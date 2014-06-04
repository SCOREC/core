#include <stdio.h>
#include <map>
#include <assert.h>

namespace parma {
  template <class T> class Associative {
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
    private:
      typedef std::map<int, T> Container;
      Container c;
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
        apf::MeshIterator* it = m->begin(dim-1);
        totalSides = 0;
        while ((s = m->iterate(it)))
          if (m->countUpward(s)==1 && m->isShared(s)) {
            const int peerId = apf::getOtherCopy(m,s).first;
            c.set(peerId, c.get(peerId)+1);
            ++totalSides;
          }
        m->end(it);
      }
  };

  double getEntWeight(apf::MeshEntity* e, apf::MeshTag* w) {
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
      sum += getEntWeight(e, w);
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
        Sides::Item* side;
        s.begin();
        while( (side = s.iterate()) ) 
          PCU_COMM_PACK(side->first, selfWeight);
        s.end();
        PCU_Comm_Send();
        while (PCU_Comm_Listen()) {
          double otherWeight;
          PCU_COMM_UNPACK(otherWeight);
          c.set(PCU_Comm_Sender(), otherWeight);
        }
      }
  };

  class Ghosts : public Associative<double> {
    public:
      Ghosts(apf::Mesh* m, GhostFinder* finder, Sides* sides) {
        init(m, finder, sides)
      }
    private:
      void init(apf::Mesh* m, GhostFinder* finder, Sides* sides) {
        Sides::Item* side;
        sides.begin();
        while( (side = sides.iterate()) )
          c.set(side->first, finder->weight(side->first));
        sides.end();
      }
  };

  class GhostFinder {
    public:
      GhostFinder(apf::Mesh* m, apf::Mesh* weight, int depth, int bridge);
      /**
       * @brief get the weight of vertices ghosted to peer
       */
      double weight(int peer);
    private:
      GhostFinder();
      int depth;
      int bridge;
      apf::Mesh* m;
      apf::MeshTag* w;
      apf::MeshTag* depth;
  };

  class ParmaGhost {
    public:
      ParmaGhost(apf::Mesh* m, apf::MeshTag* w, int layers, int bridgeDim);
      ~ParmaGhost();
      void run(double maxImb, int verbosity);
    private:
      ParmaGhost();
      double maxImb;
      int verbose;
      int bridge;
      int layers;
      apf::Mesh* m;
      void getImbalance();
      void exchangeGhosts();
      void diffuse();
      Sides* sides;
      Weights* weights;
      GhostFinder* ghostFinder;
  };

  ParmaGhost::ParmaGhost(apf::Mesh* mIn, apf::MeshTag* w, 
      int layersIn, int bridgeIn) 
    : m(mIn), layers(layersIn), bridge(bridgeIn)
  {
    sides = new Sides(m);
    weights = new Weights(m, w, &sides);
  }

  ParmaGhost::~ParmaGhost() {
    delete sides;
    delete weights;
    delete ghostFinder;
  }

  void ParmaGhost::run(double maxImb, int verbosityIn) {
    if ( imbalance() < maxImb) return;

  }

  double ParmaGhost::imbalance() { //FIXME
    double totalWeight = weights.self() + ghosts;
    PCU_Add_Doubles(&totalWeight,1);
    double averageWeight = totalWeight / PCU_Comm_Peers();
    double maxWeight = selfWeight;
    PCU_Max_Doubles(&maxWeight, 1);
    return maxWeight / averageWeight;
  }
}
