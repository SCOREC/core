#include <apf.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <stdio.h>
#include <map>
#include <assert.h>

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
      GhostFinder(apf::Mesh* m, apf::Mesh* weight, int layer, int bridge);
      /**
       * @brief get the weight of vertices ghosted to peer
       */
      double weight(int peer);
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
        init(m, finder, sides);
      }
      double total() {
        double total=0;
        const Ghosts::Item* ghost;
        begin();
        while( (ghost = iterate()) )
          total += ghost->second;
        end();
        return total;
      }
    private:
      void init(apf::Mesh* m, GhostFinder* finder, Sides* sides) {
        const Sides::Item* side;
        sides->begin();
        while( (side = sides->iterate()) )
          set(side->first, finder->weight(side->first));
        sides->end();
      }
  };

  class ParmaGhost {
    public:
      ParmaGhost(apf::Mesh* m, apf::MeshTag* w, int layers, int bridgeDim);
      ~ParmaGhost();
      void run(double maxImb, int verbosity);
    private:
      ParmaGhost();
      apf::Mesh* m;
      double maxImb;
      apf::MeshTag* w;
      int layers;
      int bridge;
      int verbose;
      double imbalance();
      void exchangeGhosts();
      void diffuse();
      Sides* sides;
      Weights* weights;
      GhostFinder* ghostFinder;
      Ghosts* ghosts;
  };

  ParmaGhost::ParmaGhost(apf::Mesh* mIn, apf::MeshTag* wIn, 
      int layersIn, int bridgeIn) 
    : m(mIn), w(wIn), layers(layersIn), bridge(bridgeIn)
  {
    sides = new Sides(m);
    weights = new Weights(m, w, sides);
  }

  ParmaGhost::~ParmaGhost() {
    delete sides;
    delete weights;
    delete ghostFinder;
  }

  void ParmaGhost::run(double maxImb, int verbosityIn) {
    if ( imbalance() < maxImb ) return;
  }

  double ParmaGhost::imbalance() { //FIXME
    double totalWeight = weights->self() + ghosts->total();
    PCU_Add_Doubles(&totalWeight,1);
    double averageWeight = totalWeight / PCU_Comm_Peers();
    double maxWeight = weights->self();
    PCU_Max_Doubles(&maxWeight, 1);
    return maxWeight / averageWeight;
  }
}
