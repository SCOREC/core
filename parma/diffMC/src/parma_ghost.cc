#include <apfPartition.h>
#include <PCU.h>
#include <stdio.h>
#include <map>
#include <assert.h>
#include <sstream>
#include <string>
#include <apfNumbering.h>

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
      bool has(int key) {
        return (c.count(key) != 0);
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
        PCU_Debug_Print("%s", str.c_str());
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
   *         of the voroni mesh the vertex weights will correspond to element 
   *         weights of the voroni mesh
   */
  double getWeight(apf::Mesh* m, apf::MeshTag* w) {
    apf::MeshIterator* it = m->begin(0);
    apf::MeshEntity* e;
    double sum = 0;
    while ((e = m->iterate(it)))
      if (m->isOwned(e))
        sum += getEntWeight(m, e, w);
    m->end(it);
    return sum;
  }

  class Weights : public Associative<double> {
    public:
      Weights(apf::Mesh* m, apf::MeshTag* w, Sides* s) {
        selfWeight = getWeight(m, w);
        init(s);
      }
      double self() {
        return selfWeight;
      }
    private:
      Weights();
      double selfWeight;
      void init(Sides* s) {
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

  void renderIntTag(apf::Mesh* m, apf::MeshTag* tag,
      const char* filename, int peer)
  {
    apf::Numbering* n = apf::createNumbering(m,
        m->getTagName(tag), m->getShape(), 1);
    apf::MeshIterator* it = m->begin(0);
    apf::MeshEntity* v;
    while ((v = m->iterate(it))) {
      if (m->hasTag(v, tag)) {
        int x;
        m->getIntTag(v, tag, &x);
        apf::number(n, v, 0, 0, x);
      } else {
        apf::number(n, v, 0, 0, 0);
      }
    }
    m->end(it);
    std::stringstream ss;
    ss << filename
      << '_' << PCU_Comm_Self()
      << '_' << peer;
    std::string s = ss.str();
    apf::writeOneVtkFile(s.c_str(), m);
    apf::destroyNumbering(n);
  }

  void runBFS(apf::Mesh* m, int layers, std::vector<apf::MeshEntity*> current,
      std::vector<apf::MeshEntity*> next, apf::MeshTag* visited) {
    int yes=1;
    for (unsigned int i=0;i<next.size();i++) {
      m->setIntTag(next[i],visited,&yes);
    }
    for (int i=1;i<=layers;i++) {
      for (unsigned int j=0;j<current.size();j++) {
        apf::MeshEntity* vertex = current[j];
        if (!m->isOwned(vertex))
          continue;
        apf::Up edges;
        m->getUp(vertex,edges);
        for (int k=0;k<edges.n;k++) {
          apf::MeshEntity* edge = edges.e[k];
          apf::Downward vertices;
          int nvertices = m->getDownward(edge,getDimension(m,edge)-1,vertices);
          assert(nvertices==2);
          apf::MeshEntity* v = vertices[0];
          if (v==vertex)
            v=vertices[1];
          if (m->hasTag(v,visited)) continue;
          next.push_back(v);
          m->setIntTag(current[i],visited,&i);
        }
        current=next;
        next.clear();
      }
    }
  }

  bool searchBFS(apf::Mesh* m,apf::MeshEntity* start, int target, int layers) {
    apf::MeshTag* visited = m->createIntTag("visited",1);
    std::vector<apf::MeshEntity*> current;
    std::vector<apf::MeshEntity*> next;
    current.push_back(start);
    int yes;
    m->setIntTag(start,visited,&yes);
    for (int i=0;i<layers;i++) {
      for (unsigned int j=0;j<current.size();j++) {
        apf::MeshEntity* vertex = current[j];
        apf::Up edges;
        m->getUp(vertex,edges);
        for (int k=0;k<edges.n;k++) {
          apf::MeshEntity* edge = edges.e[k];
          apf::Downward vertices;
          int nvertices = m->getDownward(edge,getDimension(m,edge)-1,vertices);
          assert(nvertices==2);
          apf::MeshEntity* v = vertices[0];
          if (v==vertex)
            v=vertices[1];
          if (m->getOwner(v)==target) {
            m->destroyTag(visited);
            return true;
          }
          if (m->hasTag(v,visited)) continue;
          next.push_back(v);
          m->setIntTag(v,visited,&yes);
        }
        
      }
      current=next;
      next.clear();
    }
    m->destroyTag(visited);
    return false;
  }

  class GhostFinder {
    public:
      GhostFinder(apf::Mesh* m, apf::MeshTag* w, int l, int b) 
        : mesh(m), wtag(w), layers(l), bridge(b) {
        depth = m->createIntTag("depths",1);
        apf::MeshIterator* itr = mesh->begin(0);
        apf::MeshEntity* v;
        std::vector<apf::MeshEntity*> current;
        std::vector<apf::MeshEntity*> next;
        while ((v=mesh->iterate(itr))) {
          if (mesh->isShared(v)&&mesh->isOwned(v))
            next.push_back(v);
          if (!mesh->isOwned(v))
            current.push_back(v);
        }
        runBFS(mesh,layers,current,next,depth);
      }
      /**
       * @brief get the weight of vertices ghosted to peer
       */
      double weight(int peer) {
        double totalWeight=0;
        apf::MeshIterator* itr = mesh->begin(0);
        apf::MeshEntity* v;
        while ((v=mesh->iterate(itr))) {
          if (mesh->hasTag(v,depth)&&searchBFS(mesh,v,peer,layers)) {
            assert(mesh->hasTag(v,wtag));
            double w;
            mesh->getDoubleTag(v,wtag,&w);
            totalWeight+=w;
          }
        }
        return totalWeight;
      }
    ~GhostFinder() {
      mesh->destroyTag(depth);
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
        exchange();
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
      double total() {
        return totW;
      }
    private:
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

  class Selector : public Associative<double> {
    public:
      Selector(apf::Mesh* m, apf::MeshTag* w, Targets* t) 
        : mesh(m), tgts(t), wtag(w) {}
      apf::Migration* run() {
        apf::Migration* plan = new apf::Migration(mesh);
        vtag = mesh->createIntTag("ghost_visited",1);
        const int maxBoundedElm = 6;
        double planW=0;
        for( size_t maxAdjElm=2; maxAdjElm<=maxBoundedElm; maxAdjElm+=2)
          planW += select(planW, maxAdjElm, plan);
        apf::removeTagFromDimension(mesh,vtag,0);
        mesh->destroyTag(vtag);
        return plan;
      }
    private:
      apf::Mesh* mesh;
      Targets* tgts;
      apf::MeshTag* vtag;
      apf::MeshTag* wtag;
      Selector();
      double add(apf::MeshEntity* vtx, const size_t maxAdjElm, 
          const int destPid, apf::Migration* plan) {
        apf::DynamicArray<apf::MeshEntity*> adjElms;
        mesh->getAdjacent(vtx, mesh->getDimension(), adjElms);
        if( adjElms.getSize() > maxAdjElm ) 
          return 0;
        double w = getEntWeight(mesh,vtx,wtag);
        for(size_t i=0; i<adjElms.getSize(); i++) {
          apf::MeshEntity* elm = adjElms[i];
          if ( mesh->hasTag(elm, vtag) ) continue;
          mesh->setIntTag(elm, vtag, &destPid); 
          plan->send(elm, destPid);
          set(destPid, get(destPid)+w);
        }
        return w;
      }
      double select(const double planW, 
          const size_t maxAdjElm, 
          apf::Migration* plan) {
        double planWeight = 0;
        apf::MeshEntity* vtx;
        apf::MeshIterator* itr = mesh->begin(0);
        while( (vtx = mesh->iterate(itr)) && 
               (planW + planWeight < tgts->total()) ) {
          apf::Copies rmt;
          mesh->getRemotes(vtx, rmt);
          if( 1 == rmt.size() ) {
            int destPid = (rmt.begin())->first;
            if( tgts->has(destPid) && get(destPid) < tgts->get(destPid) )
              planWeight += add(vtx, maxAdjElm, destPid, plan);
          }
        }
        mesh->end(itr);
        return planWeight;
      }
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
        ghosts = new Ghosts(ghostFinder, sides);
        targets = new Targets(sides, weights, ghosts, alpha);
        targets->print("tgts");
        selects = new Selector(m, w, targets); 
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
    if ( imb < maxImb ) 
      return false;
    apf::Migration* plan = selects->run();
    m->migrate(plan);
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
        (void) verbose; // silence!
    }
    bool runStep(apf::MeshTag* weights, double tolerance) {
      parma::ParmaGhost ghost(mesh, weights, layers, bridge, factor);
      return ghost.run(tolerance);
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
