#include <assert.h>
#include <PCU.h>
#include "parma_entWeights.h"
#include "parma_sides.h"

namespace {
  double getEntWeight(apf::Mesh* m, apf::MeshEntity* e, apf::MeshTag* w) {
    assert(m->hasTag(e,w));
    double weight;
    m->getDoubleTag(e,w,&weight);
    return weight;
  }


  apf::MeshEntity* getOtherVtx(apf::Mesh* m, 
      apf::MeshEntity* edge, apf::MeshEntity* vtx) {
    apf::Downward dwnVtx;
    int nDwnVtx = m->getDownward(edge,getDimension(m,edge)-1,dwnVtx);
    assert(nDwnVtx==2);
    return (dwnVtx[0] != vtx) ? dwnVtx[0] : dwnVtx[1];
  }


  bool isSharedWithTarget(apf::Mesh* m,apf::MeshEntity* v, int target) {
    if( ! m->isShared(v) ) return false;
    apf::Copies rmts;
    m->getRemotes(v,rmts);
    APF_ITERATE(apf::Copies, rmts, itr) 
      if (itr->first==target) 
        return true;
    return false;
  }


  double runBFS(apf::Mesh* m, int layers, std::vector<apf::MeshEntity*> current,
      std::vector<apf::MeshEntity*> next, apf::MeshTag* visited,apf::MeshTag* wtag) {
    int yes=1;
    double weight;
    apf::MeshEntity* checkVertex=NULL;
    for (unsigned int i=0;i<next.size();i++) {
      checkVertex=next[0];
      weight+=getEntWeight(m,next[i],wtag);
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
          if (m->hasTag(v,visited)) continue;
          assert(v!=checkVertex);
          next.push_back(v);
          m->setIntTag(v,visited,&i);
          weight+=getEntWeight(m,v,wtag);
        } 
      }
      current=next;
      next.clear();
    }
    return weight;
  }


  class Ghosts : public parma::Associative<double> {
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
        apf::removeTagFromDimension(mesh,depth,0);
        m->destroyTag(depth);
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


}

namespace parma {  
  GhostWeights::GhostWeights(apf::Mesh* m, apf::MeshTag* w, Sides* s, int d) 
    : Weights(m, w, s), entDim(d) 
  {
    assert(entDim >= 0 && entDim <= 3);
    weight = selfWeight(m, w);
    init(m, w, s);
  }
  double GhostWeights::self() {
    return weight;
  }
  double GhostWeights::getEntWeight(apf::Mesh* m, apf::MeshEntity* e, apf::MeshTag* w) {
    assert(m->hasTag(e,w));
    double entW = 0;
    m->getDoubleTag(e,w,&entW);
    return entW;
  }
  double GhostWeights::selfWeight(apf::Mesh* m, apf::MeshTag* w) {
    apf::MeshIterator* it = m->begin(entDim);
    apf::MeshEntity* e;
    double sum = 0;
    while ((e = m->iterate(it)))
      sum += getEntWeight(m, e, w);
    m->end(it);
    return sum;
  }
  void GhostWeights::init(apf::Mesh* m, apf::MeshTag* w, Sides* s) {
    //find from ghosts
    const Sides::Item* side;
    sides->begin();
    while( (side = sides->iterate()) )
      set(side->first, finder->weight(side->first));
    sides->end(); 
    //exchange from ghosts
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
    //exchange to ghosts

    PCU_Comm_Begin();
    const Sides::Item* side;
    s->begin();
    while( (side = s->iterate()) ) 
      PCU_COMM_PACK(side->first, weight);
    s->end();
    PCU_Comm_Send();
    while (PCU_Comm_Listen()) {
      double otherWeight;
      PCU_COMM_UNPACK(otherWeight);
      set(PCU_Comm_Sender(), otherWeight);
    }
  }
  Weights* makeGhostWeights(apf::Mesh* m, apf::MeshTag* w, Sides* s, 
      int dim, int layers, int bridge) {
    return new GhostWeights(m, w, s, dim, layers, bridge);
  }
} //end namespace
