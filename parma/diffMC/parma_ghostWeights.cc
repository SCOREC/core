#include <assert.h>
#include <PCU.h>
#include <apf.h>
#include <apfMesh.h>
#include <apfNumbering.h>
#include "parma_weights.h"
#include "parma_sides.h"
#include "parma_ghostOwner.h"

namespace {
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
      std::vector<apf::MeshEntity*> next, apf::MeshTag* visited,
      apf::MeshTag* wtag)
  {
    int yes=1;
    double weight = 0;
    apf::MeshEntity* checkVertex=NULL;
    for (unsigned int i=0;i<next.size();i++) {
      checkVertex=next[0];
      weight += parma::getEntWeight(m,next[i],wtag);
      m->setIntTag(next[i],visited,&yes);
    }
    for (int i=1;i<=layers;i++) {
      for (unsigned int j=0;j<current.size();j++) {
        apf::MeshEntity* vertex = current[j];
        apf::Up edges;
        m->getUp(vertex,edges);
        for (int k=0;k<edges.n;k++) {
          apf::MeshEntity* v = getOtherVtx(m,edges.e[k],vertex);
          if (!parma::isOwned(m, v))
            continue;
          if (m->hasTag(v,visited))
            continue;
          assert(v!=checkVertex);
          next.push_back(v);
          m->setIntTag(v,visited,&i);
          weight += parma::getEntWeight(m,v,wtag);
        }
      }
      current=next;
      next.clear();
    }
    return weight;
  }

  double ownedVtxWeight(apf::Mesh* m, apf::MeshTag* w) {
    apf::MeshIterator* it = m->begin(0);
    apf::MeshEntity* e;
    double entW = 0;
    double sum = 0;
    while ((e = m->iterate(it))) {
      assert(m->hasTag(e,w));
      if (parma::isOwned(m,e)) {
        m->getDoubleTag(e,w,&entW);
        sum += entW;
      }
    }
    m->end(it);
    return sum;
  }
}

namespace parma {
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
            if (isOwned(mesh,v))
              next.push_back(v);
            else if (getOwner(mesh,v)==peer)
              current.push_back(v);
          }
        }
        double weight = runBFS(mesh,layers,current,next,depth,wtag);
        apf::removeTagFromDimension(mesh,depth,0);
        mesh->destroyTag(depth);
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



  class GhostWeights : public Weights {
    public:
      GhostWeights(apf::Mesh* m, apf::MeshTag* wtag, Sides* s, int layers, int bridge)
        : Weights(m, wtag, s), entDim(0), weight(0)
      {
        GhostFinder finder(m, wtag, layers, bridge);
        findGhosts(&finder, s);
        exchangeGhostsFrom();
        weight += ownedVtxWeight(m, wtag);
        exchange();
      }
      ~GhostWeights() {}
      double self() {
        return weight;
      }
    private:
      GhostWeights();
      int entDim;
      double weight;
      void findGhosts(GhostFinder* finder, Sides* sides) {
        const Sides::Item* side;
        sides->begin();
        while( (side = sides->iterate()) )
          set(side->first, finder->weight(side->first));
        sides->end();
      }
      void exchangeGhostsFrom() {
        PCU_Comm_Begin();
        const GhostWeights::Item* ghost;
        begin();
        while( (ghost = iterate()) )
          PCU_COMM_PACK(ghost->first, ghost->second);
        end();
        PCU_Comm_Send();
        while (PCU_Comm_Listen()) {
          double ghostsFromPeer = 0;
          PCU_COMM_UNPACK(ghostsFromPeer);
          weight += ghostsFromPeer;
        }
      }
      void exchange() {
        PCU_Comm_Begin();
        const GhostWeights::Item* ghost;
        begin();
        while( (ghost = iterate()) )
          PCU_COMM_PACK(ghost->first, weight);
        end();
        PCU_Comm_Send();
        while (PCU_Comm_Listen()) {
          double peerWeight;
          PCU_COMM_UNPACK(peerWeight);
          int peer = PCU_Comm_Sender();
          set(peer, peerWeight);
        }
      }
  };
  Weights* makeGhostWeights(apf::Mesh* m, apf::MeshTag* w, Sides* s,
      int layers, int bridge) {
    return new GhostWeights(m, w, s, layers, bridge);
  }
} //end namespace
