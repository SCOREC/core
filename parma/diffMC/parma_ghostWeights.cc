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

  // vertex based BFS
  // return an array of vertex, edge and element weights
  // - no one needs faces... yet
  // the returned array needs to be deallocated
  double* runBFS(apf::Mesh* m, int layers, std::vector<apf::MeshEntity*> current,
      std::vector<apf::MeshEntity*> next, apf::MeshTag* visited,
      apf::MeshTag* wtag, int peer)
  {
    assert(layers>=0);
    const int elmDim = m->getDimension();
    int yes=1;
    double* weight = new double[4];
    for(unsigned int i=0; i<4; i++)
      weight[i] = 0;

    // self owned ents are zero'th layer of ghosts
    for (unsigned int i=0;i<next.size();i++) {
      apf::MeshEntity* e = next[i];
      const int type = m->getType(e);
      assert( type == apf::Mesh::VERTEX || type == apf::Mesh::EDGE );
      weight[type] += parma::getEntWeight(m,e,wtag);
      m->setIntTag(e,visited,&yes);
      if(type == apf::Mesh::VERTEX ) {
        //ghost elements
        apf::Adjacent elms;
        m->getAdjacent(e, elmDim, elms);
        for(size_t k=0; k<elms.size(); k++) {
          if (!m->hasTag(elms[k],visited)) {
            m->setIntTag(elms[k],visited,&yes);
            weight[elmDim] += parma::getEntWeight(m,elms[k],wtag);
          }
        }
      }
    }

    for (int i=1;i<=layers;i++) {
      // when i==1: current contains shared vertices
      //            owned by the peer that are not ghosted
      for (unsigned int j=0;j<current.size();j++) {
        apf::MeshEntity* vertex = current[j];
        //ghost vertices and edges
        apf::Up edges;
        m->getUp(vertex,edges);
        for (int k=0;k<edges.n;k++) {
          apf::MeshEntity* edge = edges.e[k];
          apf::MeshEntity* v = getOtherVtx(m,edge,vertex);
          if (parma::isOwned(m, v) && !m->hasTag(v,visited)) {
            next.push_back(v);
            m->setIntTag(v,visited,&i);
            weight[0] += parma::getEntWeight(m,v,wtag);
          }
          if (parma::isOwned(m, edge) && !m->hasTag(edge,visited)) {
            weight[1] += parma::getEntWeight(m,edge,wtag);
            m->setIntTag(edge,visited,&i);
          }
        }
        //ghost elements
        apf::Adjacent elms;
        m->getAdjacent(vertex, elmDim, elms);
        for(size_t k=0; k<elms.size(); k++) {
          if (!m->hasTag(elms[k],visited)) {
            m->setIntTag(elms[k],visited,&i);
            weight[elmDim] += parma::getEntWeight(m,elms[k],wtag);
          }
        }
      }
      current=next;
      next.clear();
    }
    PCU_Debug_Print("ghostW peer %d vtx %f edge %f elm %f\n",
        peer, weight[0], weight[1], weight[3]);
    return weight;
  }

  double ownedWeight(apf::Mesh* m, apf::MeshTag* w, int dim) {
    apf::MeshIterator* it = m->begin(dim);
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
      GhostFinder(apf::Mesh* m, apf::MeshTag* w, int l)
        : mesh(m), wtag(w), layers(l) {
        depth = NULL;
      }
      /**
       * @brief get the weight of vertices ghosted to peer
       */
      double* weight(int peer) {
        int lvl = 0;
        depth = mesh->createIntTag("parma_depths_ver",1);
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
        mesh->end(itr);
        itr = mesh->begin(1);
        apf::MeshEntity* edge;
        while ((edge=mesh->iterate(itr))) {
          if (isSharedWithTarget(mesh,edge,peer)) {
            if (isOwned(mesh,edge))
              next.push_back(edge);
            else if (getOwner(mesh,edge)==peer)
              mesh->setIntTag(edge,depth,&lvl);
          }
        }
        mesh->end(itr);

        // current: peer owned vtx   next: self owned vtx and edges
        double* weight = runBFS(mesh,layers,current,next,depth,wtag,peer);
        for (unsigned int i=0;i<4;i++)
          apf::removeTagFromDimension(mesh,depth,i);
        mesh->destroyTag(depth);
        return weight;
      }
    private:
      GhostFinder();
      apf::Mesh* mesh;
      apf::MeshTag* wtag;
      int layers;
      apf::MeshTag* depth;
  };


  // can the associative container be modified to hold an array of weights per
  // peer????
  class GhostWeights : public Associative<double*> {
    public:
      GhostWeights(apf::Mesh* m, apf::MeshTag* w, Sides* s, int layers)
        : wtag(w), weight(0)
      {
        weight = new double[4];
        for(int dim=0; dim<4; dim++)
          weight[dim] += ownedWeight(m, wtag,dim);
        GhostFinder finder(m, wtag, layers);
        findGhosts(&finder, s);
        exchangeGhostsFrom();
        exchange();
      }
      ~GhostWeights() {
        const GhostWeights::Item* w;
        begin();
        while( (w = iterate()) )
          delete [] w->second;
        end();
        delete [] weight;
      }
      double self(int dim) {
        return weight[dim];
      }
    private:
      GhostWeights();
      apf::MeshTag* wtag;
      double* weight;
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
          PCU_Comm_Pack(ghost->first, ghost->second, 4*sizeof(double));
        end();
        PCU_Comm_Send();
        while (PCU_Comm_Listen()) {
          double* ghostsFromPeer = new double[4];
          PCU_Comm_Unpack(ghostsFromPeer, 4*sizeof(double));
          for(int i=0; i<4; i++)
            weight[i] += ghostsFromPeer[i];
        }
      }
      void exchange() {
        PCU_Comm_Begin();
        const GhostWeights::Item* ghost;
        begin();
        while( (ghost = iterate()) )
          PCU_Comm_Pack(ghost->first, weight, 4*sizeof(double));
        end();
        PCU_Comm_Send();
        while (PCU_Comm_Listen()) {
          double* peerWeight = new double[4];
          PCU_Comm_Unpack(peerWeight, 4*sizeof(double));
          int peer = PCU_Comm_Sender();
          set(peer, peerWeight);
        }
      }
  };

  class GhostToEntWeight : public Weights {
    public:
      GhostToEntWeight(GhostWeights* gw, int dim)
        : Weights(NULL,NULL,NULL) {
        wtag = gw->getTag();
        const GhostWeights::Item* ghost;
        gw->begin();
        while( (ghost = gw->iterate()) )
          set(ghost->first, ghost->second[dim]);
        gw->end();
        weight = gw->self(dim);
      }
      double self() {
        return weight;
      }
      apf::MeshTag* getTag() {
        return wtag;
      }
    private:
      GhostToEntWeight();
      double weight;
      apf::MeshTag* wtag;
  };

  Weights* convertGhostToEntWeight(GhostWeights* gw, int dim) {
    return new GhostToEntWeight(gw,dim);
  }

  GhostWeights* makeGhostWeights(apf::Mesh* m, apf::MeshTag* w, Sides* s,
      int layers) {
    return new GhostWeights(m, w, s, layers);
  }

  void destroyGhostWeights(GhostWeights* gw) {
    delete gw;
  }

} //end namespace
