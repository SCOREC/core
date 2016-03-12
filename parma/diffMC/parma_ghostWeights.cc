#include <assert.h>
#include <PCU.h>
#include <apf.h>
#include <apfMesh.h>
#include <apfNumbering.h>
#include "parma_weights.h"
#include "parma_sides.h"
#include "parma_ghostOwner.h"

namespace {
  bool isOwnedByPeer(apf::Mesh* m,apf::MeshEntity* v, int peer) {
    if( ! m->isShared(v) ) return false;
    return (parma::getOwner(m,v) == peer);
  }

  // vertex based BFS
  // return an array of vertex, edge and element weights
  // - no one needs faces... yet
  // the returned array needs to be deallocated
  double* runBFS(apf::Mesh* m, int layers, std::vector<apf::MeshEntity*> current,
      apf::MeshTag* visited, apf::MeshTag* wtag, int peer)
  {
    assert(layers>=0);
    const int elmDim = m->getDimension();
    double* weight = new double[4];
    for(unsigned int i=0; i<4; i++)
      weight[i] = 0;

    std::vector<apf::MeshEntity*> next;
    for (int i=1;i<=layers;i++) {
      for (unsigned int j=0;j<current.size();j++) {
        apf::MeshEntity* vertex = current[j];
        apf::Adjacent elms;
        apf::Downward verts;
        apf::Downward edges;
        m->getAdjacent(vertex, elmDim, elms);
        for(size_t k=0; k<elms.size(); k++) {
          if (!m->hasTag(elms[k],visited)) {
            //ghost elements
            m->setIntTag(elms[k],visited,&i);
            weight[elmDim] += parma::getEntWeight(m,elms[k],wtag);
            //ghost edges
            const int nedges = m->getDownward(elms[k],1,edges);
            for(int l=0; l < nedges; l++)
              if (parma::isOwned(m, edges[l]) && !m->hasTag(edges[l],visited)) {
                m->setIntTag(edges[l],visited,&i);
                weight[1] += parma::getEntWeight(m,edges[l],wtag);
              }
            //ghost vertices
            const int nverts = m->getDownward(elms[k],0,verts);
            for(int l=0; l < nverts; l++)
              if (parma::isOwned(m, verts[l]) && !m->hasTag(verts[l],visited)) {
                next.push_back(verts[l]);
                m->setIntTag(verts[l],visited,&i);
                weight[0] += parma::getEntWeight(m,verts[l],wtag);
              }
          }
        }
      }
      current=next;
      next.clear();
    }
    PCU_Debug_Print("ghostW peer %d vtx %f edge %f elm %f\n",
        peer, weight[0], weight[1], weight[elmDim]);
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

      double* weight(int peer) {
        int lvl = 0;
        depth = mesh->createIntTag("parma_depths_ver",1);
        apf::MeshIterator* itr = mesh->begin(0);
        apf::MeshEntity* e;
        std::vector<apf::MeshEntity*> current;
        while( (e=mesh->iterate(itr)) )
          if( isOwnedByPeer(mesh,e,peer) )
            current.push_back(e);
        mesh->end(itr);
        //tag the un-owned boundary edges so their weights are not counted
        itr = mesh->begin(1);
        while( (e=mesh->iterate(itr)) )
          if( isOwnedByPeer(mesh,e,peer) )
            mesh->setIntTag(e,depth,&lvl);
        mesh->end(itr);

        // current: peer owned vtx
        double* weight = runBFS(mesh,layers,current,depth,wtag,peer);
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

  class GhostWeights : public Associative<double*> {
    public:
      GhostWeights(apf::Mesh* m, apf::MeshTag* w, Sides* s, int layers)
        : weight(0)
      {
        const int dim = m->getDimension();
        weight = new double[4];
        for(int d=0; d<=dim; d++)
          weight[d] = ownedWeight(m,w,d);
        for(int d=dim+1; d<=3; d++)
          weight[d] = 0;
        GhostFinder finder(m, w, layers);
        findGhosts(&finder, s);
        exchangeGhostsFrom();
        exchange();
        PCU_Debug_Print("totW vtx %f edge %f elm %f\n",
            weight[0], weight[1], weight[dim]);
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
    private:
      GhostToEntWeight();
      double weight;
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
