#ifndef PARMA_DIJKSTRA_H_
#define PARMA_DIJKSTRA_H_

#include <apfMesh.h>
#include <set>
#include "parma_distQ.h"

namespace parma {
  class DijkstraContains {
    public:
      virtual ~DijkstraContains() {}
      virtual bool has(apf::MeshEntity* e)=0;
      virtual bool bdryHas(apf::MeshEntity* e)=0;
  };

  void dijkstra(apf::Mesh* m, DijkstraContains* c,
      apf::MeshEntity* src, apf::MeshTag* d);
  void dijkstra(apf::Mesh* m, DijkstraContains* c,
      parma::DistanceQueue<parma::Less>& pq, apf::MeshTag* d);
}

#endif
