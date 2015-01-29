#ifndef PARMA_DIJKSTRA_H_
#define PARMA_DIJKSTRA_H_

#include <apfMesh.h>
#include <set>

namespace parma {
  class DijkstraContains {
    public:
      virtual ~DijkstraContains() {}
      virtual bool has(apf::MeshEntity* e)=0;
  };

  void dijkstra(apf::Mesh* m, apf::MeshTag* c, apf::MeshTag* d, 
      apf::MeshEntity* src);

  void dijkstra(apf::Mesh* m, DijkstraContains* c,
      apf::MeshEntity* src, apf::MeshTag* d);
}

#endif
