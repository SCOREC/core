#ifndef PARMA_DIJKSTRA_H_
#define PARMA_DIJKSTRA_H_

#include <apfMesh.h>

namespace parma {
  void dijkstra(apf::Mesh* m, apf::MeshTag* c, apf::MeshTag* d, 
      apf::MeshEntity* src);
}

#endif

