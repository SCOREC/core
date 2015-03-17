#ifndef PARMA_GRAPHDIST_H_
#define PARMA_GRAPHDIST_H_

#include <apfMesh.h>

namespace parma {
  apf::MeshTag* measureGraphDist(apf::Mesh* m);
  apf::MeshTag* getDistTag(apf::Mesh* m);
  void clearDistTag(apf::Mesh* m, apf::MeshTag* dist, apf::Migration* plan);
}

#endif
