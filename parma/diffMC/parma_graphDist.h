#ifndef PARMA_GRAPHDIST_H_
#define PARMA_GRAPHDIST_H_

#include <apfMesh.h>
#include "parma_distQ.h"

namespace parma {
  apf::MeshTag* measureGraphDist(apf::Mesh* m);
  DistanceQueue<Greater> * BoundaryVertices(apf::Mesh* m, apf::MeshTag* d);
}

#endif
