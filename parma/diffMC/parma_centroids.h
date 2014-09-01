#ifndef PARMA_CENTROIDS_H
#define PARMA_CENTROIDS_H
#include <apfMesh.h>
#include "parma_associative.h"

namespace parma {
  class Sides;
  class Centroids : public Associative<apf::Vector3> {
    public:
      Centroids(apf::Mesh* m, apf::MeshTag* w, Sides* s);
      ~Centroids() {}
      apf::Vector3 self();
    private:
      apf::Vector3 centroid;
      double weight;
      apf::Vector3 selfCentroid(apf::Mesh* m, apf::MeshTag* w);
      void init(apf::Mesh* m, Sides* s);
      Centroids();
  };
}

#endif

