#ifndef PARMA_SELECTOR_H
#define PARMA_SELECTOR_H
#include <apfMesh.h>
#include "parma_targets.h"

namespace parma {
  class Sides;
  class Selector {
    public:
      Selector(apf::Mesh* m, apf::MeshTag* w) : mesh(m), wtag(w) {}
      virtual ~Selector() {}
      virtual apf::Migration* run(Targets* tgts)=0;
    protected:
      apf::Mesh* mesh;
      apf::MeshTag* wtag;
    private:
      Selector();
  };
  Selector* makeVtxSelector(apf::Mesh* m, apf::MeshTag* w);
  Selector* makeEdgeSelector(apf::Mesh* m, apf::MeshTag* w);
  Selector* makeElmSelector(apf::Mesh* m, apf::MeshTag* w);
  Selector* makeElmLtVtxSelector(apf::Mesh* m, apf::MeshTag* w, double maxVtx);
  class Centroids;
  Selector* makeCentroidSelector(apf::Mesh* m, apf::MeshTag* w, Centroids* c);
  Selector* makeShapeSelector(apf::Mesh* m, apf::MeshTag* w, Centroids* c);
  Selector* makeWeldSelector(apf::Mesh* m, apf::MeshTag* w, Sides* s);
}
#endif
