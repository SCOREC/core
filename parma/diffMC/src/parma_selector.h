#ifndef PARMA_SELECTOR_H
#define PARMA_SELECTOR_H
#include <apfMesh.h>
#include "parma_targets.h"

namespace parma {
  class Selector {
    public:
      Selector(apf::Mesh* m, apf::MeshTag* w);
      virtual apf::Migration* run(Targets* tgts)=0;
    private:
      apf::Mesh* mesh;
      apf::MeshTag* vtag;
      apf::MeshTag* wtag;
      Selector();
  };
  Selector* makeVtxSelector(apf::Mesh* m, apf::MeshTag* w);
}
#endif
