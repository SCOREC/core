#ifndef PARMA_BDRYVTX_H_
#define PARMA_BDRYVTX_H_

#include <apfMesh.h>

namespace parma {
  class BdryVtxItr {
    public:
      virtual apf::MeshEntity* next() = 0;
      virtual ~BdryVtxItr() {}
  };

  BdryVtxItr* makeBdryVtxDistItr(apf::Mesh* m, apf::MeshTag* d);
}

#endif
