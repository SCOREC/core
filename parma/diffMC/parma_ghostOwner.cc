#include <apf.h>
#include <apfMesh.h>
#include "parma_ghostOwner.h"

namespace parma {
  int getOwner(apf::Mesh* m, apf::MeshEntity* v) {
    apf::Parts res;
    m->getResidence(v, res);
    return *(res.begin());
  }

  bool isOwned(apf::Mesh* m, apf::MeshEntity* v) {
    return m->getPCU()->Self() == getOwner(m,v);
  }
}
