#ifndef PARMA_GHOSTOWNER_H 
#define PARMA_GHOSTOWNER_H
 
#include <apf.h>
#include <apfMesh.h>

namespace parma {
  int getOwner(apf::Mesh* m, apf::MeshEntity* v);
  bool isOwned(apf::Mesh* m, apf::MeshEntity* v);
}
#endif
