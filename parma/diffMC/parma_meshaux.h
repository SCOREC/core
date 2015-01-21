#ifndef PARMA_MESHAUX_H_
#define PARMA_MESHAUX_H_

#include "apf.h"
#include "apfMesh.h"

inline void clearTag(apf::Mesh*& m, apf::MeshTag* t) {
  for(int d=0; d <= m->getDimension(); d++)
    apf::removeTagFromDimension(m, t, d);
}

inline apf::MeshEntity* getUpEnt(apf::Mesh* m, apf::MeshEntity* e) {
  return m->getUpward(e, 0);
}


#endif
