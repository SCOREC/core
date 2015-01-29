#ifndef PARMA_MESHAUX_H_
#define PARMA_MESHAUX_H_

#include "apf.h"
#include "apfMesh.h"

inline void clearTag(apf::Mesh*&, apf::MeshTag*);
inline apf::MeshEntity* getUpEnt(apf::Mesh*, apf::MeshEntity*);
inline void getEdgeAdjVtx(apf::Mesh*, apf::MeshEntity*, apf::Adjacent&);
inline void getDwn2ndAdj(apf::Mesh*, apf::MeshEntity*, apf::Adjacent&);
inline bool onBoundary(apf::Mesh*, apf::MeshEntity*);


void clearTag(apf::Mesh*& m, apf::MeshTag* t) {
  for(int d=0; d <= m->getDimension(); d++)
    apf::removeTagFromDimension(m, t, d);
}

apf::MeshEntity* getUpEnt(apf::Mesh* m, apf::MeshEntity* e) {
  return m->getUpward(e, 0);
}

void getEdgeAdjVtx(apf::Mesh* m, apf::MeshEntity* v, apf::Adjacent& adj) {
  int bridge = 1; int tgt = 0;
  getBridgeAdjacent(m, v, bridge, tgt, adj);
}

void getDwn2ndAdj(apf::Mesh* m, apf::MeshEntity* ent, apf::Adjacent& adj) {
  const int dim = apf::getDimension(m, ent);
  apf::getBridgeAdjacent(m, ent, dim - 1, dim, adj);
}

bool onBoundary(apf::Mesh* m, apf::MeshEntity* e) {
  int gd = m->getModelType(m->toModel(e));
  int md = m->getDimension();
  bool shared = m->isShared(e);
  return shared || (!shared && gd < md);
}

#endif
