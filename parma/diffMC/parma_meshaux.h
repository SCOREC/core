#ifndef PARMA_MESHAUX_H_
#define PARMA_MESHAUX_H_

#include "apf.h"
#include "apfMesh.h"

typedef apf::DynamicArray<apf::MeshEntity*> dynEntArr;

inline void clearTag(apf::Mesh*& m, apf::MeshTag* t) {
  for(int d=0; d <= m->getDimension(); d++)
    apf::removeTagFromDimension(m, t, d);
}

inline apf::MeshEntity* getUpEnt(apf::Mesh* m, apf::MeshEntity* e) {
  return m->getUpward(e, 0);
}

inline void getDwn2ndAdj(apf::Mesh* m, apf::MeshEntity* ent, dynEntArr& adj) {
  const int dim = apf::getDimension(m, ent);
  apf::getBridgeAdjacent(m, ent, dim - 1, dim, adj);
}

inline int getNumFaceOnPb(apf::Mesh* m, const int destPid, apf::Downward& adjF,
    const int numFaces) {
  int numPbFaces = 0;

  apf::Parts pbPid;
  pbPid.insert(m->getId());
  pbPid.insert(destPid);

  for(int i=0; i<numFaces; i++) {
    assert(m->getDimension()-1 == apf::getDimension(m, adjF[i]));
    apf::Parts resPid;
    m->getResidence(adjF[i], resPid);
    if ( resPid == pbPid )
      numPbFaces++;
  }
  return numPbFaces;
}

#endif
