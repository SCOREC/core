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

inline void getEdgeAdjVtx(apf::Mesh* m, apf::MeshEntity* vtx,
    dynEntArr& adjVtx) {
  apf::getBridgeAdjacent(m, vtx, 1, 0, adjVtx);
}

inline void getDwn2ndAdj(apf::Mesh* m, apf::MeshEntity* ent, dynEntArr& adj) {
  const int dim = apf::getDimension(m, ent);
  apf::getBridgeAdjacent(m, ent, dim - 1, dim, adj);
}

inline int getNumTaggedEnts(apf::Mesh* m, apf::MeshTag* vtag, int visited,
    apf::Downward& ents, const int numEnts) {
  int taggedEnts = 0;
  int tval = -1;
  for(int i=0; i<numEnts; i++) {
    if( m->hasTag(ents[i], vtag) ) {
      m->getIntTag(ents[i], vtag, &tval);
      if ( tval == visited )
        taggedEnts++;
    }
  }
  return taggedEnts;
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
