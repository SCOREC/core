#include "splitMeshOnGFace.h"

extern void M_splitMeshOnGFace(pUnstructuredMesh, pGFace, int keepOrig, pPList nosplit);

void ph::splitMeshOnGFace(pUnstructuredMesh pmesh, pGFace gf) {
  M_splitMeshOnGFace(pmesh, gf, 0, 0);
}
