#ifndef PH_ADJACENT_H
#define PH_ADJACENT_H

#include <apfMesh.h>

namespace ph {

/* this is used to make sure the first vtx on side 0 is the first on side 1 */
void orderForPhastaInterface(int t, apf::MeshEntity** vin, apf::MeshEntity** vout, int offset = 0);

/* returns element vertices in the order expected by PHASTA */
void getVertices(apf::Mesh* m, apf::MeshEntity* e, apf::MeshEntity** v);

/* returns element vertices in the order expected by PHASTA,
   rotated such that the face (f) comes first */
void getBoundaryVertices(apf::Mesh* m, apf::MeshEntity* e, apf::MeshEntity* f,
    apf::MeshEntity** v);

extern int const* const face_apf2ph[apf::Mesh::TYPES];

}

#endif
