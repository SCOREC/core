#ifndef PH_ADJACENT_H
#define PH_ADJACENT_H

#include <apfMesh.h>

namespace ph {

/* returns element vertices in the order expected by PHASTA */
void getVertices(apf::Mesh* m, apf::MeshEntity* e, apf::MeshEntity** v);

/* returns element vertices in the order expected by PHASTA,
   rotated such that the face (f) comes first */
void getBoundaryVertices(apf::Mesh* m, apf::MeshEntity* e, apf::MeshEntity* f,
    apf::MeshEntity** v);

extern int const* const face_apf2ph[apf::Mesh::TYPES];

}

#endif
