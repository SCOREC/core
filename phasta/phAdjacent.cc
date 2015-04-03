#include "phAdjacent.h"
#include <maMesh.h> /* rotation functions ! */

namespace ph {

void orderForPhasta(int t, apf::MeshEntity** vin, apf::MeshEntity** vout)
{
  /* ph2apf[type][i] is the apf index of phasta index (i) for that type */
  /* in all of these, we also try to preserve the "first" faces */
  /* phasta's tet orders the curl of the bottom face outward */
  /* the first face is the bottom face */
  static int const tet_ph2apf[4] = {0,2,1,3};
  /* despite what some comments in phParAdapt say, the pyramid
     quad face curls inward by their output ordering, making
     this the identity map */
  static int const pyr_ph2apf[5] = {0,1,2,3,4};
  /* prism ordering is identical, cool. */
  static int const pri_ph2apf[6] = {0,1,2,3,4,5};
  static int const* const ph2apf[apf::Mesh::TYPES] =
  {0 //vertex
  ,0 //edge
  ,0 //triangle
  ,0 //quad
  ,tet_ph2apf
  ,0 //hex
  ,pri_ph2apf
  ,pyr_ph2apf};
  assert(ph2apf[t]);
  int nv = apf::Mesh::adjacentCount[t][0];
  for (int i = 0; i < nv; ++i)
    vout[i] = vin[ph2apf[t][i]];
}

/* phasta defines its faces with the same local vertex indices as
   we define our faces, but their local vertex indices differ
   as seen above. this results in the following table of
   element face positions from apf to phasta */

static int const tet_face_apf2ph[4] = {0,3,2,1};
int const* const face_apf2ph[apf::Mesh::TYPES] = {
  0,   //vertex
  0,   //edge
  0,   //triangle
  0,   //quad
  tet_face_apf2ph,
  0,   //hex
  0,   //prism
  0,   //pyramid
};

void getVertices(apf::Mesh* m, apf::MeshEntity* e, apf::MeshEntity** v)
{
  apf::Downward v0;
  m->getDownward(e, 0, v0);
  orderForPhasta(m->getType(e), v0, v);
}

void getBoundaryVertices(apf::Mesh* m, apf::MeshEntity* e, apf::MeshEntity* f,
    apf::MeshEntity** v)
{
  /* for rotation codes, see maTables.cc */
  /* the rotation code for a tet that brings
     apf face (i) to the bottom is tet_rot[i] */
  static int const tet_rot[4] = {0,2,4,1};
  /* pyramid (quad face) rotation codes given
     that the face of interest is apf index (i) */
  static int const pyr_rot[5] = {0,0,1,2,3};
  /* prism rotation codes given
     that the face of interest is apf index (i) */
  /* by the forces of nature, they match the pyramid codes */
  static int const pri_rot[5] = {0,0,1,2,3};
  static int const* const rot[apf::Mesh::TYPES] =
  {0 //vertex
  ,0 //edge
  ,0 //triangle
  ,0 //quad
  ,tet_rot
  ,0 //hex
  ,pri_rot
  ,pyr_rot};
  apf::Downward faces;
  int nf = m->getDownward(e, 2, faces);
  int i = apf::findIn(faces, nf, f);
  int t = m->getType(e);
  apf::Downward v0;
  ma::rotateEntity(m, e, rot[t][i], v0);
  orderForPhasta(t, v0, v);
}

}
