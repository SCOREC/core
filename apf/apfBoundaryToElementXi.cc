#include "apf.h"
#include "apfMesh.h"
#include "apfShape.h"
#include <pcu_util.h>

namespace apf {

/* at some point we may want to work these tables
 * into a more central place along with the
 * shape functions...
 */
static Vector3 const edge_vert_xi[2] = {
  Vector3(-1,0,0),
  Vector3(1,0,0),
};
static Vector3 const tri_vert_xi[3] = {
  Vector3(0,0,0),
  Vector3(1,0,0),
  Vector3(0,1,0),
};
static Vector3 const tet_vert_xi[4] = {
  Vector3(0,0,0),
  Vector3(1,0,0),
  Vector3(0,1,0),
  Vector3(0,0,1),
};
static Vector3 const* const elem_vert_xi[Mesh::TYPES] = {
  0, /* vertex */
  edge_vert_xi,
  tri_vert_xi,
  0, /* quad */
  tet_vert_xi,
  0, /* hex */
  0, /* prism */
  0  /* pyramid */
};

Vector3 boundaryToElementXi(
    Mesh* m,
    MeshEntity* boundary,
    MeshEntity* element,
    Vector3 const& xi)
{
  Downward bv;
  int nbv = m->getDownward(boundary, 0, bv);
  int bt = m->getType(boundary);
  Downward ev;
  int nev = m->getDownward(element, 0, ev);
  int et = m->getType(element);
  /* the boundary parametric space
   * and parent parametric space
   * are linear functions of one another
   */
  EntityShape* shape = apf::getLagrange(1)->getEntityShape(bt);
  NewArray<double> shape_vals;
  shape->getValues(m, boundary, xi, shape_vals);
  Vector3 exi(0,0,0);
  for (int i = 0; i < nbv; ++i) {
    int evi = findIn(ev, nev, bv[i]);
    PCU_ALWAYS_ASSERT(evi >= 0);
    exi += elem_vert_xi[et][evi] * shape_vals[i];
  }
  return exi;
}

}
