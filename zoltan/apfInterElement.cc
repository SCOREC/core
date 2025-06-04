#include "apfZoltan.h"
#include <apfMesh.h>
#include <pcu_util.h>

namespace apf {

static long getOpposite(GlobalNumbering* gn, MeshTag* opposites,
    MeshEntity* element, int side_i)
{
  Downward sides;
  Mesh* m = getMesh(gn);
  int dim = getDimension(m, element);
  m->getDownward(element, dim - 1, sides);
  MeshEntity* side = sides[side_i];
  if (m->hasTag(side, opposites)) {
    long gid;
    m->getLongTag(side, opposites, &gid);
    return gid;
  }
  Up elements;
  m->getUp(side, elements);
  if (elements.n == 1)
    return -1;
  int i = findIn(elements.e, 2, element);
  return getNumber(gn, Node(elements.e[1 - i], 0));
}

int* getElementToElement(apf::Mesh* m)
{
  GlobalNumbering* gn = makeGlobal(numberElements(m, "zb_element"));
  MeshTag* opposites = tagOpposites(gn, "zb_opposite");
  int dim = m->getDimension();
  int type = getFirstType(m, dim);
  int nsides = apf::Mesh::adjacentCount[type][dim - 1];
  int* e2e = new int[m->count(dim) * nsides];
  MeshIterator* it = m->begin(dim);
  int i = 0;
  MeshEntity* e;
  while ((e = m->iterate(it)))
    for (int j = 0; j < nsides; ++j)
      e2e[i++] = getOpposite(gn, opposites, e, j);
  m->end(it);
  m->destroyTag(opposites);
  destroyGlobalNumbering(gn);
  return e2e;
}

}
