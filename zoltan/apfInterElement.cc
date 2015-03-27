#include <PCU.h>
#include "apfZoltan.h"
#include <apfMesh.h>

namespace apf {

static bool hasOtherSide(Mesh* m, MeshEntity* s)
{
  if (m->isShared(s))
    return true;
  if (!m->hasMatching())
    return false;
  Matches matches;
  m->getMatches(s, matches);
  return matches.getSize();
}

static Copy getOtherSide(Mesh* m, MeshEntity* s)
{
  if (m->isShared(s))
    return getOtherCopy(m, s);
  Matches matches;
  m->getMatches(s, matches);
  assert(matches.getSize() == 1);
  return matches[0];
}

static MeshEntity* getSideElement(Mesh* m, MeshEntity* s)
{
  return m->getUpward(s, 0);
}

static long getElementGid(GlobalNumbering* gn, MeshEntity* e)
{
  return getNumber(gn, Node(e, 0));
}

static void packOtherGid(GlobalNumbering* gn, MeshEntity* s)
{
  Mesh* m = getMesh(gn);
  Copy other = getOtherSide(m, s);
  long gid = getElementGid(gn, getSideElement(m, s));
  PCU_COMM_PACK(other.peer, other.entity);
  PCU_COMM_PACK(other.peer, gid);
}

static void unpackOtherGid(Mesh* m, MeshTag* t)
{
  MeshEntity* s;
  PCU_COMM_UNPACK(s);
  long gid;
  PCU_COMM_UNPACK(gid);
  m->setLongTag(s, t, &gid);
}

MeshTag* tagOpposites(GlobalNumbering* gn, const char* name)
{
  Mesh* m = getMesh(gn);
  MeshTag* t = m->createLongTag(name, 1);
  int sideDim = m->getDimension() - 1;
  PCU_Comm_Begin();
  MeshIterator* it = m->begin(sideDim);
  MeshEntity* e;
  while ((e = m->iterate(it)))
    if (hasOtherSide(m, e))
      packOtherGid(gn, e);
  m->end(it);
  PCU_Comm_Send();
  while (PCU_Comm_Receive())
    unpackOtherGid(m, t);
  return t;
}

static void packOtherLid(Numbering* n, MeshEntity* s, int i)
{
  Mesh* m = getMesh(n);
  Copy other = getOtherSide(m, s);
  PCU_COMM_PACK(other.peer, other.entity);
  PCU_COMM_PACK(other.peer, i);
}

static void unpackOtherLid(Mesh* m, Numbering* n, LocalCopy* e2e, int nsides)
{
  MeshEntity* s;
  PCU_COMM_UNPACK(s);
  int lid;
  PCU_COMM_UNPACK(lid);
  MeshEntity* e = getSideElement(m, s);
  int i = getNumber(n, e, 0, 0);
  int which; bool flip; int rotate;
  getAlignment(m, e, s, which, flip, rotate);
  e2e[i * nsides + which] = LocalCopy(PCU_Comm_Sender(), lid);
}

LocalCopy* getLocalElementToElement(apf::Mesh* m)
{
  int dim = m->getDimension();
  int sideDim = dim - 1;
  int type = getFirstType(m, dim);
  int nsides = apf::Mesh::adjacentCount[type][dim - 1];
  LocalCopy* e2e = new LocalCopy[m->count(dim) * nsides];
  Numbering* n = numberElements(m, "getLocalElementToElement");
  PCU_Comm_Begin();
  MeshIterator* it = m->begin(dim);
  MeshEntity* e;
  int i = 0;
  while ((e = m->iterate(it))) {
    Downward sides;
    m->getDownward(e, sideDim, sides);
    for (int j = 0; j < nsides; ++j) {
      Up up;
      m->getUp(sides[j], up);
      if (up.n == 1) {
        if (hasOtherSide(m, sides[j]))
          packOtherLid(n, sides[j], i);
        else
          e2e[i * nsides + j] = LocalCopy();
      } else {
        MeshEntity* oe = up.e[1 - findIn(up.e, 2, e)];
        int which; bool flip; int rotate;
        getAlignment(m, e, sides[j], which, flip, rotate);
        int oi = getNumber(n, oe, 0, 0);
        e2e[i * nsides + j] = LocalCopy(PCU_Comm_Self(), oi);
      }
    }
    ++i;
  }
  m->end(it);
  PCU_Comm_Send();
  while (PCU_Comm_Receive())
    unpackOtherLid(m, n, e2e, nsides);
  return e2e;
}

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
