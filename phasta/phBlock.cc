#include "phBlock.h"
#include <apf.h>
#include <gmi.h>

namespace ph {

bool BlockKey::operator<(BlockKey const& other) const
{
  if (nElementVertices != other.nElementVertices)
    return nElementVertices < other.nElementVertices;
  if (elementType != other.elementType)
    return elementType < other.elementType;
  if (nBoundaryFaceEdges != other.nBoundaryFaceEdges)
    return nBoundaryFaceEdges < other.nBoundaryFaceEdges;
  return polynomialOrder < other.polynomialOrder;
}

static int getPhastaType(apf::Mesh* m, apf::MeshEntity* e)
{
  static int const table[apf::Mesh::TYPES] = 
  {-1  //vertex
  ,-1  //edge
  ,-1  //triangle
  ,-1  //quad
  ,TETRAHEDRON
  ,HEXAHEDRON
  ,WEDGE
  ,PYRAMID};
  return table[m->getType(e)];
}

static void insertKey(Blocks& b, BlockKey const& k)
{
  if (b.keyToIndex.count(k)) {
    int idx = b.keyToIndex[k];
    ++(b.nElements[idx]);
  } else {
    int idx = b.keyToIndex.size();
    b.keyToIndex[k] = idx;
    b.nElements[idx] = 1;
    b.keys[idx] = k;
    b.nElementNodes[idx] = k.nElementVertices;
  }
}

static void getBlockKeyCommon(apf::Mesh* m, apf::MeshEntity* e, BlockKey& k)
{
  k.elementType = getPhastaType(m, e);
  k.nElementVertices =
    apf::Mesh::adjacentCount[m->getType(e)][0];
  k.polynomialOrder = 1;
}

void getInteriorBlockKey(apf::Mesh* m, apf::MeshEntity* e, BlockKey& k)
{
  getBlockKeyCommon(m, e, k);
  /* what this value is should not matter much for interior elements */
  k.nBoundaryFaceEdges = k.elementType == HEXAHEDRON ? 4 : 3;
}

static void getInteriorBlocks(apf::Mesh* m, Blocks& b)
{
  apf::MeshIterator* it = m->begin(m->getDimension());
  apf::MeshEntity* e;
  while ((e = m->iterate(it))) {
    BlockKey k;
    getInteriorBlockKey(m, e, k);
    insertKey(b, k);
  }
  m->end(it);
}

static void applyTriQuadHack(BlockKey& k)
{
  /* distinguish between WEDGE_TRI (wedge with triangle on boundary)
     and WEDGE_QUAD (wedge with quad on boundary) */
  if (WEDGE == k.elementType)
    /* WEDGE_TRI and WEDGE_QUAD are lined up just right for this */
    k.elementType = k.nBoundaryFaceEdges;
  /* same hack for pyramids */
  else if ((PYRAMID == k.elementType) && (3 == k.nBoundaryFaceEdges))
    k.elementType = PYRAMID_TRI;
}

void getBoundaryBlockKey(apf::Mesh* m, apf::MeshEntity* e,
    apf::MeshEntity* f, BlockKey& k)
{
  getBlockKeyCommon(m, e, k);
  k.nBoundaryFaceEdges =
    apf::Mesh::adjacentCount[m->getType(f)][1];
  applyTriQuadHack(k);
}

void getBoundaryBlocks(apf::Mesh* m, Blocks& b)
{
  int boundaryDim = m->getDimension() - 1;
  apf::MeshIterator* it = m->begin(boundaryDim);
  apf::MeshEntity* f;
  while ((f = m->iterate(it))) {
    apf::ModelEntity* me = m->toModel(f);
    if (m->getModelType(me) != boundaryDim)
      continue;
    if (m->countUpward(f)>1)   // don't want interior region boundaries here...
      continue;
    apf::MeshEntity* e = m->getUpward(f, 0);
    BlockKey k;
    getBoundaryBlockKey(m, e, f, k);
    insertKey(b, k);
  }
  m->end(it);
}

void getAllBlocks(apf::Mesh* m, AllBlocks& b)
{
  getInteriorBlocks(m, b.interior);
  getBoundaryBlocks(m, b.boundary);
}

std::string getBlockKeyPhrase(BlockKey& b, const char* prefix)
{
  std::string s = prefix;
  static const char* const polyTable[5] =
  {NULL
  ,"linear "
  ,"quadratic "
  ,"cubic "
  ,"quartic "};
  s += polyTable[b.polynomialOrder];
  static const char* typeTable[TYPES] =
  {NULL
  ,"tetrahedron "
  ,"hexahedron "
  ,"wedge "
  ,"wedge quadface "
  ,"pyramid "
  ,"pyramid triface "};
  s += typeTable[b.elementType];
  return s;
}

}
