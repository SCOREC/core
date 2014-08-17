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

int getPhastaType(apf::Mesh* m, apf::MeshEntity* e)
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

int getBoundaryFaceEdges(apf::Mesh* m, apf::MeshEntity* e)
{
  static int const table[apf::Mesh::TYPES] = 
  {-1  //vertex
  ,-1  //edge
  ,-1  //triangle
  ,-1  //quad
  ,3 //tet
  ,4 //hex
  ,3 //wedge
  ,3}; //pyramid
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

static void getInteriorBlocks(apf::Mesh* m, Blocks& b)
{
  apf::MeshIterator* it = m->begin(m->getDimension());
  apf::MeshEntity* e;
  while ((e = m->iterate(it))) {
    BlockKey k;
    k.nElementVertices =
      apf::Mesh::adjacentCount[m->getType(e)][0];
    k.polynomialOrder = 1;
    k.nBoundaryFaceEdges = getBoundaryFaceEdges(m, e);
    k.elementType = getPhastaType(m, e);
    insertKey(b, k);
  }
  m->end(it);
}

void applyTriQuadHack(BlockKey& k)
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
  k.nElementVertices =
    apf::Mesh::adjacentCount[m->getType(e)][0];
  k.polynomialOrder = 1;
  k.nBoundaryFaceEdges =
    apf::Mesh::adjacentCount[m->getType(f)][1];
  k.elementType = getPhastaType(m, e);
  applyTriQuadHack(k);
}

void getBoundaryBlocks(apf::Mesh* m, Blocks& b)
{
  gmi_model* gm = m->getModel();
  gmi_iter* git = gmi_begin(gm, m->getDimension() - 1);
  gmi_ent* gf;
  while ((gf = gmi_next(gm, git))) {
    apf::ModelEntity* modelFace = (apf::ModelEntity*)gf;
    apf::MeshIterator* it = m->begin(m->getDimension() - 1);
    apf::MeshEntity* f;
    while ((f = m->iterate(it))) {
      if (m->toModel(f) != modelFace)
        continue;
      apf::MeshEntity* e = m->getUpward(f, 0);
      BlockKey k;
      getBoundaryBlockKey(m, e, f, k);
      insertKey(b, k);
    }
    m->end(it);
  }
  gmi_end(gm, git);
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
