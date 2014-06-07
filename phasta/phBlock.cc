#include "phBlock.h"
#include <apf.h>

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
  int idx = b.keyToIndex[k];
/* rely on the fact that previously unmapped keys
   get mapped to zero by default to identify them,
   this is part of the reason why the indices
   are 1-based */
  if (idx) {
    --idx;
    ++(b.nElements[idx]);
  } else {
    idx = b.keyToIndex.size() - 1;
    b.keyToIndex[k] = idx + 1;
    b.nElements[idx] = 1;
    b.keys[idx] = k;
    b.nElementNodes[idx] = k.nElementVertices;
  }
}

void getInteriorBlocks(apf::Mesh* m, Blocks& b)
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

void getBoundaryBlocks(apf::Mesh* m, Blocks& b,
    ModelBounds& modelFaces)
{
  APF_ITERATE(ModelBounds, modelFaces, mit) {
    apf::ModelEntity* modelFace = *mit;
    apf::MeshIterator* it = m->begin(m->getDimension() - 1);
    apf::MeshEntity* f;
    while ((f = m->iterate(it))) {
      if (m->toModel(f) != modelFace)
        continue;
      apf::MeshEntity* e = m->getUpward(f, 0);
      BlockKey k;
      k.nElementVertices =
        apf::Mesh::adjacentCount[m->getType(e)][0];
      k.polynomialOrder = 1;
      k.nBoundaryFaceEdges =
        apf::Mesh::adjacentCount[m->getType(f)][1];
      k.elementType = getPhastaType(m, e);
      applyTriQuadHack(k);
      insertKey(b, k);
    }
    m->end(it);
  }
}

void getAllBlocks(apf::Mesh* m, AllBlocks& b,
    ModelBounds& modelFaces)
{
  getInteriorBlocks(m, b.interior);
  getBoundaryBlocks(m, b.boundary, modelFaces);
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
