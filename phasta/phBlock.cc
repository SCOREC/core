#include "phBlock.h"
#include "phBC.h"
#include <apf.h>
#include <gmi.h>
#include <pcu_util.h>

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

bool BlockKeyInterface::operator<(BlockKeyInterface const& other) const
{
  if (elementType1 != other.elementType1)
    return elementType1 < other.elementType1;
  return this->BlockKey::operator<(other);
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

static void insertKeyInterface
(
  BlocksInterface&         b,
  BlockKeyInterface const& k
)
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
  if ((WEDGE == k.elementType) && (4 == k.nBoundaryFaceEdges))
    k.elementType = WEDGE_QUAD;
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

void getBoundaryBlocks(apf::Mesh* m, BCs& bcs, Blocks& b)
{
  int boundaryDim = m->getDimension() - 1;
  apf::MeshIterator* it = m->begin(boundaryDim);
  apf::MeshEntity* f;
  while ((f = m->iterate(it))) {
    apf::ModelEntity* me = m->toModel(f);
    if (m->getModelType(me) != boundaryDim)
      continue;
    if (getBCValue(m->getModel(), bcs.fields["DG interface"], (gmi_ent*) me) != 0) {
      apf::DgCopies dgCopies;
      m->getDgCopies(f, dgCopies);
      if (dgCopies.getSize() == 1) // This prevents adding interface elements...
        continue;
    }
    if (m->countUpward(f)>1)   // don't want interior region boundaries here...
      continue;
    apf::MeshEntity* e = m->getUpward(f, 0);
    BlockKey k;
    getBoundaryBlockKey(m, e, f, k);
    insertKey(b, k);
  }
  m->end(it);
}

static void applyTriQuadHackElement
(
  int elementType,
  int nBoundaryFaceEdges
)
{
  if ((WEDGE == elementType) && (4 == nBoundaryFaceEdges))
    elementType = WEDGE_QUAD;
  else if ((PYRAMID == elementType) && (3 == nBoundaryFaceEdges))
    elementType = PYRAMID_TRI;
}

void applyTriQuadHackInterface
(
  BlockKeyInterface& k
)
{
  applyTriQuadHackElement(k.elementType,  k.nBoundaryFaceEdges);
  applyTriQuadHackElement(k.elementType1, k.nBoundaryFaceEdges);
}

void getInterfaceBlockKey
(
  apf::Mesh*         m,
  apf::MeshEntity*   e0,
  apf::MeshEntity*   e1,
  apf::MeshEntity*   f,
  BlockKeyInterface& k
)
{
  k.elementType      = getPhastaType(m, e0);
  k.elementType1     = getPhastaType(m, e1);
  k.nElementVertices =
    apf::Mesh::adjacentCount[m->getType(e0)][0];
  k.nElementVertices1 =
    apf::Mesh::adjacentCount[m->getType(e1)][0];
  k.polynomialOrder = 1;
  k.nBoundaryFaceEdges =
    apf::Mesh::adjacentCount[m->getType(f)][1];
  applyTriQuadHackInterface(k);
}

void getInterfaceBlocks(apf::Mesh* m, BCs& bcs, BlocksInterface& b)
{
  int interfaceDim = m->getDimension() - 1;
  apf::MeshIterator* it = m->begin(interfaceDim);
  apf::MeshEntity* face;

  while ((face = m->iterate(it))) {
    apf::ModelEntity* me = m->toModel(face);
    if (getBCValue(m->getModel(), bcs.fields["DG interface"], (gmi_ent*) me) == 0)
      continue;
    if (m->getModelType(me) != interfaceDim)
      continue;
    apf::DgCopies dgCopies;
    m->getDgCopies(face, dgCopies);
    if (dgCopies.getSize() != 1)
      continue;
    apf::MeshEntity* e0 = m->getUpward(face, 0);
    PCU_ALWAYS_ASSERT(dgCopies[0].peer == m->getPCU()->Self());
    apf::MeshEntity* e1 = m->getUpward(dgCopies[0].entity, 0);
    /* in order to avoid repetition of elements */
    if (e0 > e1)
      continue;

    BlockKeyInterface k;
    getInterfaceBlockKey(m, e0, e1, face, k);
    insertKeyInterface(b, k);
  }
  m->end(it);
}

void getAllBlocks(apf::Mesh* m, BCs& bcs, AllBlocks& b)
{
  getInteriorBlocks(m, b.interior);
  getBoundaryBlocks(m, bcs, b.boundary);
  getInterfaceBlocks(m, bcs, b.interface);
}

std::string getPolyOrder
(
  int polyOrder
)
{
  static const char* const polyTable[5] =
  {NULL
  ,"linear "
  ,"quadratic "
  ,"cubic "
  ,"quartic "};
  return polyTable[polyOrder];
}

std::string getElementType
(
  int elementType
)
{
  static const char* typeTable[TYPES] =
  {NULL
  ,"tetrahedron "
  ,"hexahedron "
  ,"wedge triface" // this will push this into the interior element as well but have to propagate
  ,"wedge quadface "
  ,"pyramid quadface"
  ,"pyramid triface "};
  return typeTable[elementType];
}

std::string getBlockKeyPhrase
(
  BlockKey& b,
  const char* prefix
)
{
  std::string s = prefix;
  s += getPolyOrder(b.polynomialOrder);
  s += getElementType(b.elementType);
  return s;
}

std::string getBlockKeyPhraseInterface
(
  BlockKeyInterface& b,
  const char* prefix
)
{
  std::string s = prefix;
  s += getPolyOrder(b.polynomialOrder);
  s += getElementType(b.elementType);
  s += getElementType(b.elementType1);
  return s;
}

}
