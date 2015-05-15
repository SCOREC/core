#ifndef PH_BLOCK_H
#define PH_BLOCK_H

#include <apfMesh.h>
#include <map>
#include <set>
#include <vector>

namespace ph {

enum {
  TETRAHEDRON  = 1,
  HEXAHEDRON   = 2,
  WEDGE        = 3,
  WEDGE_TRI    = 3,
  WEDGE_QUAD   = 4,
  PYRAMID      = 5,
  PYRAMID_QUAD = 5,
  PYRAMID_TRI  = 6,
  TYPES        = 7
};

struct BlockKey
{
  int nElementVertices;
  int polynomialOrder;
  int nBoundaryFaceEdges;
  int elementType;
  bool operator<(BlockKey const& other) const;
};

struct BlockKeyInterface : public BlockKey
{
  int nElementVertices1;
  int elementType1;
  bool operator<(BlockKeyInterface const& other) const;
};

enum {
  MAX_BLOCK_KEYS = 12
};

struct BlocksCommon
{
  typedef std::map<BlockKey, int> Map;
  Map keyToIndex;
  int nElements[MAX_BLOCK_KEYS];
  int nElementNodes[MAX_BLOCK_KEYS];
  int getSize()
  {
    return keyToIndex.size();
  }
};

struct Blocks : public BlocksCommon
{
  BlockKey keys[MAX_BLOCK_KEYS];
};

struct BlocksInterface : public BlocksCommon
{
  BlockKeyInterface keys[MAX_BLOCK_KEYS];
};

struct AllBlocks
{
  Blocks interior;
  Blocks boundary;
  BlocksInterface interface;
};

void getAllBlocks(apf::Mesh* m, AllBlocks& b);

std::string getBlockKeyPhrase(BlockKey& b, const char* prefix);

void getInteriorBlockKey(apf::Mesh* m, apf::MeshEntity* e, BlockKey& k);
void getBoundaryBlockKey(apf::Mesh* m, apf::MeshEntity* e,
    apf::MeshEntity* f, BlockKey& k);
void getInterfaceBlockKey
(
  apf::Mesh*         m, 
  apf::MeshEntity*   e0, 
  apf::MeshEntity*   e1, 
  apf::MeshEntity*   f, 
  BlockKeyInterface& k
);

}

#endif
