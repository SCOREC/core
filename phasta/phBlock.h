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

int getElementType(apf::Mesh* m, apf::MeshEntity* e);

enum {
  MAX_BLOCK_KEYS = 12
};

struct Blocks
{
  typedef std::map<BlockKey, int> Map;
/* indices starting from 1 ! ^ */
  Map keyToIndex;
  BlockKey keys[MAX_BLOCK_KEYS];
  int nElements[MAX_BLOCK_KEYS];
  int nElementNodes[MAX_BLOCK_KEYS];
  int getSize()
  {
    return keyToIndex.size();
  }
};

struct AllBlocks
{
  Blocks interior;
  Blocks boundary;
};

typedef std::set<apf::ModelEntity*> ModelBounds;
typedef std::vector<apf::MeshEntity*> MeshBounds;

void getAllBlocks(apf::Mesh* m, AllBlocks& b,
    ModelBounds& modelFaces);

std::string getBlockKeyPhrase(BlockKey& b, const char* prefix);

int getPhastaType(apf::Mesh* m, apf::MeshEntity* e);

}

#endif
