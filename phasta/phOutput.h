#ifndef PH_OUTPUT_H
#define PH_OUTPUT_H

#include "phInput.h"
#include "phBlock.h"

namespace apf {
class Mesh;
}

namespace ph {

struct EnsaArrays
{
  double* coordinates;
  int* ilwork; /* ??? */
/* note: int will overflow at about 2 billion total nodes */
  int* globalNodeNumbers;
/* ien[i][j][k] is the local vertex id of
   vertex k of
   element j of
   interior block i */
  int*** ien;
/* ienb[i][j][k] is the local vertex id of
   vertex k of
   element j of
   boundary block i */
  int*** ienb;
/* ibcb[i][j][k] is the natural boundary condition
   status code
   number k in [0,1] of
   element j of
   boundary block i */
  int*** ibcb;
/* bcb[i][j][k] is the natural boundary condition
   value for
   boundary condition k of
   element j of
   boundary block i */
  double*** bcb;
};

struct Output
{
  Input* in;
  apf::Mesh* mesh;
  int nGlobalEntities[4];
  int nOverlapNodes;
  int nOwnedNodes;
  int nGlobalNodes;
  int nBoundaryElements;
  int nMaxElementNodes;
  int nEssentialBCNodes;
  int nlwork; /* ??? */
  AllBlocks blocks;
  EnsaArrays arrays;
};

}

#endif
