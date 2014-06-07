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
/* describes inter-part connectivity,
   see ph::encodeLinks */
  int* ilwork;
/* periodic masters array, one per node... */  
  int* iper;
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
/* nbc[i] is the index into essential boundary condition
   arrays of local node i (probably ;) */
  int* nbc;
/* ibc[i] is the essential boundary condition
   status code of essential BC node i */
  int* ibc;
/* bc[i][j] is the essential boundary condition
   value of 
   essential boundary condition j of
   essential BC node i */
  double** bc;
};

struct Output
{
  ~Output();
  Input* in;
  apf::Mesh* mesh;
/* again, int will overflow */
  int nGlobalEntities[4];
  int nOverlapNodes;
  int nOwnedNodes;
  int nGlobalNodes;
  int nBoundaryElements;
  int nMaxElementNodes;
  int nEssentialBCNodes;
  int nlwork; /* see ph::encodeLinks */
  AllBlocks blocks;
  EnsaArrays arrays;
};

void generateOutput(Input& in, apf::Mesh* mesh, Output& o);

}

#endif
