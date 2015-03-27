#ifndef PH_OUTPUT_H
#define PH_OUTPUT_H

#include "phInput.h"
#include "phBlock.h"
#include "phBC.h"

namespace apf {
class Mesh;
}

namespace ph {

struct EnsaArrays
{
  double* coordinates;
/* describes inter-part connectivity,
   see ph::encodeILWORK */
  int* ilwork;
/* describes inter-part element connectivity,
   see ph::encodeILWORKF */
  int* ilworkf;
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
/* ibcb (part 0) has these bits:
   MF NP TV HF TW F1 F2 F3 F4
   0  1  2  3  4  5  6  7  8
   part 1 is just the value of SID */
  int*** ibcb;
/* bcb[i][j][k] is the natural boundary condition
   value for
   boundary condition k of
   element j of
   boundary block i */
/* bcb is organized as follows:
   MF NP TV     HF F1 F2 F3 F4
   0  1  2 3 4  5  6  7  8  9 */
  double*** bcb;
/* nbc[i] is the index into essential boundary condition
   arrays of local node i (probably ;) */
  int* nbc;
/* ibc[i] is the essential boundary condition
   status code of essential BC node i */
/* ibc bits are as follows:
var: rho t p u v w sc1 sc2 sc3 sc4 perio spebc
bit:  0  1 2 3 4 5  6   7   8   9    10    11 */
  int* ibc;
/* bc[i][j] is the essential boundary condition
   value of 
   essential boundary condition j of
   essential BC node i */
/* bc is organized as follows:
var:  rho t p c11 c12 c13 m1 c21 c22 c23 m2 theta sc1 sc2 sc3 sc4
idx:   0  1 2  3   4   5  6   7   8   9  10   11  12  13  14  15  */
  double** bc;
/* encodes part-local element to element connectivity */
  int* ienneigh;
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
  int nlwork; /* size of arrays.ilwork */
  int nlworkf; /* size of arrays.ilworkf */
  AllBlocks blocks;
  EnsaArrays arrays;
};

void generateOutput(Input& in, BCs& bcs, apf::Mesh* mesh, Output& o);
void writeGeomBC(Output& o, std::string path);

}

#endif
