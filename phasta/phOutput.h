#ifndef PH_OUTPUT_H
#define PH_OUTPUT_H

#include "phInput.h"
#include "phBlock.h"
#include "phBC.h"

namespace apf {
class Mesh;
}

struct GRStream;

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
/* ienif0,1[i][j][k] are the local vertex id of
   vertex k of
   element j of
   interface block i
   ienif0 and ienif1 correspond to the two elements on the interface
 */
  int*** ienif0;
  int*** ienif1;
/* mattype[i][j] is the material type of
   element j of
   interior block i */
  int** mattype;
/* mattypeb[i][j] is the material type of
   element j of
   boundary block i */
  int** mattypeb;
/* mattypeif0,1[i][j] are the material types of
   element j of
   interface blocks 0,1 i */
  int** mattypeif0;
  int** mattypeif1;
/* ibcb[i][j][k] is the natural boundary condition
   status code
   number k in [0,1] of
   element j of
   boundary block i */
/* ibcb (part 0) has these bits:
   MF NP TV HF TW F1 F2 F3 F4 TVM
   0  1  2  3  4  5  6  7  8   9
   part 1 is just the value of SID */
  int*** ibcb;
/* bcb[i][j][k] is the natural boundary condition
   value for
   boundary condition k of
   element j of
   boundary block i */
/* bcb is organized as follows:
   MF NP TV     HF F1 F2 F3 F4 --TVM---
   0  1  2 3 4  5  6  7  8  9  10 11 12 */
  double*** bcb;
/* nbc[i] is the index into essential boundary condition
   arrays of local node i (probably ;) */
  int* nbc;
/* ibc[i] is the essential boundary condition
   status code of essential BC node i */
/* ibc bits are as follows:
var: rho t p u v w sc1 sc2 sc3 sc4 perio --NULL-- eu ev ew
bit:  0  1 2 3 4 5  6   7   8   9    10  11 12 13 14 15 16 */
  int* ibc;
/* bc[i][j] is the essential boundary condition
   value of
   essential boundary condition j of
   essential BC node i */
/* bc is organized as follows:
var:  rho t p c11 c12 c13 m1 c21 c22 c23 m2 theta sc1 sc2 sc3 sc4 ec11 ec12 ec13 em1 ec21 ec22 ec23 em2
idx:   0  1 2  3   4   5  6   7   8   9  10   11  12  13  14  15   16   17   18  19   20   21   22  23  */
  double** bc;
/* encodes part-local element to element connectivity */
  int* ienneigh;
/* encodes parallel communication between mesh edges */
  int* ilworkl;
/* tetrahedra to local edges, already in FORTRAN order and from 1 */
  int* iel;
/* edges to tetrahedra offsets */
  int* ileo;
/* edges to tetrahedra adjacencies */
  int* ile;
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
  int nBoundaryElements;
  int nInterfaceElements;
  int nMaxElementNodes;
  int nEssentialBCNodes;
  int nOverlapEdges;
  int nlwork; /* size of arrays.ilwork */
  int nlworkf; /* size of arrays.ilworkf */
  int nlworkl; /* size of arrays.ilworkl */
  bool hasDGInterface;
  FILE* (*openfile_write)(Output& out, const char* path);
  GRStream* grs;
  AllBlocks blocks;
  EnsaArrays arrays;
};

void generateOutput(Input& in, BCs& bcs, apf::Mesh* mesh, Output& o);
void writeGeomBC(Output& o, std::string path, int timestep_or_dat = 0);

}

#endif
