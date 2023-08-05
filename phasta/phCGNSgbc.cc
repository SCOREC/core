#include <PCU.h>
#include "phOutput.h"
#include "phIO.h"
#include "phiotimer.h"
#include <sstream>
#include <pcu_util.h>
#include <lionPrint.h>
#include <cstdlib>

namespace ph {

// renamed, retained but not yet updated
static std::string buildCGNSgbcFileName(std::string timestep_or_dat)
{
  std::stringstream ss;
  int rank = PCU_Comm_Self() + 1;
  ss << "geombc." << timestep_or_dat << "." << rank;
  return ss.str();
}

enum {
  MAX_PARAMS = 12
};

// renamed, update is only a transpose to match CNGS.  Parallel will require mapping here or later to global numbering
void getInteriorConnectivityCGNS(Output& o, int block, apf::DynamicArray<int>& c)
{
  int nelem = o.blocks.interior.nElements[block];
  int nvert = o.blocks.interior.keys[block].nElementVertices;
  c.setSize(nelem * nvert);
  size_t i = 0;
  for (int elem = 0; elem < nelem; ++elem)
    for (int vert = 0; vert < nvert; ++vert)
      c[i++] = o.arrays.ien[block][elem][vert] + 1; /* FORTRAN indexing */
  PCU_ALWAYS_ASSERT(i == c.getSize());
}

//renamed, update is both a transpose to match CNGS and reduction to only filling the first number of vertices on the boundary whereas PHAST wanted full volume
void getBoundaryConnectivityCGNS(Output& o, int block, apf::DynamicArray<int>& c)
{
  int nelem = o.blocks.boundary.nElements[block];
// CGNS wants surface elements  int nvert = o.blocks.boundary.keys[block].nElementVertices;
  int nvert = o.blocks.boundary.keys[block].nBoundaryFaceEdges;
  c.setSize(nelem * nvert);
  size_t i = 0;
  for (int elem = 0; elem < nelem; ++elem)
    for (int vert = 0; vert < nvert; ++vert)
      c[i++] = o.arrays.ienb[block][elem][vert] + 1;
  PCU_ALWAYS_ASSERT(i == c.getSize());
}

void getInterfaceConnectivityCGNS // not extended yet other than transpose
(
  Output& o,
  int block,
  apf::DynamicArray<int>& c
)
{
  int nelem = o.blocks.interface.nElements[block];
  int nvert0 = o.blocks.interface.keys[block].nElementVertices;
  int nvert1 = o.blocks.interface.keys[block].nElementVertices1;
  c.setSize(nelem * (nvert0 + nvert1));
  size_t i = 0;
  for (int elem = 0; elem < nelem; ++elem)
    for (int vert = 0; vert < nvert0; ++vert)
      c[i++] = o.arrays.ienif0[block][elem][vert] + 1;
  for (int elem = 0; elem < nelem; ++elem)
    for (int vert = 0; vert < nvert1; ++vert)
      c[i++] = o.arrays.ienif1[block][elem][vert] + 1;
  PCU_ALWAYS_ASSERT(i == c.getSize());
}

// renamed but not updated yet
void getNaturalBCCodesCGNS(Output& o, int block, apf::DynamicArray<int>& codes)
{
  int nelem = o.blocks.boundary.nElements[block];
  codes.setSize(nelem * 2); 
  size_t i = 0;
  for (int j = 0; j < 2; ++j)
    for (int elem = 0; elem < nelem; ++elem)
      codes[i++] = o.arrays.ibcb[block][elem][j];
  PCU_ALWAYS_ASSERT(i == codes.getSize());
}

// renamed and calling the renamed functions above with output writes commented as they are PHASTA file style
void writeBlocksCGNS(FILE* f, Output& o)
{
  apf::DynamicArray<int> c;
  int params[MAX_PARAMS];
  for (int i = 0; i < o.blocks.interior.getSize(); ++i) {
    BlockKey& k = o.blocks.interior.keys[i];
    std::string phrase = getBlockKeyPhrase(k, "connectivity interior ");
    params[0] = o.blocks.interior.nElements[i];
//    fillBlockKeyParams(params, k);
    getInteriorConnectivityCGNS(o, i, c);
//    ph_write_ints(f, phrase.c_str(), &c[0], c.getSize(), 7, params);
  }
  for (int i = 0; i < o.blocks.boundary.getSize(); ++i) {
    BlockKey& k = o.blocks.boundary.keys[i];
    std::string phrase = getBlockKeyPhrase(k, "connectivity boundary ");
    params[0] = o.blocks.boundary.nElements[i];
//    fillBlockKeyParams(params, k);
    getBoundaryConnectivityCGNS(o, i, c);
//    ph_write_ints(f, phrase.c_str(), &c[0], c.getSize(), 8, params);
// this is probably the easiest path to getting the list that tells us the face (through surfID of smd) that each boundary element face is on
    phrase = getBlockKeyPhrase(k, "nbc codes ");
    apf::DynamicArray<int> codes;
    getNaturalBCCodesCGNS(o, i, codes);
//    ph_write_ints(f, phrase.c_str(), &codes[0], codes.getSize(), 8, params);
  }

}



// retaining in case it is useful but only renamed at this point
void writeCGNSgbc(Output& o, std::string path, int timestep)
{
  double t0 = PCU_Time();
  apf::Mesh* m = o.mesh;
  std::stringstream tss; 
  std::string timestep_or_dat;
  if (! timestep)
    timestep_or_dat = "dat";
  else {
    tss << timestep;   
    timestep_or_dat = tss.str();
  }
  path += buildCGNSgbcFileName(timestep_or_dat);
  phastaio_setfile(GEOMBC_WRITE);
  FILE* f = o.openfile_write(o, path.c_str());
  if (!f) {
    lion_eprint(1,"failed to open \"%s\"!\n", path.c_str());
    abort();
  }
  ph_write_preamble(f);
  int params[MAX_PARAMS];
/* all of these strings are looked for by the other programs
   reading this format, so don't fix spelling errors or
   other silliness, it has already been set in stone */
/*
  writeInt(f, "number of nodes", m->count(0));
  writeInt(f, "number of modes", o.nOverlapNodes);
  writeInt(f, "number of shapefunctions soved on processor", 0);
  writeInt(f, "number of global modes", 0);
  writeInt(f, "number of interior elements", m->count(m->getDimension()));
  writeInt(f, "number of boundary elements", o.nBoundaryElements);
  writeInt(f, "maximum number of element nodes", o.nMaxElementNodes);
  writeInt(f, "number of interior tpblocks", o.blocks.interior.getSize());
  writeInt(f, "number of boundary tpblocks", o.blocks.boundary.getSize());
  writeInt(f, "number of nodes with Dirichlet BCs", o.nEssentialBCNodes);

  params[0] = m->count(0);
  params[1] = 3;
  ph_write_doubles(f, "co-ordinates", o.arrays.coordinates,
      params[0] * params[1], 2, params);
  writeInt(f, "number of processors", PCU_Comm_Peers());
  writeInt(f, "size of ilwork array", o.nlwork);
  params[0] = m->count(0);
  writeInts(f, " mode number map from partition to global",
      o.arrays.globalNodeNumbers, m->count(0));
  writeBlocksCGNS(f, o);
  writeInts(f, "bc mapping array", o.arrays.nbc, m->count(0));
  writeInts(f, "bc codes array", o.arrays.ibc, o.nEssentialBCNodes);
  apf::DynamicArray<double> bc;
  PHASTAIO_CLOSETIME(fclose(f);)
  double t1 = PCU_Time();
  if (!PCU_Comm_Self())
    lion_oprint(1,"geombc file written in %f seconds\n", t1 - t0);
*/
}

}
