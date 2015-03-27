#include <PCU.h>
#include "phOutput.h"
#include "phIO.h"
#include <sstream>

namespace ph {

static std::string buildGeomBCFileName()
{
  std::stringstream ss;
  int rank = PCU_Comm_Self() + 1;
  ss << "geombc.dat." << rank;
  return ss.str();
}

enum {
  MAX_PARAMS = 8
};

void getInteriorConnectivity(Output& o, int block, apf::DynamicArray<int>& c)
{
  int nelem = o.blocks.interior.nElements[block];
  int nvert = o.blocks.interior.keys[block].nElementVertices;
  c.setSize(nelem * nvert);
  size_t i = 0;
  for (int vert = 0; vert < nvert; ++vert)
    for (int elem = 0; elem < nelem; ++elem)
      c[i++] = o.arrays.ien[block][elem][vert] + 1; /* FORTRAN indexing */
  assert(i == c.getSize());
}

void getBoundaryConnectivity(Output& o, int block, apf::DynamicArray<int>& c)
{
  int nelem = o.blocks.boundary.nElements[block];
  int nvert = o.blocks.boundary.keys[block].nElementVertices;
  c.setSize(nelem * nvert);
  size_t i = 0;
  for (int vert = 0; vert < nvert; ++vert)
    for (int elem = 0; elem < nelem; ++elem)
      c[i++] = o.arrays.ienb[block][elem][vert] + 1;
  assert(i == c.getSize());
}

void getNaturalBCCodes(Output& o, int block, apf::DynamicArray<int>& codes)
{
  int nelem = o.blocks.boundary.nElements[block];
  codes.setSize(nelem * 2);
  size_t i = 0;
  for (int j = 0; j < 2; ++j)
    for (int elem = 0; elem < nelem; ++elem)
      codes[i++] = o.arrays.ibcb[block][elem][j];
  assert(i == codes.getSize());
}

void getNaturalBCValues(Output& o, int block, apf::DynamicArray<double>& values)
{
  int nelem = o.blocks.boundary.nElements[block];
  int nbc = countNaturalBCs(*o.in);
  values.setSize(nelem * nbc);
  size_t i = 0;
  for (int bc = 0; bc < nbc; ++bc)
    for (int elem = 0; elem < nelem; ++elem)
      values[i++] = o.arrays.bcb[block][elem][bc];
  assert(i == values.getSize());
}

void getEssentialBCValues(Output& o, apf::DynamicArray<double>& values)
{
  int nnode = o.nEssentialBCNodes;
  int nbc = countEssentialBCs(*o.in);
  values.setSize(nnode * nbc);
  size_t i = 0;
  for (int bc = 0; bc < nbc; ++bc)
    for (int node = 0; node < nnode; ++node)
      values[i++] = o.arrays.bc[node][bc];
  assert(i == values.getSize());
}

void fillBlockKeyParams(int* params, BlockKey& k)
{
  params[1] = k.nElementVertices;
  params[2] = k.polynomialOrder;
  params[3] = k.nElementVertices; /* num nodes */
  params[4] = k.nBoundaryFaceEdges; /* num boundary nodes */
  params[5] = k.nBoundaryFaceEdges;
  params[6] = k.elementType;
}

void writeBlocks(FILE* f, Output& o)
{
  apf::DynamicArray<int> c;
  int params[MAX_PARAMS];
  for (int i = 0; i < o.blocks.interior.getSize(); ++i) {
    BlockKey& k = o.blocks.interior.keys[i];
    std::string phrase = getBlockKeyPhrase(k, "connectivity interior ");
    params[0] = o.blocks.interior.nElements[i];
    fillBlockKeyParams(params, k);
    getInteriorConnectivity(o, i, c);
    ph_write_ints(f, phrase.c_str(), &c[0], c.getSize(), 7, params);
  }
  for (int i = 0; i < o.blocks.boundary.getSize(); ++i) {
    BlockKey& k = o.blocks.boundary.keys[i];
    std::string phrase = getBlockKeyPhrase(k, "connectivity boundary ");
    params[0] = o.blocks.boundary.nElements[i];
    fillBlockKeyParams(params, k);
    params[7] = countNaturalBCs(*o.in);
    getBoundaryConnectivity(o, i, c);
    ph_write_ints(f, phrase.c_str(), &c[0], c.getSize(), 8, params);
    phrase = getBlockKeyPhrase(k, "nbc codes ");
    apf::DynamicArray<int> codes;
    getNaturalBCCodes(o, i, codes);
    ph_write_ints(f, phrase.c_str(), &codes[0], codes.getSize(), 8, params);
    phrase = getBlockKeyPhrase(k, "nbc values ");
    apf::DynamicArray<double> values;
    getNaturalBCValues(o, i, values);
    ph_write_doubles(f, phrase.c_str(), &values[0], values.getSize(), 8, params);
  }
}

static void writeInt(FILE* f, const char* name, int i)
{
  ph_write_header(f, name, 0, 1, &i);
}

static void writeInts(FILE* f, const char* name, int* i, int n)
{
  ph_write_ints(f, name, i, n, 1, &n);
}

static void writeDoubles(FILE* f, const char* name, double* d, int n)
{
  ph_write_doubles(f, name, d, n, 1, &n);
}

static void writeElementGraph(Output& o, FILE* f)
{
  if (o.in->formElementGraph) {
    apf::Mesh* m = o.mesh;
    int dim = m->getDimension();
    int type = getFirstType(m, dim);
    int nsides = apf::Mesh::adjacentCount[type][dim - 1];
    size_t nelem = m->count(dim);
    writeInt(f, "size of ilworkf array", o.nlworkf);
    writeInts(f, "ilworkf", o.arrays.ilworkf, o.nlworkf);
    writeInts(f, "ienneigh", o.arrays.ienneigh, nelem * nsides);
  }
}

void writeGeomBC(Output& o, std::string path)
{
  double t0 = PCU_Time();
  apf::Mesh* m = o.mesh;
  path += buildGeomBCFileName();
  FILE* f = fopen(path.c_str(), "w");
  if (!f) {
    fprintf(stderr,"failed to open \"%s\"!\n", path.c_str());
    abort();
  }
  ph_write_preamble(f);
  int params[MAX_PARAMS];
/* all of these strings are looked for by the other programs
   reading this format, so don't fix spelling errors or
   other silliness, it has already been set in stone */
  writeInt(f, "number of nodes", m->count(0));
  writeInt(f, "number of nodes in the mesh", o.nGlobalEntities[0]);
  writeInt(f, "number of edges in the mesh", o.nGlobalEntities[1]);
  writeInt(f, "number of faces in the mesh", o.nGlobalEntities[2]);
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
  if (o.nlwork)
    writeInts(f, "ilwork ", o.arrays.ilwork, o.nlwork);
  params[0] = m->count(0);
  writeInts(f, " mode number map from partition to global",
      o.arrays.globalNodeNumbers, m->count(0));
  writeBlocks(f, o);
  writeInts(f, "bc mapping array", o.arrays.nbc, m->count(0));
  writeInts(f, "bc codes array", o.arrays.ibc, o.nEssentialBCNodes);
  apf::DynamicArray<double> bc;
  getEssentialBCValues(o, bc);
  writeDoubles(f, "boundary condition array", &bc[0], bc.getSize());
  writeInts(f, "periodic masters array", o.arrays.iper, m->count(0));
  writeElementGraph(o, f);
  fclose(f);
  double t1 = PCU_Time();
  if (!PCU_Comm_Self())
    printf("geombc file written in %f seconds\n", t1 - t0);
}

}
