#include "phOutput.h"
#include "phIO.h"
#include <sstream>
#include <PCU.h>

namespace ph {

static std::string buildGeomBCFileName(Input& in)
{
  std::stringstream ss;
  int rank = PCU_Comm_Self() + 1;
  ss << "geombc.dat." << rank;
  return ss.str();
}

/* note: doing this as an "int" will almost certainly
   overflow for billion-element meshes.
   but, thats the format used in previous geombc. */
int getTotalNodes(apf::Mesh* m)
{
  int n = apf::countOwned(m, 0);
  PCU_Add_Ints(&n, 1);
  return n;
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
      c[i++] = o.arrays.ien[block][elem][vert];
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
      c[i++] = o.arrays.ienb[block][elem][vert];
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

int countNaturalBCs(Input& in)
{
  return in.ensa_dof + 1;
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
    params[0] = o.blocks.interior.nElements[i];
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

void writeGeomBC(Output& o, std::string path)
{
  apf::Mesh* m = o.mesh;
  path += buildGeomBCFileName(*o.in);
  FILE* f = fopen(path.c_str(), "w");
  ph_write_preamble(f);
  int params[MAX_PARAMS];
  params[0] = m->count(0);
  ph_write_header(f, "number of nodes", 0, 1, params);
  params[0] = o.nGlobalEntities[0];
  ph_write_header(f, "number of nodes in the mesh", 0, 1, params);
  params[0] = o.nGlobalEntities[1];
  ph_write_header(f, "number of edges in the mesh", 0, 1, params);
  params[0] = o.nGlobalEntities[2];
  ph_write_header(f, "number of faces in the mesh", 0, 1, params);
  params[0] = o.nOverlapNodes;
  ph_write_header(f, "number of modes", 0, 1, params);
  params[0] = o.nOwnedNodes;
  ph_write_header(f, "number of shapefunctions soved on processor", 0, 1, params);
  /* keep the bad spelling ! ------------------- ^ */
  params[0] = o.nGlobalNodes;
  ph_write_header(f, "number of global modes", 0, 1, params);
  params[0] = m->count(m->getDimension());
  ph_write_header(f, "number of interior elements", 0, 1, params);
  params[0] = o.nBoundaryElements;
  ph_write_header(f, "number of boundary elements", 0, 1, params);
  params[0] = o.nMaxElementNodes;
  ph_write_header(f, "maximum number of element nodes", 0, 1, params);
  params[0] = o.blocks.interior.getSize();
  ph_write_header(f, "number of interior tpblocks", 0, 1, params);
  params[0] = o.blocks.boundary.getSize();
  ph_write_header(f, "number of boundary tpblocks", 0, 1, params);
  params[0] = o.nDirichletNodes;
  ph_write_header(f, "number of nodes with Dirichlet BCs", 0, 1, params);
  params[0] = m->count(0);
  params[1] = 3;
  ph_write_doubles(f, "co-ordinates", o.arrays.coordinates,
      params[0] * params[1], 2, params);
  params[0] = PCU_Comm_Peers();
  ph_write_header(f, "number of processors", 0, 1, params);
  params[0] = o.nlwork;
  ph_write_header(f, "size of ilwork array", 0, 1, params);
  if (o.nlwork)
    ph_write_ints(f, "ilwork", o.arrays.ilwork, o.nlwork, 1, params);
  params[0] = m->count(0);
  ph_write_ints(f, " mode number map from partition to global",
      o.arrays.globalNodeNumbers, m->count(0), 1, params);
  writeBlocks(f, o);
  fclose(f);
}

}
