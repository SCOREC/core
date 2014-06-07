#include "phOutput.h"
#include "phLinks.h"
#include "phAdjacent.h"
#include <PCU.h>

namespace ph {

static void getCounts(Output& o)
{
  for (int i = 0; i < 4; ++i)
    o.nGlobalEntities[i] = apf::countOwned(o.mesh, i);
  PCU_Add_Ints(o.nGlobalEntities, 4);
  o.nOwnedNodes = apf::countOwned(o.mesh, 0);
  o.nOverlapNodes = o.mesh->count(0);
  o.nGlobalNodes = o.nGlobalEntities[0];
}

static void getCoordinates(Output& o)
{
  apf::Mesh* m = o.mesh;
  double* x = new double[m->count(0) * 3];
  apf::MeshEntity* v;
  size_t i = 0;
  apf::MeshIterator* it = m->begin(0);
  while ((v = m->iterate(it))) {
    apf::Vector3 p;
    m->getPoint(v, 0, p);
    p.toArray(&x[i]);
    i += 3;
  }
  m->end(it);
  assert(i == m->count(0) * 3);
  o.arrays.coordinates = x;
}

static void getGlobal(Output& o)
{
  apf::Mesh* m = o.mesh;
  apf::Numbering* n = apf::numberOwnedNodes(o.mesh, "ph_owned");
  apf::GlobalNumbering* gn = apf::makeGlobal(n);
  apf::synchronize(gn);
  o.arrays.globalNodeNumbers = new int[m->count(0)];
  apf::MeshEntity* v;
  size_t i = 0;
  apf::MeshIterator* it = m->begin(0);
  while ((v = m->iterate(it)))
    o.arrays.globalNodeNumbers[i++] = apf::getNumber(gn, apf::Node(v, 0));
  m->end(it);
  assert(i == m->count(0));
  apf::destroyGlobalNumbering(gn);
}

static void getBlocks(Output& o, ModelBounds& modelFaces)
{
  getAllBlocks(o.mesh, o.blocks, modelFaces);
}

static void getLinks(Output& o, apf::Numbering* n)
{
  Links links;
  getVertexLinks(o.mesh, links);
  size_t size;
  encodeLinks(n, links, size, o.arrays.ilwork);
  o.nlwork = size;
}

static void getInterior(Output& o, apf::Numbering* n)
{
  apf::Mesh* m = o.mesh;
  Blocks& bs = o.blocks.interior;
  int*** ien = new int**[bs.getSize()];
  for (int i = 0; i < bs.getSize(); ++i) {
    ien[i] = new int*[bs.nElements[i]];
    int t = bs.keys[i].elementType;
    int nv = bs.keys[i].nElementVertices;
    apf::MeshEntity* e;
    int j = 0;
    apf::MeshIterator* it = m->begin(m->getDimension());
    while ((e = m->iterate(it))) {
      if (getPhastaType(m, e) != t)
        continue;
      ien[i][j] = new int[nv];
      apf::Downward v;
      getVertices(m, e, v);
      for (int k = 0; k < nv; ++k)
        ien[i][j][k] = apf::getNumber(n, v[k], 0, 0);
      ++j;
    }
    m->end(it);
  }
  o.arrays.ien = ien;
}

static void getBoundary(Output& o, ModelBounds& modelFaces, apf::Numbering* n)
{
  apf::Mesh* m = o.mesh;
  int nbc = countNaturalBCs(*o.in);
  Blocks& bs = o.blocks.boundary;
  int*** ienb = new int**[bs.getSize()];
  int*** ibcb = new int**[bs.getSize()];
  double*** bcb = new double**[bs.getSize()];
  int i = 0;
  APF_ITERATE(ModelBounds, modelFaces, mit) {
    apf::ModelEntity* mf = *mit;
    ienb[i] = new int*[bs.nElements[i]];
    ibcb[i] = new int*[bs.nElements[i]];
    bcb[i] = new double*[bs.nElements[i]];
    int t = bs.keys[i].elementType;
    int nv = bs.keys[i].nElementVertices;
    apf::MeshEntity* f;
    int j = 0;
    apf::MeshIterator* it = m->begin(m->getDimension() - 1);
    while ((f = m->iterate(it))) {
      if (m->toModel(f) != mf)
        continue;
      apf::MeshEntity* e = m->getUpward(f, 0);
      if (getPhastaType(m, e) != t)
        continue;
      ienb[i][j] = new int[nv];
      apf::Downward v;
      getBoundaryVertices(m, e, f, v);
      for (int k = 0; k < nv; ++k)
        ienb[i][j][k] = getNumber(n, v[k], 0, 0);
      ibcb[i][j] = new int[2];
/* fake, fixme, todo */
      ibcb[i][j][0] = 42;
      ibcb[i][j][1] = 42;
      bcb[i][j] = new double[nbc];
      for (int k = 0; k < nbc; ++k)
        bcb[i][j][k] = 42.0;
      ++j;
    }
    m->end(it);
    assert(j == bs.nElements[i]);
    ++i;
  }
  assert(i == bs.getSize());
  o.arrays.ienb = ienb;
  o.arrays.ibcb = ibcb;
  o.arrays.bcb = bcb;
}

static void getBoundaryElements(Output& o)
{
  Blocks& bs = o.blocks.boundary;
  int n = 0;
  for (int i = 0; i < bs.getSize(); ++i)
    n += bs.nElements[i];
  o.nBoundaryElements = n;
}

static void getMaxElementNodes(Output& o)
{
  int n = 0;
  Blocks& ibs = o.blocks.interior;
  for (int i = 0; i < ibs.getSize(); ++i)
    n = std::max(n, ibs.keys[i].nElementVertices);
  Blocks& bbs = o.blocks.boundary;
  for (int i = 0; i < bbs.getSize(); ++i)
    n = std::max(n, bbs.keys[i].nElementVertices);
  o.nMaxElementNodes = n;
}

static void getFakePeriodicMasters(Output& o)
{
  apf::Mesh* m = o.mesh;
  int* iper = new int[m->count(0)];
  for (size_t i = 0; i < m->count(0); ++i)
    iper[i] = 42;
  o.arrays.iper = iper;
}

static void getFakeEssentialBCs(Output& o)
{
  apf::Mesh* m = o.mesh;
  o.nEssentialBCNodes = 0;
  o.arrays.nbc = new int[m->count(0)];
  for (size_t i = 0; i < m->count(0); ++i)
    o.arrays.nbc[i] = 42;
  o.arrays.ibc = 0;
  o.arrays.bc = 0;
}

Output::~Output()
{
  delete [] arrays.coordinates;
  delete [] arrays.ilwork;
  delete [] arrays.iper;
  delete [] arrays.globalNodeNumbers;
  Blocks& ibs = blocks.interior;
  for (int i = 0; i < ibs.getSize(); ++i) {
    for (int j = 0; j < ibs.nElements[i]; ++j)
      delete [] arrays.ien[i][j];
    delete [] arrays.ien[i];
  }
  delete [] arrays.ien;
  Blocks& bbs = blocks.boundary;
  for (int i = 0; i < bbs.getSize(); ++i) {
    for (int j = 0; j < bbs.nElements[i]; ++j) {
      delete [] arrays.ienb[i][j];
      delete [] arrays.ibcb[i][j];
      delete [] arrays.bcb[i][j];
    }
    delete [] arrays.ienb[i];
    delete [] arrays.ibcb[i];
    delete [] arrays.bcb[i];
  }
  delete [] arrays.ienb;
  delete [] arrays.ibcb;
  delete [] arrays.bcb;
  delete [] arrays.nbc;
  delete [] arrays.ibc;
  delete [] arrays.bc;
}

void generateOutput(Input& in, apf::Mesh* mesh, Output& o)
{
  ModelBounds modelFaces; //FIXME: BC application not done yet
  o.in = &in;
  o.mesh = mesh;
  getCounts(o);
  getCoordinates(o);
  getGlobal(o);
  getBlocks(o, modelFaces);
  apf::Numbering* n = apf::numberOverlapNodes(mesh, "ph_local");
  getLinks(o, n);
  getInterior(o, n);
  getBoundary(o, modelFaces, n);
  apf::destroyNumbering(n);
  getBoundaryElements(o);
  getMaxElementNodes(o);
  getFakePeriodicMasters(o);
  getFakeEssentialBCs(o);
}

}
