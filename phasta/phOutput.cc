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
  int n = m->count(0);
  double* x = new double[n * 3];
  apf::MeshEntity* v;
  int i = 0;
  apf::MeshIterator* it = m->begin(0);
  while ((v = m->iterate(it))) {
    apf::Vector3 p;
    m->getPoint(v, 0, p);
    for (int j = 0; j < 3; ++j)
      x[j * n + i] = p[j]; /* FORTRAN indexing */
    ++i;
  }
  m->end(it);
  assert(i == n);
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

static void getBlocks(Output& o, BCs& bcs)
{
  apf::Mesh* m = o.mesh;
  ModelBounds modelFaces;
  getBCFaces(m, bcs, modelFaces);
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

static void getBoundary(Output& o, BCs& bcs, apf::Numbering* n)
{
  apf::Mesh* m = o.mesh;
  ModelBounds modelFaces;
  getBCFaces(m, bcs, modelFaces);
  int nbc = countNaturalBCs(*o.in);
  Blocks& bs = o.blocks.boundary;
  int*** ienb = new int**[bs.getSize()];
  int*** ibcb = new int**[bs.getSize()];
  double*** bcb = new double**[bs.getSize()];
  apf::NewArray<int> js(bs.getSize());
  for (int i = 0; i < bs.getSize(); ++i) {
    ienb[i] = new int*[bs.nElements[i]];
    ibcb[i] = new int*[bs.nElements[i]];
    bcb[i] = new double*[bs.nElements[i]];
    js[i] = 0;
  }
  APF_ITERATE(ModelBounds, modelFaces, mit) {
    apf::ModelEntity* mf = *mit;
    apf::MeshEntity* f;
    apf::MeshIterator* it = m->begin(m->getDimension() - 1);
    while ((f = m->iterate(it))) {
      if (m->toModel(f) != mf)
        continue;
      BlockKey k;
      apf::MeshEntity* e = m->getUpward(f, 0);
      getBoundaryBlockKey(m, e, f, k);
      assert(bs.keyToIndex.count(k));
      int i = bs.keyToIndex[k];
      int j = js[i];
      int nv = k.nElementVertices;
      apf::Downward v;
      getBoundaryVertices(m, e, f, v);
      ienb[i][j] = new int[nv];
      for (int k = 0; k < nv; ++k)
        ienb[i][j][k] = apf::getNumber(n, v[k], 0, 0);
      ibcb[i][j] = new int[2](); /* <- parens initialize to zero */
      bcb[i][j] = new double[nbc]();
      applyNaturalBCs(m, f, bcs, bcb[i][j], ibcb[i][j]);
      ++js[i];
    }
    m->end(it);
  }
  for (int i = 0; i < bs.getSize(); ++i)
    assert(js[i] == bs.nElements[i]);
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
    iper[i] = 0;
  o.arrays.iper = iper;
}

static void getEssentialBCs(BCs& bcs, Output& o)
{
  Input& in = *o.in;
  apf::Mesh* m = o.mesh;
  int* nbc = new int[m->count(0)];
  int* ibc = new int[m->count(0)]();
  double** bc = new double*[m->count(0)];
  int nec = countEssentialBCs(in);
  double* bc_j = new double[nec]();
  size_t i = 0;
  size_t j = 0;
  apf::MeshEntity* v;
  apf::MeshIterator* it = m->begin(0);
  while ((v = m->iterate(it))) {
    int ibc_j;
    bool did = applyEssentialBCs(m, v, bcs, bc_j, &ibc_j);
    if (did) {
      nbc[i] = j;
      ibc[j] = ibc_j;
      bc[j] = bc_j;
      bc_j = new double[nec]();
      ++j;
    } else {
      nbc[i] = -1;
    }
    ++i;
  }
  m->end(it);
  delete [] bc_j;
  o.arrays.nbc = nbc;
  o.arrays.ibc = ibc;
  o.arrays.bc = bc;
  o.nEssentialBCNodes = j;
}

static void applyInitialConditions(BCs& bcs, Output& o)
{
  Input& in = *o.in;
  apf::Mesh* m = o.mesh;
  apf::MeshEntity* v;
  apf::NewArray<double> s(in.ensa_dof);
  apf::Field* f = m->findField("solution");
  apf::MeshIterator* it = m->begin(0);
  while ((v = m->iterate(it))) {
    apf::getComponents(f, v, 0, &s[0]);
    applySolutionBCs(m, v, bcs, &s[0]);
    apf::setComponents(f, v, 0, &s[0]);
  }
  m->end(it);
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
  for (int i = 0; i < nEssentialBCNodes; ++i)
    delete [] arrays.bc[i];
  delete [] arrays.bc;
}

void generateOutput(Input& in, BCs& bcs, apf::Mesh* mesh, Output& o)
{
  double t0 = MPI_Wtime();
  o.in = &in;
  o.mesh = mesh;
  getCounts(o);
  getCoordinates(o);
  getGlobal(o);
  getBlocks(o, bcs);
  apf::Numbering* n = apf::numberOverlapNodes(mesh, "ph_local");
  getLinks(o, n);
  getInterior(o, n);
  getBoundary(o, bcs, n);
  apf::destroyNumbering(n);
  getBoundaryElements(o);
  getMaxElementNodes(o);
  getFakePeriodicMasters(o);
  getEssentialBCs(bcs, o);
  applyInitialConditions(bcs, o);
  double t1 = MPI_Wtime();
  if (!PCU_Comm_Self())
    printf("generated output structs in %f seconds\n",t1 - t0);
}

}
