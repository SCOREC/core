#include <PCU.h>
#include "phOutput.h"
#include "phLinks.h"
#include "phAdjacent.h"
#include "phBubble.h"

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

/* so apparently old phParAdapt just used EN_id,
   and the id generator from pumi would do things
   like this. I guess PHASTA is ok with a unique
   number for each copy, regardless of part boundary
   sharing...
 
update: Michel says these global numbers are ignored
        by phasta. get rid of them when you can.
 */
static void getGlobal(Output& o)
{
  apf::Mesh* m = o.mesh;
  int n = m->count(0);
  int self = PCU_Comm_Self();
  int peers = PCU_Comm_Peers();
  int id = self + 1;
  o.arrays.globalNodeNumbers = new int[n];
  for (int i = 0; i < n; ++i) {
    o.arrays.globalNodeNumbers[i] = id;
    id += peers;
  }
}

static void getVertexLinks(Output& o, apf::Numbering* n)
{
  Links links;
  getLinks(o.mesh, 0, links);
  encodeILWORK(n, links, o.nlwork, o.arrays.ilwork);
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
  gmi_model* gm = m->getModel();
  gmi_ent* gf;
  gmi_iter* git = gmi_begin(gm, m->getDimension() - 1);
  while ((gf = gmi_next(gm, git))) {
    apf::ModelEntity* mf = (apf::ModelEntity*)gf;
    int* ibcbMaster = new int[2]();
    double* bcbMaster = new double[nbc]();
    applyNaturalBCs(gm, gf, bcs, bcbMaster, ibcbMaster);
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
      bcb[i][j] = new double[nbc]();
      for (int k = 0; k < nbc; ++k)
        bcb[i][j][k] = bcbMaster[k];
      ibcb[i][j] = new int[2](); /* <- parens initialize to zero */
      for (int k = 0; k < 2; ++k)
        ibcb[i][j][k] = ibcbMaster[k];
      ++js[i];
    }
    delete [] ibcbMaster;
    delete [] bcbMaster;
    m->end(it);
  }
  gmi_end(gm, git);
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

static bool matchLess(apf::Copy& a, apf::Copy& b)
{
  if (a.peer != b.peer)
    return a.peer < b.peer;
  return a.entity < b.entity;
}

/* returns the global periodic master iff it is on this
   part, otherwise returns e */
static apf::MeshEntity* getPeriodicMaster(apf::Mesh* m, apf::MeshEntity* e)
{
  if ( ! m->hasMatching())
    return e;
  apf::Matches matches;
  m->getMatches(e, matches);
  if (!matches.getSize())
    return e;
  int self = PCU_Comm_Self();
  apf::Copy master(self, e);
  for (size_t i = 0; i < matches.getSize(); ++i)
    if (matchLess(matches[i], master))
      master = matches[i];
  if (master.peer == self)
    return master.entity;
  return e;
}

static void getPeriodicMasters(Output& o, apf::Numbering* n)
{
  apf::Mesh* m = o.mesh;
  int* iper = new int[m->count(0)];
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* e;
  int i = 0;
  while ((e = m->iterate(it))) {
    apf::MeshEntity* master = getPeriodicMaster(m, e);
    if (master == e)
      iper[i] = 0;
    else
      iper[i] = apf::getNumber(n, master, 0, 0) + 1;
    ++i;
  }
  m->end(it);
  o.arrays.iper = iper;
}

static void getEssentialBCsOn(BCs& bcs, Output& o, gmi_ent* ge)
{
  Input& in = *o.in;
  apf::Mesh* m = o.mesh;
  int ibcMaster = 0;
  int nec = countEssentialBCs(in);
  double* bcMaster = new double[nec]();
  gmi_model* gm = m->getModel();
  bool did = applyEssentialBCs(gm, ge, bcs, bcMaster, &ibcMaster);
  /* matching introduces an iper bit which in our system
     is really a per-entity thing, not specifically dictated
     by classification, so in that case we have to look
     at all the vertices anyway to see if they are periodic slaves */
  if (did || m->hasMatching()) {
    apf::MeshEntity* v;
    apf::MeshIterator* it = m->begin(0);
    int i = 0;
    int& ei = o.nEssentialBCNodes;
    while ((v = m->iterate(it))) {
      if (m->toModel(v) == (apf::ModelEntity*)ge) {
        apf::MeshEntity* master = getPeriodicMaster(m, v);
        if (did || (master != v)) {
          o.arrays.nbc[i] = ei + 1;
          o.arrays.ibc[ei] = ibcMaster;
          if (master != v)
            o.arrays.ibc[ei] |= (1<<10); //yes, hard coded...
          double* bc_ei = new double[nec]();
          for (int j = 0; j < nec; ++j)
            bc_ei[j] = bcMaster[j];
          o.arrays.bc[ei] = bc_ei;
          ++ei;
        }
      }
      ++i;
    }
    m->end(it);
  }
  delete [] bcMaster;
}

static void getEssentialBCs(BCs& bcs, Output& o)
{
  apf::Mesh* m = o.mesh;
  int nv = m->count(0);
  o.arrays.nbc = new int[nv];
  for (int i = 0; i < nv; ++i)
    o.arrays.nbc[i] = 0;
  o.arrays.ibc = new int[nv]();
  o.arrays.bc = new double*[nv];
  o.nEssentialBCNodes = 0;
  gmi_model* gm = m->getModel();
  for (int d = 3; d >= 0; --d) {
    gmi_iter* it = gmi_begin(gm, d);
    gmi_ent* ge;
    while ((ge = gmi_next(gm, it)))
      getEssentialBCsOn(bcs, o, ge);
    gmi_end(gm, it);
  }
}

static void getInitialConditions(BCs& bcs, Output& o)
{
  Input& in = *o.in;
  apf::Mesh* m = o.mesh;
  apf::MeshEntity* v;
  apf::NewArray<double> s(in.ensa_dof);
  apf::Field* f = m->findField("solution");
  apf::MeshIterator* it = m->begin(0);
  gmi_model* gm = m->getModel();
  while ((v = m->iterate(it))) {
    gmi_ent* ge = (gmi_ent*)m->toModel(v);
    apf::getComponents(f, v, 0, &s[0]);
/* unfortunately, there is no way to know which
   components this overwrites without creating
   bad coupling between pieces of code, so we
   call this for every vertex even though it
   only depends on classification */
    applySolutionBCs(gm, ge, bcs, &s[0]);
    apf::setComponents(f, v, 0, &s[0]);
  }
  m->end(it);
}

static void getElementGraph(Output& o)
{
  if (o.in->formElementGraph) {
    apf::Numbering* n = apf::numberElements(o.mesh, "ph::getElementGraph");
    o.arrays.ienneigh = formIENNEIGH(n);
    Links links;
    getLinks(o.mesh, o.mesh->getDimension() - 1, links);
    encodeILWORKF(n, links, o.nlworkf, o.arrays.ilworkf);
  } else {
    o.arrays.ilworkf = 0;
    o.arrays.ienneigh = 0;
  }
}

Output::~Output()
{
  delete [] arrays.coordinates;
  delete [] arrays.ilwork;
  delete [] arrays.ilworkf;
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
  delete [] arrays.ienneigh;
}

void generateOutput(Input& in, BCs& bcs, apf::Mesh* mesh, Output& o)
{
  double t0 = PCU_Time();
  o.in = &in;
  o.mesh = mesh;
  getCounts(o);
  getCoordinates(o);
  getGlobal(o);
  getAllBlocks(o.mesh, o.blocks);
  apf::Numbering* n = apf::numberOverlapNodes(mesh, "ph_local");
  getVertexLinks(o, n);
  getInterior(o, n);
  getBoundary(o, bcs, n);
  getPeriodicMasters(o, n);
  apf::destroyNumbering(n);
  getBoundaryElements(o);
  getMaxElementNodes(o);
  getEssentialBCs(bcs, o);
  getInitialConditions(bcs, o);
  getElementGraph(o);
  if (in.initBubbles)
    initBubbles(o.mesh, in);
  double t1 = PCU_Time();
  if (!PCU_Comm_Self())
    printf("generated output structs in %f seconds\n",t1 - t0);
}

}
