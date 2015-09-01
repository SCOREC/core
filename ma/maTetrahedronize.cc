#include <PCU.h>
#include "maTetrahedronize.h"
#include "maCrawler.h"
#include "maAdapt.h"
#include "maRefine.h"
#include "maLayer.h"
#include <apfNumbering.h>
#include <apfShape.h>
#include <apfCavityOp.h>
#include "maShape.h"

namespace ma {

int getDiagonalFromFlag(Adapt* a, Entity* e)
{
  if (getFlag(a,e,DIAGONAL_1))
    return 0;
  if (getFlag(a,e,DIAGONAL_2))
    return 1;
  return -1;
}

int getFlagFromDiagonal(int diagonal)
{
  if (diagonal==0)
    return DIAGONAL_1;
  assert(diagonal==1);
  return DIAGONAL_2;
}

/* when triangulating a stack of quads that grows
   from one surface edge, the direction of all
   the diagonals for those quads is dictated
   by the global numbers of the surface edge vertices.
   this prevents cyclic diagonals around prisms */
static bool getEdgeDirection(apf::GlobalNumbering* n, Entity* e)
{
  apf::Mesh* m = getMesh(n);
  Entity* v[2];
  m->getDownward(e,0,v);
  return apf::getNumber(n, apf::Node(v[0], 0))
         <
         apf::getNumber(n, apf::Node(v[1], 0));
}

/* globally number the layer base vertices,
   then choose directions based on the numbering.
   this prevents cycles around layer base triangles
   and, by extension, layer prisms */
static void chooseBaseDiagonals(Adapt* a)
{
  Mesh* m = a->mesh;
  apf::Numbering* local =
    apf::numberOwnedDimension(m, "layer_base_number", 0);
  apf::GlobalNumbering* global = apf::makeGlobal(local);
  apf::synchronize(global);
  Entity* e;
  Iterator * it = m->begin(1);
  while ((e = m->iterate(it)))
    if (getFlag(a,e,LAYER_BASE))
    {
      if (getEdgeDirection(global, e))
        setFlag(a,e,DIAGONAL_1);
      else
        setFlag(a,e,DIAGONAL_2);
    }
  m->end(it);
  apf::destroyGlobalNumbering(global);
}

/* Crawler helper: given a layer edge, find the
   adjacent not-yet-visited quad */
static Entity* getOtherQuad(Adapt* a, Entity* e, Predicate& visited)
{
  Mesh* m = a->mesh;
  apf::Up up;
  m->getUp(e,up);
  for (int i=0; i < up.n; ++i)
  {
    Entity* of = up.e[i];
    if ((m->getType(of) == apf::Mesh::QUAD)&&
        ( ! visited(of)))
      return of;
  }
  return 0;
}

/* retrieve a bit that indicates the
   relative orientation of a layer
   edge and its adjacent quad, as it
   applies to propagating diagonal
   flags up the stack (either the
   diagonal direction is negated or not)
   this function is given the quadEdges
   and directions directly to avoid
   re-computing them */
static int getQuadEdgeDiagonalBit(
    Entity* edge,
    Entity** quadEdges,
    int* directions)
{
  int i = findIn(quadEdges,4,edge);
  int i_bit = i & 1;
  int dir_bit = directions[i];
  return i_bit ^ dir_bit;
}

/* given a layer edge with a diagonal
   flag on it and the "previous" quad,
   propagate the flag onto the next
   quad and the next edge up the stack */
static Entity* flagQuad(Adapt* a, Entity* q, Entity* e)
{
  Mesh* m = a->mesh;
  int diagonal = getDiagonalFromFlag(a,e);
  assert(diagonal != -1);
  Entity* es[4];
  int ds[4];
  getFaceEdgesAndDirections(m,q,es,ds);
  diagonal ^= getQuadEdgeDiagonalBit(e,es,ds);
  setFlag(a, q, getFlagFromDiagonal(diagonal));
  e = getQuadEdgeOppositeEdge(m,q,e);
  /* bit flip going out is the opposite of bit flip
     going in   V   */
  diagonal ^= 1 ^ getQuadEdgeDiagonalBit(e,es,ds);
  setFlag(a, e, getFlagFromDiagonal(diagonal));
  return e;
}

/* this Crawler propagates diagonal flags
   from layer base edges onto all horizontal
   layer edges and all layer quads */
struct QuadFlagger : public Crawler
{
  QuadFlagger(Adapt* a):Crawler(a) {}
  void begin(Layer& first)
  {
    getDimensionBase(adapter, 1, first);
  }
  Entity* crawl(Entity* e)
  {
    HasFlag p(adapter, DIAGONAL_1 | DIAGONAL_2);
    Entity* q = getOtherQuad(adapter, e, p);
    Entity* e2 = 0;
    if (q)
      e2 = flagQuad(adapter, q, e);
    clearFlag(adapter, e, DIAGONAL_1 | DIAGONAL_2);
    return e2;
  }
  void send(Entity* e, int to)
  {
    int diagonal = getDiagonalFromFlag(adapter, e);
    PCU_COMM_PACK(to, diagonal);
  }
  bool recv(Entity* e, int)
  {
    int diagonal;
    PCU_COMM_UNPACK(diagonal);
    if (getFlag(adapter, e, DIAGONAL_1 | DIAGONAL_2))
      return false;
    setFlag(adapter, e, getFlagFromDiagonal(diagonal));
    return true;
  }
};

static void flagQuadDiagonals(Adapt* a)
{
  QuadFlagger op(a);
  crawlLayers(&op);
}

/* as it turns out, the layer refinement templates
   can create elements which are unsafe to tetrahedronize
   even if the input elements are safe.

   this seems to be a fundamental geometric property
   that is hard to prevent, so we'll start cleaning up after it.

   this class identifies pyramids which are unsafe to
   tetrahedronize and tries to change the diagonal
   flag on their quad to be the safe one, if there is one.

   it also prints warnings when it identifies cases
   that it can't help with.

   the CHECKED flag is used both to prevent visiting
   elements twice in the CavityOp and to mark quads
   whose diagonals have been overridden, to identify
   cases of conflicting overrides. */

static int align_quad(apf::Mesh* mesh, Entity* quad, Entity* canonical_0_vert)
{
  Entity* qv[4];
  mesh->getDownward(quad, 0, qv);
  int dvi = apf::findIn(qv, 4, canonical_0_vert);
  return dvi % 2;
}

static void overrideQuadDiagonal(Adapt* a, Entity* quad, int diagonal)
{
  assert(diagonal == 1 || diagonal == 0);
  Mesh* mesh = a->mesh;
  int old_diagonal = getDiagonalFromFlag(a, quad);
  if (old_diagonal != diagonal) {
    if (getFlag(a, quad, CHECKED)) {
      std::stringstream ss;
      ss << "quad at " << apf::getLinearCentroid(mesh, quad)
        << " has conflicting overrides on its diagonal.\n";
      ss << "a negative tet WILL get produced here.\n";
      std::string s = ss.str();
      fprintf(stderr,"%s",s.c_str());
    } else {
      std::stringstream ss;
      ss << "overriding diagonal at " << apf::getLinearCentroid(mesh, quad) << '\n';
      std::string s = ss.str();
      fprintf(stderr,"%s",s.c_str());
      int flag = getFlagFromDiagonal(diagonal);
      int old_flag = getFlagFromDiagonal(old_diagonal);
      clearFlag(a, quad, old_flag);
      setFlag(a, quad, flag);
      setFlag(a, quad, CHECKED);
      assert(getFlag(a, quad, flag));
      assert(!getFlag(a, quad, old_flag));
    }
  } else {
    std::stringstream ss;
    ss << "diagonal at " << apf::getLinearCentroid(mesh, quad)
       << " had a consistent override.\n";
    std::string s = ss.str();
    fprintf(stderr,"%s",s.c_str());
  }
}

struct UnsafePyramidOverride : public apf::CavityOp
{
  UnsafePyramidOverride(Adapt* a_):
    apf::CavityOp(a_->mesh)
  {
    a = a_;
  }
  Outcome setEntity(Entity* r)
  {
    if ((mesh->getType(r) != apf::Mesh::PYRAMID) ||
        getFlag(a, r, CHECKED))
      return SKIP;
    pyramid = r;
    bool isOk = isPyramidOk(mesh, pyramid, &good_rotation);
    if (isOk) {
      setFlag(a, pyramid, CHECKED);
      return SKIP;
    }
    if (good_rotation == -1) {
      setFlag(a, pyramid, CHECKED);
      std::stringstream ss;
      ss << "pyramid at " << apf::getLinearCentroid(mesh, pyramid)
         << " has no good rotation!\n";
      ss << "a negative tet WILL get produced here.\n";
      std::string s = ss.str();
      fprintf(stderr,"%s",s.c_str());
      return SKIP;
    }
    Entity* faces[5];
    mesh->getDownward(pyramid, 2, faces);
    quad = faces[0];
    if (!requestLocality(&quad, 1))
      return REQUEST;
    return OK;
  }
  void apply()
  {
    Entity* pv[5];
    mesh->getDownward(pyramid, 0, pv);
    Entity* dv = 0;
    if (good_rotation == 0)
      dv = pv[0];
    else
      dv = pv[1];
    int diagonal = align_quad(mesh, quad, dv);
    setFlag(a, pyramid, CHECKED);
    overrideQuadDiagonal(a, quad, diagonal);
  }
  Adapt* a;
  int good_rotation;
  Entity* pyramid;
  Entity* quad;
};

struct UnsafePrismOverride : public apf::CavityOp
{
  UnsafePrismOverride(Adapt* a_):
    apf::CavityOp(a_->mesh)
  {
    a = a_;
  }
  Outcome setEntity(Entity* r)
  {
    if ((mesh->getType(r) != apf::Mesh::PRISM) ||
        getFlag(a, r, CHECKED))
      return SKIP;
    prism = r;
    bool isOk = isPrismOk(mesh, prism, &good_diagonal_codes);
    if (isOk) {
      setFlag(a, prism, CHECKED);
      return SKIP;
    }
    if (good_diagonal_codes == 0) {
      setFlag(a, prism, CHECKED);
      std::stringstream ss;
      ss << "prism at " << apf::getLinearCentroid(mesh, prism)
         << " has no good diagonals!\n";
      ss << "a negative tet WILL get produced here.\n";
      std::string s = ss.str();
      fprintf(stderr,"%s",s.c_str());
      return SKIP;
    }
    Entity* faces[5];
    mesh->getDownward(prism, 2, faces);
    for (int i = 0; i < 3; ++i)
      quads[i] = faces[i + 1];
    if (!requestLocality(quads, 3))
      return REQUEST;
    return OK;
  }
  int getAllowedDiagonals()
  {
    int allowed_diagonals = 0;
    Entity* pv[6];
    mesh->getDownward(prism, 0, pv);
    for (int i = 0; i < 3; ++i) {
      if (!getFlag(a, quads[i], CHECKED)) {
        allowed_diagonals |= (3 << (2 * i));
        continue;
      }
      int quad_diagonal = getDiagonalFromFlag(a, quads[i]);
      Entity* dv = pv[i];
      int quad_alignment = align_quad(mesh, quads[i], dv);
      int diagonal = quad_diagonal ^ quad_alignment;
      allowed_diagonals |= (1 << (2 * i + diagonal));
    }
    assert(allowed_diagonals > 0);
    assert(allowed_diagonals <= ((1<<6)-1));
    return allowed_diagonals;
  }
  bool areDiagonalsAllowed(int diagonals, int allowed_diagonals)
  {
    int is_good = 1 & (good_diagonal_codes >> diagonals);
    if (!is_good)
      return false;
    for (int i = 0; i < 3; ++i) {
      int diagonal = 1 & (diagonals >> i);
      int allowed = 1 & (allowed_diagonals >> (2 * i + diagonal));
      if (!allowed)
        return false;
    }
    return true;
  }
  void enforceDiagonals(int diagonals)
  {
    Entity* pv[6];
    mesh->getDownward(prism, 0, pv);
    for (int i = 0; i < 3; ++i) {
      Entity* dv = pv[i];
      int quad_diagonal = 1 & (diagonals >> i);
      Entity* quad = quads[i];
      int quad_alignment = align_quad(mesh, quad, dv);
      int diagonal = quad_diagonal ^ quad_alignment;
      overrideQuadDiagonal(a, quad, diagonal);
    }
  }
  void apply()
  {
    int allowed_diagonals = getAllowedDiagonals();
    for (int diagonals = 1; diagonals < 7; ++diagonals)
      if (areDiagonalsAllowed(diagonals, allowed_diagonals)) {
        enforceDiagonals(diagonals);
        return;
      }
    { 
      std::stringstream ss;
      ss << "prism at " << apf::getLinearCentroid(mesh, prism)
         << " has no safe acyclic diagonals\n";
      ss << "will try cyclic diagonals\n";
      std::string s = ss.str();
      fprintf(stderr, "%s", s.c_str());
    }
    if (areDiagonalsAllowed(0, allowed_diagonals)) {
      enforceDiagonals(0);
      return;
    }
    if (areDiagonalsAllowed(7, allowed_diagonals)) {
      enforceDiagonals(7);
      return;
    } 
    { 
      std::stringstream ss;
      ss << "prism at " << apf::getLinearCentroid(mesh, prism)
         << " has no safe diagonals!\n";
      ss << "A negative tet WILL get made here\n";
      std::string s = ss.str();
      fprintf(stderr, "%s", s.c_str());
    }
  }
  Adapt* a;
  Entity* prism;
  Entity* quads[3];
  int good_diagonal_codes;
};

static void overrideDiagonalsForUnsafeElements(Adapt* a)
{
  UnsafePyramidOverride op(a);
  op.applyToDimension(3);
  UnsafePrismOverride op2(a);
  op2.applyToDimension(3);
  clearFlagFromDimension(a, CHECKED, 2);
  clearFlagFromDimension(a, CHECKED, 3);
}

static void prepareLayerToTets(Adapt* a)
{
  findLayerBase(a);
  chooseBaseDiagonals(a);
  flagQuadDiagonals(a);
  overrideDiagonalsForUnsafeElements(a);
}

/* a tetrahedronization operation re-uses
   the refinement machinery (it is actually
   refinement if you think about it, parent
   elements are split into child elements).
   this function fills the refinement containers
   will all layer elements */
static void addAllLayerElements(Refine* r)
{
  Adapt* a = r->adapt;
  Mesh* m = a->mesh;
  int quadCount = apf::countEntitiesOfType(m, apf::Mesh::QUAD);
  int prismCount = apf::countEntitiesOfType(m, apf::Mesh::PRISM);
  int pyramidCount = apf::countEntitiesOfType(m, apf::Mesh::PYRAMID);
  r->toSplit[2].setSize(quadCount);
  r->toSplit[3].setSize(prismCount + pyramidCount);
  Entity* e;
  Iterator* it = m->begin(2);
  int nf = 0;
  while ((e = m->iterate(it)))
    if (m->getType(e) == apf::Mesh::QUAD) {
      m->setIntTag(e, r->numberTag, &nf);
      r->toSplit[2][nf] = e;
      ++nf;
    }
  m->end(it);
  assert(static_cast<size_t>(nf) == r->toSplit[2].getSize());
  it = m->begin(3);
  int nr = 0;
  while ((e = m->iterate(it)))
    if ((m->getType(e) == apf::Mesh::PRISM)||
        (m->getType(e) == apf::Mesh::PYRAMID)) {
      m->setIntTag(e, r->numberTag, &nr);
      r->toSplit[3][nr] = e;
      ++nr;
    }
  m->end(it);
  assert(static_cast<size_t>(nr) == r->toSplit[3].getSize());
}

/* the commonly reused part of the
   tetrahedronization driver.
   as mentioned above, this uses refinement
   machinery, so these calls are copied
   from maRefine.cc */
void tetrahedronizeCommon(Refine* r)
{
  resetCollection(r);
  collectForTransfer(r);
  collectForMatching(r);
  splitElements(r);
  processNewElements(r);
  destroySplitElements(r);
  cleanupAfter(r);
}

/* the normal tetrahedronize-the-whole-layer function.
   this is where shouldTurnLayerToTets takes control
   and user reports are made */
void tetrahedronize(Adapt* a)
{
  if ( ! a->input->shouldTurnLayerToTets)
    return;
  assert(a->hasLayer);
  double t0 = PCU_Time();
  prepareLayerToTets(a);
  Refine* r = a->refine;
  addAllLayerElements(r);
  tetrahedronizeCommon(r);
  double t1 = PCU_Time();
  print("boundary layer converted to tets in %f seconds",t1-t0);
}

/* start of the island pyramid cleanup system,
   the goal of which is to remove pyramids which are
   adjacent to one another via their quadrilateral faces
   and which do not cap prismatic layer stacks.
   this allows us to later assume no such things
   exist, which simplifies the definition of LAYER.
*/

/* this Crawler sets the CHECKED flag on all layer quads
   which are reachable by crawling up from the base edges.
   this helps identify quads which are not reachable this way.
   recall from findLayerBase that LAYER_BASE is defined
   as edges on the closure of prisms, so the edges
   of quads between island pyramids will not be LAYER_BASE */
struct QuadMarker : public Crawler
{
  QuadMarker(Adapt* a_):
    Crawler(a_)
  {
    a = a_;
    m = a->mesh;
  }
  void begin(Layer& first)
  {
    getDimensionBase(a, 1, first);
    for (size_t i = 0; i < first.size(); ++i)
      setFlag(a, first[i], CHECKED);
  }
  void end()
  {
    clearFlagFromDimension(a, CHECKED, 1);
  }
  Entity* crawl(Entity* e)
  {
    HasFlag p(a, CHECKED);
    Entity* q = getOtherQuad(a, e, p);
    if (!q)
      return 0;
    setFlag(a, q, CHECKED);
    Entity* oe = getQuadEdgeOppositeEdge(m, q, e);
    setFlag(a, oe, CHECKED);
    return oe;
  }
  void send(Entity*, int)
  {
  }
  bool recv(Entity* e, int)
  {
    if (getFlag(a, e, CHECKED))
      return false;
    setFlag(a, e, CHECKED);
    return true;
  }
  Adapt* a;
  Mesh* m;
};

static void markNonIslandQuads(Adapt* a)
{
  QuadMarker op(a);
  crawlLayers(&op);
  syncFlag(a, 2, CHECKED);
}

/* quads which are not CHECKED at this
   point are between island pyramids,
   and are marked with SPLIT and DIAGONAL_1
   to prepare their tetrahedronization */
static void markIslandQuads(Adapt* a)
{
  Mesh* m = a->mesh;
  Entity* e;
  Iterator* it = m->begin(2);
  while ((e = m->iterate(it)))
    if (m->getType(e) == apf::Mesh::QUAD) {
      if ( ! getFlag(a, e, CHECKED)) {
        setFlag(a, e, SPLIT);
        setFlag(a, e, DIAGONAL_1);
      }
    }
  m->end(it);
  clearFlagFromDimension(a, CHECKED, 2);
  assert(checkFlagConsistency(a, 2, SPLIT));
  assert(checkFlagConsistency(a, 2, DIAGONAL_1));
}

/* this function marks pyramids adjacent
   to island quads with the SPLIT flag
   to set them up for tetrahedronization.
   it also checks that only pyramids
   are adjacent to the identified quads,
   verifying our assumption about how to
   identify island pyramids.
   the global total number of island
   pyramids is returned */
static long markIslandPyramids(Adapt* a)
{
  Mesh* m = a->mesh;
  Entity* e;
  long n = 0;
  Iterator* it = m->begin(2);
  while ((e = m->iterate(it)))
    if (getFlag(a, e, SPLIT)) {
      apf::Up up;
      m->getUp(e, up);
      for (int i = 0; i < up.n; ++i) {
        Entity* elem = up.e[i];
        assert(m->getType(elem) == apf::Mesh::PYRAMID);
        setFlag(a, elem, SPLIT);
        ++n;
      }
    }
  m->end(it);
  PCU_Add_Longs(&n, 1);
  return n;
}

static int countEntitiesWithFlag(Adapt* a, int flag, int dim)
{
  Mesh* m = a->mesh;
  Iterator* it = m->begin(dim);
  Entity* e;
  int n = 0;
  while ((e = m->iterate(it)))
    if (getFlag(a, e, flag))
      ++n;
  m->end(it);
  return n;
}

/* fills the refinement containers with
   the identified island quads and pyramids */
static void addIslandPyramids(Refine* r)
{
  Adapt* a = r->adapt;
  Mesh* m = a->mesh;
  for (int d = 2; d <= 3; ++d)
    r->toSplit[d].setSize(countEntitiesWithFlag(a, SPLIT, d));
  size_t n[4] = {};
  for (int d = 2; d <= 3; ++d) {
    Iterator* it = m->begin(d);
    Entity* e;
    while ((e = m->iterate(it)))
      if (getFlag(a, e, SPLIT))
        r->toSplit[d][n[d]++] = e;
    m->end(it);
    assert(r->toSplit[d].getSize() == n[d]);
  }
}

static long prepareIslandCleanup(Adapt* a)
{
  findLayerBase(a);
  markNonIslandQuads(a);
  markIslandQuads(a);
  return markIslandPyramids(a);
}

/* the main layer cleanup driver,
   currently responsible for removing
   island pyramids via tetrahedronization */
void cleanupLayer(Adapt* a)
{
  if (!a->hasLayer)
    return;
  if (!a->input->shouldCleanupLayer)
    return;
  double t0 = PCU_Time();
  long n = prepareIslandCleanup(a);
  if (!n) {
    print("no island pyramids found");
    return;
  }
  Refine* r = a->refine;
  addIslandPyramids(r);
  tetrahedronizeCommon(r);
  double t1 = PCU_Time();
  print("tetrahedronized %ld island pyramids in %f seconds", n, t1-t0);
}

}
