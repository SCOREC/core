#include "maTemplates.h"
#include "maAdapt.h"
#include "maLayer.h"
#include "maSnap.h"

#include <cstdio>

namespace ma {

/* this is the template used on quads during tetrahedronization.
   the choice of diagonal has been made by the algorithms in
   maTetrahedronize.cc, we just retrieve it and split the quad */
void splitQuad_0(Refine* r, Entity* q, Entity** v)
{
  quadToTrisChoice(r, q, v, getDiagonalFromFlag(r->adapt, q));
}

/* splits the quad in half along edges v[0]-v[1] and v[2]-v[3] */
void splitQuad_2(Refine* r, Entity* q, Entity** v)
{
  Entity* sv[2];
  sv[0] = findSplitVert(r,v[0],v[1]);
  sv[1] = findSplitVert(r,v[2],v[3]);
  Entity* qv[4];
  qv[0] = v[0]; qv[1] = sv[0]; qv[2] = sv[1]; qv[3] = v[3];
  buildSplitElement(r,q,QUAD,qv);
  qv[0] = sv[0]; qv[1] = v[1]; qv[2] = v[2]; qv[3] = sv[1];
  buildSplitElement(r,q,QUAD,qv);
}

/* splits the quad into 4 along all edges */
void splitQuad_4(Refine* r, Entity* q, Entity** v)
{
  Adapt* a = r->adapt;
  Mesh* m = a->mesh;
  Entity* sv[4];
  double places[4];
  sv[0] = findPlacedSplitVert(r,v[0],v[1],places[0]);
  sv[1] = findPlacedSplitVert(r,v[1],v[2],places[1]);
  sv[2] = findPlacedSplitVert(r,v[3],v[2],places[2]);
  sv[3] = findPlacedSplitVert(r,v[0],v[3],places[3]);
  double x = (places[0] + places[2])/2;
  double y = (places[1] + places[3])/2;
  Vector xi(x*2-1,y*2-1,0);
/* since no rotation should have been applied, we actually don't
   need to unrotate xi */
  apf::MeshElement* me = apf::createMeshElement(m,q);
  Vector point;
  apf::mapLocalToGlobal(me,xi,point);
  Vector param(0,0,0); //prevents uninitialized values
  if (a->input->shouldTransferParametric)
    transferParametricOnQuadSplit(m, q, sv[0] ,sv[2], y, param);
  Entity* cv = buildVertex(a, m->toModel(q), point, param);
  a->solutionTransfer->onVertex(me,xi,cv);
  a->sizeField->interpolate(me,xi,cv);
  apf::destroyMeshElement(me);
  for (int i=0; i < 4; ++i)
  {
    Entity* v2[4];
    Entity* sv2[4];
    rotateQuad(v,i,v2);
    rotateQuad(sv,i,sv2);
    Entity* qv[4];
    qv[0] = v2[0]; qv[1] = sv2[0]; qv[2] = cv; qv[3] = sv2[3];
    buildSplitElement(r,q,QUAD,qv);
  }
}

SplitFunction quad_templates[quad_edge_code_count] = 
{splitQuad_0,
 splitQuad_2,
 splitQuad_4
};

/* this is the template used on prisms during tetrahedronization.
   the algorithms in maLayer.cc are required to ensure no centroids,
   so we check that and use the good case template */
void splitPrism_0(Refine* r, Entity* p, Entity** v)
{
  Adapt* a = r->adapt;
  Mesh* m = a->mesh;
  int code = getPrismDiagonalCode(m,v);
  assert(checkPrismDiagonalCode(code));
  prismToTetsGoodCase(r,p,v,code);
}

/* split the prism into two prisms separated by a quad face.
   edges v[0]-v[1] and v[3]-v[4] are split */
void splitPrism_2(Refine* r, Entity* p, Entity** v)
{
  Entity* sv[2];
  sv[0] = findSplitVert(r,v[0],v[1]);
  sv[1] = findSplitVert(r,v[3],v[4]);
  Entity* pv[6];
  pv[3] = v[3]; pv[4] = sv[1]; pv[5] = v[5];
  pv[0] = v[0]; pv[1] = sv[0]; pv[2] = v[2];
  buildSplitElement(r,p,PRISM,pv);
  pv[3] = sv[1]; pv[4] = v[4]; pv[5] = v[5];
  pv[0] = sv[0]; pv[1] = v[1]; pv[2] = v[2];
  buildSplitElement(r,p,PRISM,pv);
}

/* cuts a hex-shaped area in half along the
   v[0]-v[5] and v[3]-v[6] diagonals
   (a hex is ordered as top and bottom quads which curl up)
   resulting in two prisms */
static void hexToPrisms(Refine* r, Entity* p, Entity** v)
{
  Entity* pv[6];
  pv[3] = v[3]; pv[4] = v[6]; pv[5] = v[2];
  pv[0] = v[0]; pv[1] = v[5]; pv[2] = v[1];
  buildSplitElement(r,p,PRISM,pv);
  pv[3] = v[3]; pv[4] = v[7]; pv[5] = v[6];
  pv[0] = v[0]; pv[1] = v[4]; pv[2] = v[5];
  buildSplitElement(r,p,PRISM,pv);
}

static void simpleRotateHex(Entity** v, int n, Entity** v2)
{
  rotateQuad(v, n, v2); /* bottom quad */
  rotateQuad(v + 4, n, v2 + 4); /* top quad */
}

/* assumes either v[0]-v[5] or v[1]-v[4],
   rotates to the right one */
static void hexToPrisms2(Refine* r, Entity* p, Entity** v)
{
  int n = 0;
  if (!edgeExists(r->adapt->mesh, v[0], v[5]))
    n = 2;
  Entity* v2[8];
  simpleRotateHex(v, n, v2);
  hexToPrisms(r, p, v2);
}

/* extruded version of splitTri2 */
void splitPrism_4(Refine* r, Entity* p, Entity** v)
{
  Entity* sv[4];
  sv[0] = findSplitVert(r,v[0],v[1]);
  sv[1] = findSplitVert(r,v[3],v[4]);
  sv[2] = findSplitVert(r,v[4],v[5]);
  sv[3] = findSplitVert(r,v[1],v[2]);
  Entity* pv[6];
  pv[0] = sv[0]; pv[1] = v[1]; pv[2] = sv[3];
  pv[3] = sv[1]; pv[4] = v[4]; pv[5] = sv[2];
  buildSplitElement(r,p,PRISM,pv);
  Entity* hv[8];
  hv[0] = sv[1]; hv[1] = sv[2]; hv[2] = sv[3]; hv[3] = sv[0];
  hv[4] = v[3];  hv[5] =  v[5]; hv[6] =  v[2]; hv[7] =  v[0];
  hexToPrisms2(r,p,hv);
}

/* split the prism into four when all 6 triangular
   face edges have been split.
   This has to be told the split vertices since
   it is called from splitPrism_9 which has to deal
   with quad-associated vertices
 */
void splitPrism_6_sv(Refine* r, Entity* p, Entity** v, Entity** sv)
{
/* make center prism */
  buildSplitElement(r,p,PRISM,sv);
/* make the three corner prisms through rotation */
  for (int i=0; i < 3; ++i)
  {
    Entity* v2[6];
    Entity* sv2[6];
    rotatePrism(v,i,v2);
    rotatePrism(sv,i,sv2);
    Entity* pv[6];
    pv[3] = sv2[3]; pv[4] = v2[4]; pv[5] = sv2[4];
    pv[0] = sv2[0]; pv[1] = v2[1]; pv[2] = sv2[1];
    buildSplitElement(r,p,PRISM,pv);
  }
}

void splitPrism_6(Refine* r, Entity* p, Entity** v)
{
  Entity* sv[6];
  sv[0] = findSplitVert(r, v[0], v[1]);
  sv[1] = findSplitVert(r, v[1], v[2]);
  sv[2] = findSplitVert(r, v[2], v[0]);
  sv[3] = findSplitVert(r, v[3], v[4]);
  sv[4] = findSplitVert(r, v[4], v[5]);
  sv[5] = findSplitVert(r, v[5], v[3]);
  splitPrism_6_sv(r, p, v, sv);
}

/* uniform refinement of a prism: start by dividing
   the prism across the middle, splitting the vertical
   edges. give the two resulting sub-problems to 
   splitPrism_6_sv */
void splitPrism_9(Refine* r, Entity* p, Entity** v)
{
  Adapt* a = r->adapt;
  Mesh* m = a->mesh;
  Entity* botv[3];
  botv[0] = findSplitVert(r,v[0],v[1]);
  botv[1] = findSplitVert(r,v[1],v[2]);
  botv[2] = findSplitVert(r,v[2],v[0]);
  Entity* midv[3];
  midv[0] = findSplitVert(r,v[0],v[3]);
  midv[1] = findSplitVert(r,v[1],v[4]);
  midv[2] = findSplitVert(r,v[2],v[5]);
  Entity* topv[3];
  topv[0] = findSplitVert(r,v[3],v[4]);
  topv[1] = findSplitVert(r,v[4],v[5]);
  topv[2] = findSplitVert(r,v[5],v[3]);
  Entity* quads[3];
  Entity* qv[4];
  qv[3] = v[3]; qv[2] = v[4];
  qv[0] = v[0]; qv[1] = v[1];
  quads[0] = apf::findElement(m,QUAD,qv);
  qv[3] = v[4]; qv[2] = v[5];
  qv[0] = v[1]; qv[1] = v[2];
  quads[1] = apf::findElement(m,QUAD,qv);
  qv[3] = v[5]; qv[2] = v[3];
  qv[0] = v[2]; qv[1] = v[0];
  quads[2] = apf::findElement(m,QUAD,qv);
  Entity* cenv[3];
  cenv[0] = findSplitVert(r,quads[0]);
  cenv[1] = findSplitVert(r,quads[1]);
  cenv[2] = findSplitVert(r,quads[2]);
  Entity* pv[6];
  Entity* sv[6];
  pv[3] =    v[3]; pv[4] =    v[4]; pv[5] =    v[5];
  pv[0] = midv[0]; pv[1] = midv[1]; pv[2] = midv[2];
  sv[3] = topv[0]; sv[4] = topv[1]; sv[5] = topv[2];
  sv[0] = cenv[0]; sv[1] = cenv[1]; sv[2] = cenv[2];
  splitPrism_6_sv(r,p,pv,sv);
  pv[3] = midv[0]; pv[4] = midv[1]; pv[5] = midv[2];
  pv[0] =    v[0]; pv[1] =    v[1]; pv[2] =    v[2];
  sv[3] = cenv[0]; sv[4] = cenv[1]; sv[5] = cenv[2];
  sv[0] = botv[0]; sv[1] = botv[1]; sv[2] = botv[2];
  splitPrism_6_sv(r,p,pv,sv);
}

SplitFunction prism_templates[prism_edge_code_count] = 
{splitPrism_0
,splitPrism_2
,splitPrism_4
,splitPrism_6
,splitPrism_9
};

/* tetrahedronizes a tet-shaped sub-region given the edge split code */
static void splitSubTet(Refine* r, Entity* p, Entity** v, int code)
{
  if (!code) {
    buildSplitElement(r, p, TET, v);
    return;
  }
  Entity* v2[4];
  int tmpl_id = matchToTemplate(TET, v, code, v2);
  tet_templates[tmpl_id](r, p, v2);
}

/* any pyramid template that does not split the base will
   create a new pyramid using the base.
   Usually this sub-pyramid lines up along one edge of the parent.
   This function deals with the remaining sub-problem:
   a region consisting of the subtraction of the old and new pyramids.
   This is broken into two tets, and we recursively call
   tet templates depending on which of the remaining 3
   vertical edges are split */
static void splitPyramidTop(Refine* r, Entity* p, Entity** v)
{
  Adapt* a = r->adapt;
  Mesh* m = a->mesh;
  Entity* e[3];
  e[0] = findEdge(m, v[1], v[4]);
  e[1] = findEdge(m, v[2], v[4]);
  e[2] = findEdge(m, v[3], v[4]);
  int split[3];
  for (int i = 0; i < 3; ++i)
    split[i] = getFlag(a, e[i], SPLIT);
  Entity* tv[4];
  int code = (split[0] << 4) | (split[1] << 5);
  tv[0] = v[0]; tv[1] = v[1]; tv[2] = v[2]; tv[3] = v[4];
  splitSubTet(r, p, tv, code);
  code = (split[1] << 4) | (split[2] << 5);
  tv[0] = v[0]; tv[1] = v[2]; tv[2] = v[3]; tv[3] = v[4];
  splitSubTet(r, p, tv, code);
}

/* given a pyramid with no base edges split, tries
   to solve it as a sub-pyramid along a split edge v[0]-v[4],
   using splitPyramidTop */
static void splitPyramidCommon(Refine* r, Entity* p, Entity** v)
{
  Entity* sv = findSplitVert(r, v[0], v[4]);
  Entity* pv[5];
  pv[0] = v[0]; pv[1] = v[1]; pv[2] = v[2]; pv[3] = v[3];
  pv[4] = sv;
  buildSplitElement(r, p, PYRAMID, pv);
  pv[0] = sv;
  pv[4] = v[4];
  splitPyramidTop(r, p, pv);
}

/* checks whether the splitting thats happened around
   this pyramid allows a sub-pyramid to exist along v[0]-v[4] */
static bool canSplitPyramidCommon(Refine* r, Entity** v)
{
  Adapt* a = r->adapt;
  Mesh* m = a->mesh;
  Entity* e = findEdge(m, v[0], v[4]);
  if (!getFlag(a, e, SPLIT))
    return false;
  Entity* sv = findSplitVert(r, e);
  return edgeExists(m, v[1], sv) && edgeExists(m, v[3], sv);
}

/* given a pyramid with no base edges split, tries
   to solve it as a sub-pyramid along a split edge v[0]-v[4],
   using splitPyramidTop */
static bool splitPyramidSearch(Refine* r, Entity* p, Entity** v)
{
  for (int i = 0; i < 4; ++i) {
    Entity* pv[5];
    rotatePyramid(v, i, pv);
    if (canSplitPyramidCommon(r, pv)) {
      splitPyramidCommon(r, p, pv);
      return true;
    }
  }
  return false;
}

/* vertical edge v[0]-v[4] split */
static void splitPyramid_1_b0(Refine* r, Entity* p, Entity** v)
{
  splitPyramidCommon(r, p, v);
}

/* vertical edges v[0]-v[4] and v[1]-v[4] split */
static void splitPyramid_2_b0(Refine* r, Entity* p, Entity** v)
{
  /* ambiguous diagonal on face 0,1,4 requires searching */
  bool ok = splitPyramidSearch(r, p, v);
  assert(ok);
}

/* vertical edges v[0]-v[4] and v[2]-v[4] split */
static void splitPyramid_3_b0(Refine* r, Entity* p, Entity** v)
{
  splitPyramidCommon(r, p, v);
}

/* vertical edges v[1]-v[4], v[2]-v[4] and v[3]-v[4] split */
static void splitPyramid_4_b0(Refine* r, Entity* p, Entity** v)
{
  /* two ambiguous diagonals here, but all cases should
     reduce to the common case */
  bool ok = splitPyramidSearch(r, p, v);
  assert(ok);
}

static Entity* makePyramidCentroid(Adapt* a, Entity* p)
{
  Mesh* m = a->mesh;
  Vector param(0,0,0); //will be in geometric region
  Model* c = m->toModel(p);
  apf::MeshElement* me = apf::createMeshElement(m, p);
  Vector xi(0,0,-3./5.); //parametric centroid for degenerate hex pyramid
  Vector point;
  apf::mapLocalToGlobal(me, xi, point);
  Entity* v = buildVertex(a, c, point, param);
  a->solutionTransfer->onVertex(me, xi, v);
  a->sizeField->interpolate(me, xi, v);
  apf::destroyMeshElement(me);
  return v;
}

/* all four vertical edges are split and all ambiguous
   diagonals curl the same way. This is the bad case for the
   pyramid, in which we must introduce a central vertex */
static void splitPyramid_5_b0_bad(Refine* r, Entity* p, Entity** v)
{
  Entity* pv[5];
  Entity* cv = makePyramidCentroid(r->adapt, p);
  for (int i = 0; i < 4; ++i)
    pv[i] = v[i];
  pv[4] = cv;
  buildSplitElement(r, p, PYRAMID, pv);
  Entity* sv[4];
  for (int i = 0; i < 4; ++i)
    sv[i] = findSplitVert(r, v[i], v[4]);
  Entity* v2[5];
  Entity* sv2[4];
  for (int i = 0; i < 4; ++i) {
    rotatePyramid(v, i, v2);
    rotateQuad(sv, i, sv2);
    pv[0] = v2[0]; pv[1] = sv2[0]; pv[2] = sv2[1]; pv[3] = v2[1];
    pyramidToTets(r, p, pv);
  }
  Entity* ov[6];
  ov[0] = cv;
  for (int i = 0; i < 4; ++i)
    ov[i + 1] = sv[i];
  ov[5] = v[4];
  octToTetsGeometric(r, p, ov);
}

/* all four vertical edges are split */
static void splitPyramid_5_b0(Refine* r, Entity* p, Entity** v)
{
  /* there are four ambiguous diagonals and thus 2^4 = 16
     combinations. two of those are bad */
  bool good = splitPyramidSearch(r, p, v);
  if (!good)
    splitPyramid_5_b0_bad(r, p, v);
}

/* split the pyramid into two new ones when edges
   v[0]-v[1] and v[2]-v[3] have been split */
static void splitPyramid_b1_gen(Refine* r, Entity* p, Entity** v,
    SplitFunction subproblems[2])
{
  Entity* sv[2];
  sv[0] = findSplitVert(r,v[0],v[1]);
  sv[1] = findSplitVert(r,v[2],v[3]);
  Entity* pv[5];
  pv[0] = v[0]; pv[1] = sv[0];
  pv[3] = v[3]; pv[2] = sv[1];
  pv[4] = v[4];
  subproblems[0](r, p, pv);
  pv[0] = sv[0]; pv[1] = v[1];
  pv[3] = sv[1]; pv[2] = v[2];
  pv[4] = v[4];
  subproblems[1](r, p, pv);
}

static void buildPyramid(Refine* r, Entity* p, Entity** v)
{
  buildSplitElement(r, p, PYRAMID, v);
}

/* b1 means the base quad is split along v[0]-v[1] and v[2]-v[3].
   refer to the b0 templates for the state of the vertical
   edges. All cases where the base quad is split in half are
   handled by cutting the pyramid in half and handing those
   subproblems to b0 templates */
static void splitPyramid_0_b1(Refine* r, Entity* p, Entity** v)
{
  SplitFunction subproblems[2] = {buildPyramid, buildPyramid};
  splitPyramid_b1_gen(r, p, v, subproblems);
}

/* case 1_b1, ambiguous diagonal goes from base to cap, allowing
   a separation of subproblems */
static void splitPyramid_1_b1_a(Refine* r, Entity* p, Entity** v)
{
  SplitFunction subproblems[2] = {splitPyramid_1_b0, buildPyramid};
  splitPyramid_b1_gen(r, p, v, subproblems);
}

/* case 1_b1, ambiguous diagonal from center of split vertical */
static void splitPyramid_1_b1_b(Refine* r, Entity* p, Entity** v)
{
  Entity* sv[3];
  sv[0] = findSplitVert(r, v[0], v[1]);
  sv[1] = findSplitVert(r, v[2], v[3]);
  sv[2] = findSplitVert(r, v[0], v[4]);
  Entity* pv[5];
  pv[0] = v[0]; pv[1] = sv[0]; pv[2] = sv[1]; pv[3] = v[3];
  pv[4] = sv[2];
  buildSplitElement(r, p, PYRAMID, pv);
  pv[0] = sv[0]; pv[1] = v[1]; pv[2] = v[2]; pv[3] = sv[1];
  buildPyramid(r, p, pv);
  pv[0] = sv[2]; pv[1] = v[1]; pv[2] = v[2]; pv[3] = sv[1];
  pv[4] = v[4];
  splitPyramidTop(r, p, pv);
  Entity* tv[4];
  tv[0] = v[3]; tv[1] = sv[2]; tv[2] = sv[1]; tv[3] = v[4];
  buildSplitElement(r, p, TET, tv);
}

static void splitPyramid_1_b1(Refine* r, Entity* p, Entity** v)
{
  Entity* ev[2];
  ev[0] = findSplitVert(r, v[0], v[1]);
  ev[1] = v[4];
  if (findUpward(r->adapt->mesh, EDGE, ev))
    splitPyramid_1_b1_a(r, p, v);
  else
    splitPyramid_1_b1_b(r, p, v);
}

static void splitPyramid_2_b1(Refine* r, Entity* p, Entity** v)
{
  Entity* sv[4];
  sv[0] = findSplitVert(r, v[0], v[1]);
  sv[1] = findSplitVert(r, v[2], v[3]);
  sv[2] = findSplitVert(r, v[0], v[4]);
  sv[3] = findSplitVert(r, v[1], v[4]);
  Entity* pv[5];
  pv[0] = v[0]; pv[1] = sv[0]; pv[2] = sv[1]; pv[3] = v[3];
  pv[4] = sv[2];
  buildSplitElement(r, p, PYRAMID, pv);
  pv[0] = sv[0]; pv[1] = v[1]; pv[2] = v[2]; pv[3] = sv[1];
  pv[4] = sv[3];
  buildSplitElement(r, p, PYRAMID, pv);
  Entity* tv[4];
  tv[0] = sv[1]; tv[1] = sv[2]; tv[2] = sv[3];
  tv[3] = v[4];
  buildSplitElement(r, p, TET, tv);
  tv[0] = v[3]; tv[1] = sv[2]; tv[2] = sv[1];
  buildSplitElement(r, p, TET, tv);
  tv[0] = sv[1]; tv[1] = sv[3]; tv[2] = v[2];
  buildSplitElement(r, p, TET, tv);
}

static void splitPyramid_3_b1(Refine* r, Entity* p, Entity** v)
{
  Entity* sv[4];
  sv[0] = findSplitVert(r, v[0], v[1]);
  sv[1] = findSplitVert(r, v[2], v[3]);
  sv[2] = findSplitVert(r, v[0], v[4]);
  sv[3] = findSplitVert(r, v[2], v[4]);
  Entity* pv[5];
  pv[0] = v[0]; pv[1] = sv[0]; pv[2] = sv[1]; pv[3] = v[3];
  pv[4] = sv[2];
  buildSplitElement(r, p, PYRAMID, pv);
  pv[0] = sv[0]; pv[1] = v[1]; pv[2] = v[2]; pv[3] = sv[1];
  pv[4] = sv[3];
  buildSplitElement(r, p, PYRAMID, pv);
  buildSplitElement(r, p, TET, sv);
  pv[0] = sv[0]; pv[1] = sv[2]; pv[2] = v[4]; pv[3] = v[1];
  pv[4] = sv[3];
  pyramidToTets(r, p, pv);
  pv[0] = sv[1]; pv[1] = sv[3]; pv[2] = v[4]; pv[3] = v[3];
  pv[4] = sv[2];
  pyramidToTets(r, p, pv);
}

/* splits a prism area into a tet and a pyramid, the pyramid
   base is the prism's first face and its cap is v[2] */
static void prismToPyramidAndTet(Refine* r, Entity* p, Entity** v)
{
  Entity* pv[5] = {v[0],v[3],v[4],v[1],v[2]};
  buildSplitElement(r, p, PYRAMID, pv);
  Entity* tv[4] = {v[2],v[3],v[4],v[5]};
  buildSplitElement(r, p, TET, tv);
}

/* invokes the above either as-is or upside-down,
   taking its cues from one or more edges which
   may already exist along the side of the prism */
static void smartPrismToPyramidAndTet(Refine* r, Entity* p, Entity** v)
{
  Mesh* m = r->adapt->mesh;
  bool isUpsideDown = edgeExists(m, v[0], v[5]) || edgeExists(m, v[1], v[5]);
  int rotation = isUpsideDown ? 4 : 0;
  Entity* v2[6];
  rotatePrism(v, rotation, v2);
  prismToPyramidAndTet(r, p, v2);
}

static void splitPyramid_4_b1(Refine* r, Entity* p, Entity** v)
{
  Entity* sv[5];
  sv[0] = findSplitVert(r, v[0], v[1]);
  sv[1] = findSplitVert(r, v[2], v[3]);
  sv[2] = findSplitVert(r, v[1], v[4]);
  sv[3] = findSplitVert(r, v[2], v[4]);
  sv[4] = findSplitVert(r, v[3], v[4]);
  Entity* wv[6] = {v[1],sv[0],sv[2],v[2],sv[1],sv[3]};
  smartPrismToPyramidAndTet(r, p, wv);
  Entity* pv[5];
  pv[0] = sv[0]; pv[1] = sv[2]; pv[2] = sv[3]; pv[3] = sv[1];
  pv[4] = sv[4];
  pyramidToTets(r, p, pv);
  pv[0] = v[0]; pv[1] = sv[0]; pv[2] = sv[2]; pv[3] = v[4];
  pyramidToTets(r, p, pv);
  pv[0] = v[3]; pv[1] = v[0]; pv[2] = sv[0]; pv[3] = sv[1];
  buildSplitElement(r, p, PYRAMID, pv);
  Entity* tv[4] = {sv[2],sv[3],sv[4],v[4]};
  buildSplitElement(r, p, TET, tv);
}

static void splitPyramid_5_b1(Refine* r, Entity* p, Entity** v)
{
  Entity* sv[6];
  sv[0] = findSplitVert(r, v[0], v[1]);
  sv[1] = findSplitVert(r, v[2], v[3]);
  for (int i = 0; i < 4; ++i)
    sv[2 + i] = findSplitVert(r, v[i], v[4]);
  Entity* wv[6];
  wv[0] = sv[0]; wv[1] = v[0]; wv[2] = sv[2];
  wv[3] = sv[1]; wv[4] = v[3]; wv[5] = sv[5];
  smartPrismToPyramidAndTet(r, p, wv);
  wv[0] = v[1]; wv[1] = sv[0]; wv[2] = sv[3];
  wv[3] = v[2]; wv[4] = sv[1]; wv[5] = sv[4];
  smartPrismToPyramidAndTet(r, p, wv);
  wv[0] = sv[2]; wv[1] = sv[3]; wv[2] = sv[0];
  wv[3] = sv[5]; wv[4] = sv[4]; wv[5] = sv[1];
  prismAndPyramidToTets(r, p, wv, v[4]);
}

static void splitPyramid_1_b2_a(Refine* r, Entity* p, Entity** v)
{
  Entity* sv[2];
  sv[0] = findSplitVert(r, v[1], v[2]);
  sv[1] = findSplitVert(r, v[3], v[0]);
  Entity* pv[5];
  pv[4] = v[4];
  pv[0] = v[0]; pv[1] = v[1]; pv[2] = sv[0]; pv[3] = sv[1];
  splitPyramid_1_b0(r, p, pv);
  pv[0] = sv[0]; pv[1] = v[2]; pv[2] = v[3]; pv[4] = sv[1];
  buildSplitElement(r, p, PYRAMID, pv);
}

static void splitPyramid_1_b2_b(Refine* r, Entity* p, Entity** v)
{
  Entity* sv[3];
  sv[0] = findSplitVert(r, v[1], v[2]);
  sv[1] = findSplitVert(r, v[3], v[0]);
  sv[2] = findSplitVert(r, v[0], v[4]);
  Entity* pv[5];
  pv[0] = v[0]; pv[1] = v[1]; pv[2] = sv[0]; pv[3] = sv[1];
  pv[4] = sv[2];
  buildSplitElement(r, p, PYRAMID, pv);
  pv[0] = sv[1]; pv[1] = sv[0]; pv[2] = v[2]; pv[3] = v[3];
  buildSplitElement(r, p, PYRAMID, pv);
  pv[0] = sv[2]; pv[1] = sv[0]; pv[2] = v[2]; pv[3] = v[3];
  pv[4] = v[4];
  splitPyramidTop(r, p, pv);
  Entity* tv[4] = {sv[2],v[1],sv[0],v[4]};
  buildSplitElement(r, p, TET, tv);
}

static void splitPyramid_1_b2(Refine* r, Entity* p, Entity** v)
{
  Entity* sv = findSplitVert(r, v[0], v[3]);
  if (edgeExists(r->adapt->mesh, sv, v[4]))
    splitPyramid_1_b2_a(r, p, v);
  else
    splitPyramid_1_b2_b(r, p, v);
}

typedef void (*SplitPyramid_2_b2_sub)(Refine* r, Entity* p, Entity** v,
    Entity** sv);

void splitPyramid_2_b2_sub_0(Refine* r, Entity* p, Entity** v,
    Entity** sv)
{
  Entity* wv[6] = {v[2], sv[0], sv[3], v[3], sv[1], sv[2]};
  smartPrismToPyramidAndTet(r, p, wv);
  Entity* pv[5] = {v[3], sv[2], sv[3], v[2], v[4]};
  pyramidToTets(r, p, pv);
}

void splitPyramid_2_b2_bad(Refine* r, Entity* p, Entity** v,
    Entity** sv)
{
  Entity* pv[5];
  Entity* cv = makePyramidCentroid(r->adapt, p);
  pv[4] = cv;
  pv[0] = sv[0]; pv[1] = v[2]; pv[2] = v[3]; pv[3] = sv[1];
  buildSplitElement(r, p, PYRAMID, pv);
  for (int i = 0; i < 4; ++i)
    pv[i] = sv[i];
  pyramidToTets(r, p, pv);
  pv[0] = sv[1]; pv[1] = v[3]; pv[2] = v[4]; pv[3] = sv[2];
  pyramidToTets(r, p, pv);
  pv[0] = sv[0]; pv[1] = sv[3]; pv[2] = v[4]; pv[3] = v[2];
  pyramidToTets(r, p, pv);
  Entity* tv[4] = {v[2],v[4],v[3],cv};
  buildSplitElement(r, p, TET, tv);
}

void splitPyramid_2_b2_sub_1(Refine* r, Entity* p, Entity** v,
    Entity** sv)
{
  Entity* pv[5];
  pv[0] = sv[0]; pv[1] = v[2]; pv[2] = v[3]; pv[3] = sv[1];
  pv[4] = sv[2];
  buildSplitElement(r, p, PYRAMID, pv);
  pv[0] = sv[0]; pv[1] = sv[3]; pv[2] = v[4]; pv[3] = v[2];
  pyramidToTets(r, p, pv);
  Entity* tv[4] = {v[2],v[4],v[3],sv[2]};
  buildSplitElement(r, p, TET, tv);
}

void splitPyramid_2_b2_sub_2(Refine* r, Entity* p, Entity** v,
    Entity** sv)
{
  Entity* pv[5];
  pv[0] = sv[0]; pv[1] = v[2]; pv[2] = v[3]; pv[3] = sv[1];
  pv[4] = sv[3];
  buildSplitElement(r, p, PYRAMID, pv);
  pv[0] = sv[1]; pv[1] = v[3]; pv[2] = v[4]; pv[3] = sv[2];
  pyramidToTets(r, p, pv);
  Entity* tv[4] = {v[2],v[4],v[3],sv[3]};
  buildSplitElement(r, p, TET, tv);
}

void splitPyramid_2_b2_sub_3(Refine* r, Entity* p, Entity** v,
    Entity** sv)
{
  Entity* pv[5];
  pv[0] = sv[0]; pv[1] = v[2]; pv[2] = v[3]; pv[3] = sv[1];
  pv[4] = v[4];
  buildSplitElement(r, p, PYRAMID, pv);
  pv[0] = sv[0]; pv[1] = sv[1]; pv[2] = sv[2]; pv[3] = sv[3];
  pyramidToTets(r, p, pv);
}

static void splitPyramid_2_b2(Refine* r, Entity* p, Entity** v)
{
  Mesh* m = r->adapt->mesh;
  Entity* sv[4];
  sv[0] = findSplitVert(r, v[1], v[2]);
  sv[1] = findSplitVert(r, v[3], v[0]);
  sv[2] = findSplitVert(r, v[0], v[4]);
  sv[3] = findSplitVert(r, v[1], v[4]);
  Entity* pv[6];
  pv[0] = v[0]; pv[1] = sv[1]; pv[2] = sv[2];
  pv[3] = v[1]; pv[4] = sv[0]; pv[5] = sv[3];
  smartPrismToPyramidAndTet(r, p, pv);
  int diagonals = 0;
  if (edgeExists(m, sv[0], sv[2]))
    diagonals |= (1<<0);
  if (edgeExists(m, sv[0], v[4]))
    diagonals |= (1<<1);
  if (edgeExists(m, sv[1], v[4]))
    diagonals |= (1<<2);
  static SplitPyramid_2_b2_sub const table[1<<3] =
  {splitPyramid_2_b2_sub_0//000 - sub_0 case 0 (lower prism)
  ,splitPyramid_2_b2_sub_0//001 - sub_0 case 1 (lower prism)
  ,splitPyramid_2_b2_bad  //010 - bad case 0 (centroid)
  ,splitPyramid_2_b2_sub_1//011 - sub_1 (skew pyramid)
  ,splitPyramid_2_b2_sub_2//100 - sub_2 (skew pyramid mirror)
  ,splitPyramid_2_b2_bad  //101 - bad case 1 (centroid other curl)
  ,splitPyramid_2_b2_sub_3//110 - sub_3 (symmetric pyramid)
  ,splitPyramid_2_b2_sub_3//111 - sub_3 (symmetric pyramid)
  };
  table[diagonals](r, p, v, sv);
}

static void splitPyramid_3_b2(Refine* r, Entity* p, Entity** v)
{/* mirror of 3_b1 */
  Entity* sv[4];
  sv[0] = findSplitVert(r, v[1], v[2]);
  sv[1] = findSplitVert(r, v[3], v[0]);
  sv[2] = findSplitVert(r, v[0], v[4]);
  sv[3] = findSplitVert(r, v[2], v[4]);
  Entity* pv[5];
  pv[0] = v[1]; pv[1] = sv[0]; pv[2] = sv[1]; pv[3] = v[0];
  pv[4] = sv[2];
  buildSplitElement(r, p, PYRAMID, pv);
  pv[0] = sv[0]; pv[1] = v[2]; pv[2] = v[3]; pv[3] = sv[1];
  pv[4] = sv[3];
  buildSplitElement(r, p, PYRAMID, pv);
  buildSplitElement(r, p, TET, sv);
  pv[0] = sv[0]; pv[1] = v[1]; pv[2] = v[4]; pv[3] = sv[3];
  pv[4] = sv[2];
  pyramidToTets(r, p, pv);
  pv[0] = sv[1]; pv[1] = v[3]; pv[2] = v[4]; pv[3] = sv[2];
  pv[4] = sv[3];
  pyramidToTets(r, p, pv);
}

static void splitPyramid_4_b2(Refine* r, Entity* p, Entity** v)
{ /* mirror of 4_b1 */
  Entity* sv[5];
  sv[0] = findSplitVert(r, v[1], v[2]);
  sv[1] = findSplitVert(r, v[3], v[0]);
  sv[2] = findSplitVert(r, v[1], v[4]);
  sv[3] = findSplitVert(r, v[2], v[4]);
  sv[4] = findSplitVert(r, v[3], v[4]);
  Entity* wv[6] = {v[2],sv[0],sv[3],v[3],sv[1],sv[4]};
  smartPrismToPyramidAndTet(r, p, wv);
  Entity* pv[5];
  pv[0] = sv[0]; pv[1] = sv[3]; pv[2] = sv[4]; pv[3] = sv[1];
  pv[4] = sv[2];
  pyramidToTets(r, p, pv);
  pv[0] = v[0]; pv[1] = sv[1]; pv[2] = sv[4]; pv[3] = v[4];
  pyramidToTets(r, p, pv);
  pv[0] = v[0]; pv[1] = v[1]; pv[2] = sv[0]; pv[3] = sv[1];
  buildSplitElement(r, p, PYRAMID, pv);
  Entity* tv[4] = {sv[2],sv[3],sv[4],v[4]};
  buildSplitElement(r, p, TET, tv);
}

/* uniform refinement of a pyramid */
static void splitPyramidUniform(Refine* r, Entity* p, Entity** v)
{
  Mesh* m = r->adapt->mesh;
  Entity* botv[4];
  botv[0] = findSplitVert(r,v[0],v[1]);
  botv[1] = findSplitVert(r,v[1],v[2]);
  botv[2] = findSplitVert(r,v[2],v[3]);
  botv[3] = findSplitVert(r,v[3],v[0]);
  Entity* midv[4];
  midv[0] = findSplitVert(r,v[0],v[4]);
  midv[1] = findSplitVert(r,v[1],v[4]);
  midv[2] = findSplitVert(r,v[2],v[4]);
  midv[3] = findSplitVert(r,v[3],v[4]);
/* the first four entries of v also specify the bottom quad */
  Entity* quad = apf::findElement(m,QUAD,v);
  Entity* cv = findSplitVert(r,quad);
/* make the four new pyramids and four new tets by rotation */
  for (int i=0; i < 4; ++i)
  {
    Entity* midv2[4];
    rotateQuad(midv,i,midv2);
    Entity* botv2[4];
    rotateQuad(botv,i,botv2);
    Entity* v2[4];
    rotateQuad(v,i,v2);
    Entity* pv[5];
    pv[3] = cv;       pv[2] = botv2[1];
    pv[0] = botv2[0]; pv[1] = v2[1];
    pv[4] = midv2[1];
    buildSplitElement(r,p,PYRAMID,pv);
    Entity* tv[4];
    tv[0] = midv2[0]; tv[1] = midv2[1];
    tv[2] = botv2[0]; tv[3] = cv;
    buildSplitElement(r,p,TET,tv);
  }
/* tetrahedronize the central octahedron by choosing the
   shortest diagonal */
  Entity* octv[6];
  octv[5] = v[4];
  octv[4] = midv[3]; octv[3] = midv[2];
  octv[1] = midv[0]; octv[2] = midv[1];
  octv[0] = cv;
  octToTetsGeometric(r,p,octv);
}

SplitFunction pyramid_templates[pyramid_edge_code_count] = 
{pyramidToTets
,splitPyramid_1_b0
,splitPyramid_2_b0
,splitPyramid_3_b0
,splitPyramid_4_b0
,splitPyramid_5_b0
,splitPyramid_0_b1
,splitPyramid_1_b1
,splitPyramid_2_b1
,splitPyramid_3_b1
,splitPyramid_4_b1
,splitPyramid_5_b1
,splitPyramid_1_b2
,splitPyramid_2_b2
,splitPyramid_3_b2
,splitPyramid_4_b2
,splitPyramidUniform
};

}
