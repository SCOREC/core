#include "maTemplates.h"
#include "maAdapt.h"
#include "maLayer.h"

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
/* TODO: in truth, we should be transferring parametric
   coordinates here. since we are a ways away from boundary
   layer snapping and the transfer logic is non-trivial,
   this is left alone for now */
  Entity* cv = buildVertex(a,m->toModel(q),point,Vector(0,0,0));
  a->sizeField->interpolate(me,xi,cv);
  a->solutionTransfer->onVertex(me,xi,cv);
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
  Mesh* m = r->adapt->mesh;
  Entity* ev[2];
  ev[0] = v[0]; ev[1] = v[5];
  int n = 0;
  if (!findUpward(m, EDGE, ev))
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

/* split the pyramid into two new ones when edges
   v[0]-v[1] and v[2]-v[3] have been split */
void splitPyramid_2(Refine* r, Entity* p, Entity** v)
{
  Entity* sv[2];
  sv[0] = findSplitVert(r,v[0],v[1]);
  sv[1] = findSplitVert(r,v[2],v[3]);
  Entity* pv[5];
  pv[0] = v[0]; pv[1] = sv[0];
  pv[3] = v[3]; pv[2] = sv[1];
  pv[4] = v[4];
  buildSplitElement(r,p,PYRAMID,pv);
  pv[0] = sv[0]; pv[1] = v[1];
  pv[3] = sv[1]; pv[2] = v[2];
  pv[4] = v[4];
  buildSplitElement(r,p,PYRAMID,pv);
}

/* uniform refinement of a pyramid */
void splitPyramid_4(Refine* r, Entity* p, Entity** v)
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
   best metric space diagonal */
  Entity* octv[6];
  octv[5] = v[4];
  octv[4] = midv[3]; octv[3] = midv[2];
  octv[1] = midv[0]; octv[2] = midv[1];
  octv[0] = cv;
  octToTetsGeometric(r,p,octv);
}

SplitFunction pyramid_templates[pyramid_edge_code_count] = 
{pyramidToTets
,splitPyramid_2
,splitPyramid_4
};

}
