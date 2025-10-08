/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "maCollapse.h"
#include "maAdapt.h"
#include "maShape.h"
#include "maShapeHandler.h"
#include <apfCavityOp.h>
#include <pcu_util.h>

namespace ma {

void Collapse::Init(Adapt* a)
{
  adapt = a;
  cavity.init(a);
  rebuildCallback = 0;
  vertToCollapse = 0;
  vertToKeep = 0;
}

bool Collapse::requestLocality(apf::CavityOp* o)
{
/* get vertices again since this is sometimes used
   before setVerts can work */
  Entity* v[2];
  Mesh* m = adapt->mesh;
  m->getDownward(edge,0,v);
  return o->requestLocality(v,2);
}

// bool Collapse::collapseValid(double qualityToBeat)
// {
//   Vector prev = getPosition(adapt->mesh, vertToCollapse);
//   Vector target = getPosition(adapt->mesh, vertToKeep);
//   adapt->mesh->setPoint(vertToCollapse, 0, target);
//   bool output = true;
//   for (Entity* shape : newElements) {
//     // if (isLowInHigh(m, adjacent.e[i], vertToCollapse)) continue;
//     if (adapt->shape->getQuality(shape) < qualityToBeat) {
//       output = false;
//       break;
//     }
//   }
//   adapt->mesh->setPoint(vertToCollapse, 0, prev);
//   return output;
// }


bool Collapse::isValid(double qualityToBeat)
{
  Vector prev = getPosition(adapt->mesh, vertToCollapse);
  Vector target = getPosition(adapt->mesh, vertToKeep);
  adapt->mesh->setPoint(vertToCollapse, 0, target);
  EntityArray elements;
  elements.setSize(elementsToKeep.size());
  int i=0;
  for (Entity* e : elementsToKeep)
    elements[i++] = e;

  bool valid = areTetsValid(adapt->mesh, elements) 
            && !hasWorseQuality(adapt, elements, qualityToBeat);

  adapt->mesh->setPoint(vertToCollapse, 0, prev);
  return valid;
}

bool Collapse::tryThisDirectionNoCancel(double qualityToBeat)
{
  PCU_ALWAYS_ASSERT( ! adapt->mesh->isShared(vertToCollapse));
  rebuildElements();
  // check quality of linear t before fitting
  if ((adapt->mesh->getDimension()==2)
    &&( ! isGood2DMesh()))
    return false;
  // check the linear quality of tets before fitting
  // reject if one is negative
  if(adapt->mesh->getDimension()==3 && cavity.shouldFit
      && !areTetsValid(adapt->mesh,newElements))
    return false;
  // since they are okay in a linear sense, now fit and do a quality assessment
  fitElements();
  if (hasWorseQuality(adapt,newElements,qualityToBeat))
    return false;

  return true;
}

/*Sometimes the snapping procedure will attempt to collapse a edge that was just 
refined. This function will prevent that when the growth of the edge is above a 
thresh hold where it might not reach the target edge length*/
bool Collapse::edgesGoodSize() {
  PCU_ALWAYS_ASSERT(elementsToKeep.size());
  PCU_ALWAYS_ASSERT(elementsToKeep.size() == newElements.size());
  size_t i=0;
  double maxSize=0;
  double ratioAtMaxSize=0;
  APF_ITERATE(EntitySet,elementsToKeep,it) {
    Entity* edgesBefore[6];
    adapt->mesh->getDownward(*it,1,edgesBefore);
    Entity* edgesAfter[6];
    adapt->mesh->getDownward(newElements[i++],1,edgesAfter);
    for (int j=0; j<6; j++) {
      double sizeBefore = adapt->sizeField->measure(edgesBefore[j]);
      double sizeAfter = adapt->sizeField->measure(edgesAfter[j]);
      double ratio = sizeAfter/sizeBefore;
      if (sizeAfter > maxSize) {
        maxSize = sizeAfter;
        ratioAtMaxSize = ratio;
      }
    }
  }
  if (maxSize > MAXLENGTH && ratioAtMaxSize > MAXLENGTHRATIO) return false;
  return true;
}

bool Collapse::tryThisDirection(double qualityToBeat)
{
  if (!tryThisDirectionNoCancel(qualityToBeat)) {
    cancel();
    return false;
  }
  return true;
}

bool Collapse::tryBothDirections(double qualityToBeat)
{
  computeElementSets();
  // if adaptation is not forced, only accept if it
  // beats the smaller of good quality and the old quality
  // and is still of valid quality
  if (!adapt->input->shouldForceAdaptation)
    qualityToBeat = std::min(adapt->input->goodQuality,
        std::max(getOldQuality(),adapt->input->validQuality));

  if (tryThisDirectionNoCancel(qualityToBeat))
    return true;
  else destroyNewElements();

  if (getFlag(adapt,vertToKeep,COLLAPSE)) {
    std::swap(vertToKeep,vertToCollapse);
    computeElementSets();
    if (tryThisDirectionNoCancel(qualityToBeat))
      return true;
    else destroyNewElements();
  }

  unmark();
  return false;
}

bool Collapse::setEdge(Entity* e)
{
  if (getFlag(adapt,e,DONT_COLLAPSE))
    return false;
  edge = e;
  vertToCollapse = 0;
  vertToKeep = 0;
  elementsToCollapse.clear();
  elementsToKeep.clear();
  return true;
}

/* this routine ensures we have "mesh material" to pull from
   on at least one side of a collapsing entity.
   "a side (edge or tri) adjacent to a vertex to collapse must have the same classification
    as the center entity (tri or tet) being collapsed" */
static bool checkRingSide(Adapt* a, Entity* side, Entity* vert, Model* centerModel)
{
  Mesh* m = a->mesh;
  if (getFlag(a,vert,COLLAPSE)) {
    if (m->toModel(side)==centerModel)
      return true;
    else {
/* this function is responsible for disabling collapse of a
   vertex that has no material to pull from ! */
      clearFlagMatched(a,vert,COLLAPSE);
      return false;
    }
  }
  else
    return false;
}

/* if there exists a ring of three edges containing
   the collapsing edge, those edges must bound a triangle */
bool checkEdgeCollapseEdgeRings(Adapt* a, Entity* edge)
{
  Mesh* m = a->mesh;
  Entity* v[2];
  m->getDownward(edge,0,v);
  if (!getFlag(a, v[0], DONT_COLLAPSE)) //Allow collapse in one direction
    PCU_ALWAYS_ASSERT( ! m->isShared(v[0]));
  if (!getFlag(a, v[1], DONT_COLLAPSE)) //Allow collapse in one direction
    PCU_ALWAYS_ASSERT( ! m->isShared(v[1]));
  apf::Up ve[2];
  m->getUp(v[0],ve[0]);
  m->getUp(v[1],ve[1]);
  Entity* ring[3];
  ring[2] = edge;
  for (int i=0; i < ve[0].n; ++i)
  {
    ring[0] = ve[0].e[i];
    Entity* midv = getEdgeVertOppositeVert(m,ring[0],v[0]);
    for (int j=0; j < ve[1].n; ++j)
    {
      ring[1] = ve[1].e[j];
      if (midv != getEdgeVertOppositeVert(m,ring[1],v[1]))
        continue;
      //at this point a ring is found
      Entity* face = findUpward(m, apf::Mesh::TRIANGLE, ring);
      if (!face) return false; //the ring doesn't bound a face
      Model* c = m->toModel(face);
/* the face must have the same classification
   as one edge adjacent to a collapsing vertex.
   explicit booleans because checkRingSide has side
   effects and should be run on both vertices. */
      bool ok0 = checkRingSide(a,ring[0],v[0],c);
      bool ok1 = checkRingSide(a,ring[1],v[1],c);
      if ( ! (ok0 || ok1))
        return false;
    }
  }
  return true;
}

/* if there exist two edge-adjacent faces connecting the vertices,
   those faces and the edge must bound a tet */
bool checkEdgeCollapseFaceRings(Adapt* a, Entity* edge)
{
  Mesh* m = a->mesh;
  Entity* v[2];
  m->getDownward(edge,0,v);
  Upward vf[2];
  m->getAdjacent(v[0],2,vf[0]);
  m->getAdjacent(v[1],2,vf[1]);
/* this is here to speed up what would otherwise be an n^2 operation
   where n is the number of faces per vertex, n=~30 on average */
  std::map<Entity*,Entity*> oppositeEdgesToFaces;
  Entity* f[2];
  for (size_t i=0; i < vf[0].getSize(); ++i)
  {
    f[0] = vf[0][i];
/* we filter out non-tri faces; collapses involving the boundary
   layer should ensure topological correctness by other methods */
    if (( m->getType(f[0]) == apf::Mesh::TRIANGLE)
      &&( ! isInClosure(m,f[0],edge)))
    {
      Entity* oppositeEdge = getTriEdgeOppositeVert(m,f[0],v[0]);
      oppositeEdgesToFaces[oppositeEdge]=f[0];
    }
  }
  for (size_t i=0; i < vf[1].getSize(); ++i)
  {
    f[1] = vf[1][i];
    if ((m->getType(f[1]) != apf::Mesh::TRIANGLE)||
        (isInClosure(m, f[1], edge)))
      continue;
    Entity* oppositeEdge = getTriEdgeOppositeVert(m,f[1],v[1]);
    std::map<Entity*,Entity*>::iterator found =
      oppositeEdgesToFaces.find(oppositeEdge);
    if (found==oppositeEdgesToFaces.end()) continue;
    //at this point a ring is found
    f[0] = found->second;
    Entity* tet = findTetByTwoTris(m,f);
    if (!tet) return false;
    Model* c = m->toModel(tet);
/* at least one face adjacent to a collapsing vertex
   must have the same classification as the tet.
   explicit booleans because checkRingSide has side
   effects and should be run on both vertices. */
    bool ok0 = checkRingSide(a,f[0],v[0],c);
    bool ok1 = checkRingSide(a,f[1],v[1],c);
    if ( ! (ok0 || ok1))
      return false;
  }
  return true;
}

bool checkEdgeCollapseTopology(Adapt* a, Entity* edge)
{
  if ( ! checkEdgeCollapseEdgeRings(a,edge))
    return false;
  if ( ! checkEdgeCollapseFaceRings(a,edge))
    return false;
  return true;
}

static bool setVertexToCollapse(Adapt* a, Entity* v)
{
  if (getFlag(a,v,DONT_COLLAPSE))
    return false;
  setFlagMatched(a,v,COLLAPSE);
  return true;
}

/* this function checks the geometric classification
   requirements on an edge collapse
   and marks vertices with the COLLAPSE flag */
bool checkEdgeCollapseClassification(Adapt* a, Entity* edge)
{
  Mesh* m = a->mesh;
  Entity* v[2];
  m->getDownward(edge,0,v);
  Model* c[3];
  c[0] = m->toModel(v[0]);
  c[1] = m->toModel(edge);
  c[2] = m->toModel(v[1]);
  int cd[3];
  cd[0] = m->getModelType(c[0]);
  cd[1] = m->getModelType(c[1]);
  cd[2] = m->getModelType(c[2]);
  /* if the vertices are classified on equal-order model entities,
     then both vertices and the edge must have the same classification */
  if (cd[0] == cd[2])
  {
    if ( ! ((c[0]==c[1])&&(c[1]==c[2])))
      return false;
    bool ok[2];
    ok[0] = setVertexToCollapse(a,v[0]);
    ok[1] = setVertexToCollapse(a,v[1]);
    return ok[0] || ok[1];
  }
  //by now cd[2] != cd[0]
  if (cd[2] > cd[0])
  {
    std::swap(cd[0],cd[2]);
    std::swap(c[0],c[2]);
    std::swap(v[0],v[1]);
  }
  //now cd[0] > cd[2]
  /* if the vertices are classified on different orders, then the
     vertex to be collapsed and the edge must be classified on
     the higher-order model entity */
  if (c[0]!=c[1])
    return false;
  return setVertexToCollapse(a,v[0]);
}

bool Collapse::checkClass()
{
  if (checkEdgeCollapseClassification(adapt,edge))
    return true;
  clearFlagMatched(adapt,edge,COLLAPSE);
  return false;
}

bool Collapse::checkTopo()
{
  if (checkEdgeCollapseTopology(adapt,edge))
  {
    setVerts();
    return true;
  }
  unmark();
  return false;
}

void Collapse::unmark()
{ /* flags must be properly cleared to prevent misinterpretation
     later in the algorithms */
  clearFlagMatched(adapt,edge,COLLAPSE);
  Entity* v[2];
  Mesh* m = adapt->mesh;
  m->getDownward(edge,0,v);
  for (int i=0; i < 2; ++i)
    if ( ! isRequiredForAnEdgeCollapse(adapt,v[i]))
      clearFlagMatched(adapt,v[i],COLLAPSE);
}

void Collapse::setVerts()
{
  Entity* v[2];
  Mesh* m = adapt->mesh;
  m->getDownward(edge,0,v);
  if (getFlag(adapt,v[0],COLLAPSE))
    vertToCollapse = v[0];
  else
    vertToCollapse = v[1];
  PCU_ALWAYS_ASSERT(getFlag(adapt,vertToCollapse,COLLAPSE));
  vertToKeep = getEdgeVertOppositeVert(m,edge,vertToCollapse);
}

void Collapse::computeElementSets()
{
  Upward adjacent;
  Mesh* m = adapt->mesh;
  m->getAdjacent(edge,m->getDimension(),adjacent);
  elementsToCollapse.clear();
  APF_ITERATE(Upward,adjacent,it)
    elementsToCollapse.insert(*it);
  m->getAdjacent(vertToCollapse,m->getDimension(),adjacent);
  elementsToKeep.clear();
  APF_ITERATE(Upward,adjacent,it)
    if ( ! elementsToCollapse.count(*it))
      elementsToKeep.insert(*it);
  PCU_ALWAYS_ASSERT(elementsToKeep.size());
}

void Collapse::rebuildElements()
{
  PCU_ALWAYS_ASSERT(elementsToKeep.size());
  newElements.setSize(elementsToKeep.size());
  cavity.beforeBuilding();
  size_t ni=0;
  APF_ITERATE(EntitySet,elementsToKeep,it)
    newElements[ni++]=
        rebuildElement(adapt->mesh, *it, vertToCollapse, vertToKeep,
            adapt->buildCallback, rebuildCallback);
  cavity.afterBuilding();
}

void Collapse::fitElements()
{
  if(cavity.shouldFit){
    EntityArray oldElements;
    getOldElements(oldElements);
    cavity.fit(oldElements);
  }
}

bool Collapse::isGood2DMesh()
{
  /* in 2D we want to prevent "inverted" triangles
     we check that each rebuilt triangle has not had
     its normal changed by more than 90 degrees. */
  Mesh* m = adapt->mesh;
  size_t i=0;
  APF_ITERATE(EntitySet,elementsToKeep,it)
    if ( ! isTwoTriAngleAcute(m,*it,newElements[i++]))
      return false;
  return true;
}

void Collapse::destroyNewElements()
{
  for (size_t i=0; i < newElements.getSize(); ++i)
    destroyElement(adapt,newElements[i]);
  newElements.setSize(0);
}

void Collapse::cancel()
{
  destroyNewElements();
  unmark();
}

void Collapse::getOldElements(EntityArray& oldElements)
{
  EntitySet& toCollapse = elementsToCollapse;
  PCU_ALWAYS_ASSERT(toCollapse.size());
  EntitySet& toKeep = elementsToKeep;
  PCU_ALWAYS_ASSERT(toKeep.size());
  oldElements.setSize(toCollapse.size() + toKeep.size());
  size_t k=0;
  APF_ITERATE(EntitySet,toCollapse,it)
    oldElements[k++] = *it;
  APF_ITERATE(EntitySet,toKeep,it)
    oldElements[k++] = *it;
  PCU_ALWAYS_ASSERT(k==oldElements.getSize());
}

double Collapse::getOldQuality()
{
  EntityArray oldElements;
  getOldElements(oldElements);
  return getWorstQuality(adapt,oldElements);
}

void Collapse::destroyOldElements()
{
  EntityArray oldElements;
  getOldElements(oldElements);
  cavity.transfer(oldElements);
  for (size_t i=0; i < oldElements.getSize(); ++i)
    destroyElement(adapt,oldElements[i]);
}

bool isRequiredForAnEdgeCollapse(Adapt* adapt, Entity* vertex)
{
  Mesh* m = adapt->mesh;
  Model* c = m->toModel(vertex);
  apf::Up ve;
  m->getUp(vertex,ve);
  for (int i=0; i < ve.n; ++i)
  {
    Entity* edge = ve.e[i];
    /* we collapse one model dimension at a time, so
       only edges on the same model entity as the
       vertex being collapsed are relevant */
    if (c != m->toModel(edge))
      continue;
    if (getFlag(adapt,edge,COLLAPSE))
    {
      Entity* v2 = getEdgeVertOppositeVert(m,edge,vertex);
      if ( ! getFlag(adapt,v2,COLLAPSE))
        return true;
    }
  }
  return false;
}

bool isRequiredForMatchedEdgeCollapse(Adapt* adapt, Entity* vertex)
{
  if (adapt->mesh->hasMatching()) {
    apf::Matches matches;
    adapt->mesh->getMatches(vertex, matches);
    APF_ITERATE(apf::Matches, matches, it)
      if (isRequiredForAnEdgeCollapse(adapt, it->entity))
        return true;
  }
  return isRequiredForAnEdgeCollapse(adapt, vertex);
}

bool collapseEdgeVertex(Collapse& collapse, Entity* edge, Entity* vert)
{
  Adapt* adapt = collapse.adapt;
  double qualityToBeat = adapt->input->shouldForceAdaptation ? adapt->input->validQuality 
                                                  : adapt->input->goodQuality;
  if (!adapt->input->shouldForceAdaptation)
    qualityToBeat = std::min(adapt->input->goodQuality,
        std::max(collapse.getOldQuality(),adapt->input->validQuality));

  if (!setupCollapse(collapse, edge, vert))
    return false;
  if (!collapse.isValid(qualityToBeat)) {
    collapse.unmark();
    return false;
  }

  collapse.rebuildElements();
  collapse.fitElements();
  return true;
}

bool setupCollapse(Collapse& collapse, Entity* edge, Entity* vert)
{
  Adapt* adapter = collapse.adapt;
  PCU_ALWAYS_ASSERT(adapter->mesh->getType(edge) == apf::Mesh::EDGE);
  PCU_ALWAYS_ASSERT(adapter->mesh->getType(vert) == apf::Mesh::VERTEX);
  if ( ! collapse.setEdge(edge))
    return false;
  if ( ! collapse.checkClass())
    return false;
  if ( ! collapse.checkTopo())
    return false;
  if ( ! getFlag(adapter, vert, COLLAPSE)) {
    collapse.unmark();
    return false;
  }
  if (collapse.vertToCollapse != vert)
    std::swap(collapse.vertToCollapse, collapse.vertToKeep);
  PCU_ALWAYS_ASSERT(collapse.vertToCollapse == vert);
  collapse.computeElementSets();
  return true;
}

}
