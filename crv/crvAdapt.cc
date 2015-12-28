/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "crv.h"
#include "crvAdapt.h"
#include <apf.h>
#include <apfMesh.h>
#include <maBalance.h>
#include <maCoarsen.h>
#include <maShape.h>
#include <maSnap.h>
#include <PCU.h>
#include <cassert>

namespace crv {

Adapt::Adapt(ma::Input* in)
: ma::Adapt(in)
{
  validityTag = mesh->createIntTag("crv_flags",1);
}

static void clearFlags(Adapt* a)
{
  ma::Mesh* m = a->mesh;
  ma::Entity* e;
  for (int d=0; d <= 3; ++d)
  {
    ma::Iterator* it = m->begin(d);
    while ((e = m->iterate(it)))
      if (m->hasTag(e,a->validityTag))
        m->removeTag(e,a->validityTag);
    m->end(it);
  }
  m->destroyTag(a->validityTag);
}

static int getFlags(Adapt* a, ma::Entity* e)
{
  ma::Mesh* m = a->mesh;
  if ( ! m->hasTag(e,a->validityTag))
    return 0; //we assume 0 is the default value for all flags
  int flags;
  m->getIntTag(e,a->validityTag,&flags);
  return flags;
}

static void setFlags(Adapt* a, ma::Entity* e, int flags)
{
  a->mesh->setIntTag(e,a->validityTag,&flags);
}

void snapRefineToBoundary(ma::Adapt* a){

  if ( ! a->input->shouldSnap)
    return;

  ma::Mesh* m = a->mesh;
  ma::Refine* r = a->refine;
  int P = m->getShape()->getOrder();
  apf::FieldShape* fs = m->getShape();
  for (int d=1; d <= m->getDimension(); ++d){
    for (size_t i=0; i < r->newEntities[d].getSize(); ++i){
      ma::EntityArray& a = r->newEntities[d][i];
      for (size_t i=0; i < a.getSize(); ++i)

        if(m->getModelType(m->toModel(a[i])) < m->getDimension()){
          snapToInterpolate(m,a[i]);
        }
    }
  }
  for (int d = m->getDimension(); d >=1; --d){
    for (size_t i=0; i < r->newEntities[d].getSize(); ++i){
      ma::EntityArray& a = r->newEntities[d][i];
      for (size_t i=0; i < a.getSize(); ++i)
        if (m->getType(a[i]) != apf::Mesh::VERTEX &&
            m->getModelType(m->toModel(a[i])) < m->getDimension()){
          int n = fs->getEntityShape(apf::Mesh::simplexTypes[d])->countNodes();
          int ni = fs->countNodesOn(d);
          apf::NewArray<double> c;
          crv::getBezierTransformationCoefficients(P,d,c);
          convertInterpolationPoints(m,a[i],n,ni,c);
        }
    }
  }
}

void splitEdges(ma::Adapt* a)
{
  assert(ma::checkFlagConsistency(a,1,ma::SPLIT));
  ma::Refine* r = a->refine;
  ma::resetCollection(r);
  ma::collectForTransfer(r);
  ma::addAllMarkedEdges(r);
  ma::splitElements(r);
  crv::snapRefineToBoundary(a);
  ma::processNewElements(r);
  ma::destroySplitElements(r);
  crv::repositionInterior(r);
  ma::forgetNewEntities(r);
}

static void refine(ma::Adapt* a)
{
  double t0 = PCU_Time();
  --(a->refinesLeft);
  long count = ma::markEdgesToSplit(a);
  if ( ! count) {
    return;
  }
  splitEdges(a);
  double t1 = PCU_Time();
  ma::print("split %li edges in %f seconds",count,t1-t0);
}

int getQualityTag(ma::Mesh* m, ma::Entity* e,
    ma::Entity* bdry)
{
  if (bdry == e) return 18;
  m->getType(bdry);
  int dim = apf::getDimension(m,bdry);
  apf::Downward down;
  int n = m->getDownward(e,dim,down);
  int index = apf::findIn(down,n,bdry);
  // set up the tag here;
  switch (dim) {
    case 0:
      return index+2;
    case 1:
      return index+8;
    case 2:
      return index+14;
    default:
      fail("invalid lower entity in quality check\n");
      break;
  }
  return -1;
}

long markBadQuality(Adapt* a)
{
  ma::Entity* e;
  long count = 0;
  ma::Mesh* m = a->mesh;
  int dimension = m->getDimension();
  ma::Iterator* it = m->begin(dimension);
  while ((e = m->iterate(it)))
  {
    /* this skip conditional is powerful: it affords us a
       3X speedup of the entire adaptation in some cases */
    if (crv::getFlag(a,e) > 0)
      continue;
    // going to use check validity here instead
    int qualityFlag = 1; // okay!
    // change this later
    int numInvalid = 0;
    ma::Entity* invalidEntity = 0;
    if(dimension == 2){
      ma::Entity* entities[6];
      numInvalid = checkTriValidity(m,e,entities,4);
      if(numInvalid) invalidEntity = entities[0];
    } else {
      ma::Entity* entities[14];
      numInvalid = checkTetValidity(m,e,entities,4);
      if(numInvalid) invalidEntity = entities[0];
    }
    if(numInvalid){
      qualityFlag = getQualityTag(m,e,invalidEntity);
    }
    if (qualityFlag > 1)
    {
      crv::setFlag(a,e,qualityFlag);
      if (m->isOwned(e))
        ++count;
    }
    else
      crv::setFlag(a,e,0);
  }
  m->end(it);
  return PCU_Add_Long(count);
}

int getFlag(Adapt* a, ma::Entity* e)
{
  return getFlags(a,e);
}

void setFlag(Adapt* a, ma::Entity* e, int flag)
{
  setFlags(a,e,flag);
}

void clearFlag(Adapt* a, ma::Entity* e)
{
  setFlags(a,e,0);
}

void adapt(ma::Input* in)
{
  in->shouldFixShape = true;
  in->shapeHandler = crv::getShapeHandler;
  ma::print("version 2.0 !");
  double t0 = PCU_Time();
  ma::validateInput(in);
  Adapt* a = new Adapt(in);
  ma::preBalance(a);
  crv::fixLargeBoundaryAngles(a);
  crv::fixInvalidEdges(a);
  for (int i=0; i < in->maximumIterations; ++i)
  {
    ma::print("iteration %d",i);
    ma::coarsen(a);
    ma::midBalance(a);
    crv::refine(a);
  }
  ma::postBalance(a);
  double t1 = PCU_Time();
  ma::print("mesh adapted in %f seconds",t1-t0);
  apf::printStats(a->mesh);
  crv::clearFlags(a);
  delete a;
  delete in;

}

}
