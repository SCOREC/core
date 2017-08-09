/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "crvAdapt.h"
#include "crvShape.h"
#include <apf.h>
#include <apfMesh.h>
#include <maBalance.h>
#include <maCoarsen.h>
#include <maShape.h>
#include <maSnap.h>
#include <maLayer.h>
#include <PCU.h>
#include <pcu_util.h>

namespace crv {

Adapt::Adapt(ma::Input* in)
: ma::Adapt(in)
{
  validityTag = mesh->createIntTag("crv_tags",1);
}

// rather than use the destructor to delete validityTag,
// this function takes care of it (since ~ma::Adapt() isn't virtual)
static void clearTags(Adapt* a)
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

static int getTags(Adapt* a, ma::Entity* e)
{
  ma::Mesh* m = a->mesh;
  if ( ! m->hasTag(e,a->validityTag))
    return 0; //we assume 0 is the default (unset) value for all tags
  int tags;
  m->getIntTag(e,a->validityTag,&tags);
  return tags;
}

static void setTags(Adapt* a, ma::Entity* e, int tags)
{
  a->mesh->setIntTag(e,a->validityTag,&tags);
}

void splitEdges(ma::Adapt* a)
{
  PCU_ALWAYS_ASSERT(ma::checkFlagConsistency(a,1,ma::SPLIT));
  ma::Refine* r = a->refine;
  ma::resetCollection(r);
  ma::collectForTransfer(r);
  ma::addAllMarkedEdges(r);
  ma::splitElements(r);
  ma::processNewElements(r);
  ma::destroySplitElements(r);
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

int getValidityTag(ma::Mesh* m, ma::Entity* e,
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

int markInvalidEntities(Adapt* a)
{
  ma::Entity* e;
  int count = 0;
  ma::Mesh* m = a->mesh;
  int dimension = m->getDimension();
  ma::Iterator* it = m->begin(dimension);
  Quality* qual = makeQuality(m,2);
  while ((e = m->iterate(it)))
  {
    /* this skip conditional is powerful: it affords us a
       3X speedup of the entire adaptation in some cases */
    int qualityTag = crv::getTag(a,e);
    if (qualityTag) continue;
    qualityTag = qual->checkValidity(e);
    if (qualityTag >= 2)
    {
      crv::setTag(a,e,qualityTag);
      if (m->isOwned(e))
        ++count;
    }
  }
  m->end(it);
  delete qual;
  return PCU_Add_Int(count);
}

int getTag(Adapt* a, ma::Entity* e)
{
  return getTags(a,e);
}

void setTag(Adapt* a, ma::Entity* e, int tag)
{
  setTags(a,e,tag);
}

void clearTag(Adapt* a, ma::Entity* e)
{
  setTags(a,e,0);
}
// use an identity configuration but with default fixing values
ma::Input* configureShapeCorrection(
    ma::Mesh* m, ma::SizeField* f,
    ma::SolutionTransfer* s)
{
  ma::Input* in = ma::configureIdentity(m,f,s);
  in->shouldFixShape = true;
  in->shouldSnap = in->mesh->canSnap();
  in->shouldTransferParametric = in->mesh->canSnap();
  return in;
}

static int fixInvalidElements(crv::Adapt* a)
{
  a->input->shouldForceAdaptation = true;
  int count = crv::fixLargeBoundaryAngles(a)
            + crv::fixInvalidEdges(a);
  int originalCount = count;
  int prev_count;
  int i = 0;
  do {
    if ( ! count)
      break;
    prev_count = count;
    count = crv::fixLargeBoundaryAngles(a)
          + crv::fixInvalidEdges(a);
    ++i;
  } while(count < prev_count);

  crv::fixLargeBoundaryAngles(a);
  ma::clearFlagFromDimension(a,ma::COLLAPSE | ma::BAD_QUALITY,1);
  a->input->shouldForceAdaptation = false;
  return originalCount - count;
}

static void flagCleaner(crv::Adapt* a)
{
  int dim = a->mesh->getDimension();

  for (int d = 0; d <= dim; d++) {
    ma::clearFlagFromDimension(a, ma::BAD_QUALITY, d);
    ma::clearFlagFromDimension(a, ma::OK_QUALITY, d);
  }
}

static bool isSimplexMesh(apf::Mesh2* m)
{
  int count = 0;
  if (m->getDimension() == 2)
    count = apf::countEntitiesOfType(m, apf::Mesh2::QUAD);
  else
    count = apf::countEntitiesOfType(m, apf::Mesh2::HEX) +
            apf::countEntitiesOfType(m, apf::Mesh2::PRISM) +
            apf::countEntitiesOfType(m, apf::Mesh2::PYRAMID);
  return (count == 0) ? true : false;
}

void adapt(ma::Input* in)
{
  PCU_ALWAYS_ASSERT_VERBOSE(isSimplexMesh(in->mesh),
      "Curved Adaptation expects an all simplex mesh!");
  std::string name = in->mesh->getShape()->getName();
  if(name != std::string("Bezier"))
    fail("mesh must be bezier to adapt\n");

  in->shapeHandler = crv::getShapeHandler;
  ma::print("Curved Adaptation Version 2.0 !");
  double t0 = PCU_Time();
  ma::validateInput(in);
  Adapt* a = new Adapt(in);
  ma::preBalance(a);

  fixInvalidElements(a);

  for (int i=0; i < in->maximumIterations; ++i)
  {
    ma::print("iteration %d",i);
    ma::coarsen(a);
    ma::midBalance(a);
    crv::refine(a);
    allowSplitCollapseOutsideLayer(a);
    flagCleaner(a); // all true-flags must be false before using markEntities
    fixCrvElementShapes(a);
  }

  allowSplitCollapseOutsideLayer(a);

  if (in->maximumIterations > 0) {
    fixInvalidElements(a);
    flagCleaner(a); // all true-flags must be false before using markEntities
    fixCrvElementShapes(a);
  }
  cleanupLayer(a);
  ma::printQuality(a);
  ma::postBalance(a);
  double t1 = PCU_Time();
  ma::print("mesh adapted in %f seconds",t1-t0);
  apf::printStats(a->mesh);
  crv::clearTags(a);
  delete a;
  delete in;
}


static void getStatsInMetricSpace(ma::Mesh* m, ma::SizeField* sf,
    std::vector<double> &edgeLengths,
    std::vector<double> &linearQualities,
    std::vector<double> &curvedQualities)
{
  int order = m->getShape()->getOrder();

  ma::Entity* e;
  ma::Iterator* it;

  // linear qualities
  it = m->begin(m->getDimension());
  while( (e = m->iterate(it)) )
    linearQualities.push_back(ma::measureElementQuality(m, sf, e));
  m->end(it);

  // curved qualities
  if (order == 1)
    curvedQualities = std::vector<double>(linearQualities.size(), 0.0);
  else {
    crv::Quality* qual = makeQuality(m, 2);
    it = m->begin(m->getDimension());
    while( (e = m->iterate(it)) ) {
      curvedQualities.push_back(qual->getQuality(e));
    }
    m->end(it);
  }

  // edge lengths
  it = m->begin(1);
  while( (e = m->iterate(it)) ) {
    apf::MeshElement* me = createMeshElement(m, e);
    ma::Matrix Q;
    sf->getTransform(me, ma::Vector(0., 0., 0.), Q);
    apf::destroyMeshElement(me);
    edgeLengths.push_back(ma::qMeasure(m, e, Q));
  }
  m->end(it);
}

static void getStatsInPhysicalSpaces(ma::Mesh* m, ma::SizeField* sf,
    std::vector<double> &edgeLengths,
    std::vector<double> &linearQualities,
    std::vector<double> &curvedQualities)
{
  int order = m->getShape()->getOrder();

  ma::Entity* e;
  ma::Iterator* it;

  // linear qualities
  it = m->begin(m->getDimension());
  while( (e = m->iterate(it)) ) {
    if (m->getType(e) == apf::Mesh::TRIANGLE) {
      ma::Vector p[3];
      ma::getVertPoints(m, e, p);
      double l[3];
      for (int i = 0; i < 3; i++)
	l[i] = (p[(i+1)%3] - p[i]).getLength();
      double A = 0.5 * apf::cross(p[1] - p[0], p[2] - p[0])[2];
      double s = 0;
      for (int i = 0; i < 3; i++)
	s += l[i] * l[i];
      double lq;
      if (A < 0)
	lq = -48. * (A*A) / (s*s);
      else
	lq =  48. * (A*A) / (s*s);
      linearQualities.push_back(lq);
    }
    else if (m->getType(e) == apf::Mesh::TET){
      ma::Vector p[4];
      ma::getVertPoints(m, e, p);
      linearQualities.push_back(ma::measureLinearTetQuality(p));
    }
  }
  m->end(it);

  // curved qualities
  if (order == 1)
    curvedQualities = std::vector<double>(linearQualities.size(), 0.0);
  else {
    crv::Quality* qual = makeQuality(m, 2);
    it = m->begin(m->getDimension());
    while( (e = m->iterate(it)) ) {
      curvedQualities.push_back(qual->getQuality(e));
    }
    m->end(it);
  }

  // edge lengths
  it = m->begin(1);
  while( (e = m->iterate(it)) ) {
    apf::MeshElement* me = createMeshElement(m, e);
    ma::Matrix Q = ma::Matrix(1.0, 0.0, 0.0,
			      0.0, 1.0, 0.0,
			      0.0, 0.0, 1.0);;
    apf::destroyMeshElement(me);
    edgeLengths.push_back(ma::qMeasure(m, e, Q));
  }
  m->end(it);
}

/** \brief Measures entity related quantities for a given mesh
  \details  quantities include normalized edge length, linear quality
  and curved quality. The values can be computed in both metric (if
  inMetric = true) and physical (if inMetric = false) spaces.*/
void stats(ma::Input* in,
    std::vector<double> &edgeLengths,
    std::vector<double> &linearQualities,
    std::vector<double> &curvedQualities,
    bool inMetric)
{
  PCU_ALWAYS_ASSERT_VERBOSE(isSimplexMesh(in->mesh),
      "crv::stats expects an all simplex mesh!");
  std::string name = in->mesh->getShape()->getName();
  int order = in->mesh->getShape()->getOrder();
  /* if(name != std::string("Bezier")) */
  /*   fail("mesh must be bezier to adapt\n"); */

  edgeLengths.clear();
  linearQualities.clear();
  curvedQualities.clear();

  /* if (order > 1) */
  /*   in->shapeHandler = crv::getShapeHandler(a); */
  /* else */

  ma::validateInput(in);
  Adapt* a = new Adapt(in);

  ma::Mesh* m = a->mesh;
  ma::SizeField* sf = a->sizeField;

  if (inMetric)
    getStatsInMetricSpace(m, sf, edgeLengths, linearQualities, curvedQualities);
  else
    getStatsInPhysicalSpaces(m, sf, edgeLengths, linearQualities, curvedQualities);

  delete a;
  delete in;
}

}
