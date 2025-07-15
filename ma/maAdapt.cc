/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include <lionPrint.h>
#include "maAdapt.h"
#include "maInput.h"
#include "maTables.h"
#include "maCoarsen.h"
#include "maRefine.h"
#include "maSolutionTransfer.h"
#include "maShape.h"
#include "maShapeHandler.h"
#include "maLayer.h"
#include <apf.h>
#include <cfloat>
#include <pcu_util.h>
#include <stdarg.h>

namespace ma {

Adapt::Adapt(Input* in)
{
  input = in;
  mesh = in->mesh;
  setupFlags(this);
  setupQualityCache(this);
  deleteCallback = 0;
  buildCallback = 0;
  sizeField = in->sizeField;
  solutionTransfer = in->solutionTransfer;
  refine = new Refine(this);
  if (in->shapeHandler){
    shape = in->shapeHandler(this);
  } else
    shape = getShapeHandler(this);
  if (in->shouldCoarsen)
    coarsensLeft = in->maximumIterations;
  else
    coarsensLeft = 0;
  refinesLeft = in->maximumIterations;
  resetLayer(this);
  if (hasLayer)
    checkLayerShape(mesh, "input mesh");
}

Adapt::~Adapt()
{
  clearFlags(this);
  clearQualityCache(this);
  delete refine;
  delete shape;
}

void setupFlags(Adapt* a)
{
  a->flagsTag = a->mesh->createIntTag("ma_flags",1);
}

void clearFlags(Adapt* a)
{
  Mesh* m = a->mesh;
  Entity* e;
  for (int d=0; d <= 3; ++d)
  {
    Iterator* it = m->begin(d);
    while ((e = m->iterate(it)))
      if (m->hasTag(e,a->flagsTag))
        m->removeTag(e,a->flagsTag);
    m->end(it);
  }
  m->destroyTag(a->flagsTag);
}

int getFlags(Adapt* a, Entity* e)
{
  Mesh* m = a->mesh;
  if ( ! m->hasTag(e,a->flagsTag))
    return 0; //we assume 0 is the default value for all flags
  int flags;
  m->getIntTag(e,a->flagsTag,&flags);
  return flags;
}

void setFlags(Adapt* a, Entity* e, int flags)
{
  a->mesh->setIntTag(e,a->flagsTag,&flags);
}

bool getFlag(Adapt* a, Entity* e, int flag)
{
  return flag & getFlags(a,e);
}

void setFlag(Adapt* a, Entity* e, int flag)
{
  int flags = getFlags(a,e);
  flags |= flag;
  setFlags(a,e,flags);
}

void setFlagMatched(Adapt* a, Entity* e, int flag)
{
  if (a->mesh->hasMatching()) {
    apf::Matches matches;
    a->mesh->getMatches(e, matches);
    APF_ITERATE(apf::Matches, matches, it) {
      PCU_ALWAYS_ASSERT(it->peer == a->mesh->getPCU()->Self());
      setFlag(a, it->entity, flag);
    }
  }
  setFlag(a, e, flag);
}

void clearFlag(Adapt* a, Entity* e, int flag)
{
  int flags = getFlags(a,e);
  flags &= ~flag;
  setFlags(a,e,flags);
}

void clearFlagMatched(Adapt* a, Entity* e, int flag)
{
  if (a->mesh->hasMatching()) {
    apf::Matches matches;
    a->mesh->getMatches(e, matches);
    APF_ITERATE(apf::Matches, matches, it) {
      PCU_ALWAYS_ASSERT(it->peer == a->mesh->getPCU()->Self());
      clearFlag(a, it->entity, flag);
    }
  }
  clearFlag(a, e, flag);
}

void clearFlagFromDimension(Adapt* a, int flag, int dimension)
{
  Mesh* m = a->mesh;
  Iterator* it = m->begin(dimension);
  Entity* e;
  while ((e = m->iterate(it)))
    clearFlag(a,e,flag);
  m->end(it);
}

void setupQualityCache(Adapt* a)
{
  a->qualityCache = a->mesh->createDoubleTag("ma_qual_cache",1);
}

void clearQualityCache(Adapt* a)
{
  Mesh* m = a->mesh;
  Entity* e;
  // only faces and regions can have the quality tag
  for (int d=2; d <= 3; ++d)
  {
    Iterator* it = m->begin(d);
    while ((e = m->iterate(it)))
      if (m->hasTag(e,a->qualityCache))
        m->removeTag(e,a->qualityCache);
    m->end(it);
  }
  m->destroyTag(a->qualityCache);
}

double getCachedQuality(Adapt* a, Entity* e)
{
  Mesh* m = a->mesh;
  int type = m->getType(e);
  int ed = apf::Mesh::typeDimension[type];
  PCU_ALWAYS_ASSERT(ed == 2 || ed == 3);
  if ( ! m->hasTag(e,a->qualityCache))
    return 0.0; //we assume 0.0 is the default value for all qualities
  double qual;
  m->getDoubleTag(e,a->qualityCache,&qual);
  return qual;
}

void setCachedQuality(Adapt* a, Entity* e, double q)
{
  Mesh* m = a->mesh;
  int type = m->getType(e);
  int ed = apf::Mesh::typeDimension[type];
  PCU_ALWAYS_ASSERT(ed == 2 || ed == 3);
  m->setDoubleTag(e,a->qualityCache,&q);
}

void destroyElement(Adapt* a, Entity* e)
{
  Mesh* m = a->mesh;
  int dim = getDimension(m,e);
  if (dim < m->getDimension())
  { //destruction is a no-op if this entity still supports
    //higher-order ones
    if (m->hasUp(e))
      return;
  }
  Downward down;
  int nd = 0;
  if (dim > 0)
    nd = m->getDownward(e,dim-1,down);
  if (a->deleteCallback) a->deleteCallback->call(e);
  m->destroy(e);
  /* destruction applies recursively to the closure of the entity */
  if (dim > 0)
    for (int i=0; i < nd; ++i)
      destroyElement(a,down[i]);
}

DeleteCallback::DeleteCallback(Adapt* a)
{
  this->adapt = a;
  a->deleteCallback = this;
}

DeleteCallback::~DeleteCallback()
{
  this->adapt->deleteCallback = 0;
}

bool checkFlagConsistency(Adapt* a, int dimension, int flag)
{
  Mesh* m = a->mesh;
  apf::Sharing* sh = apf::getSharing(m);
  m->getPCU()->Begin();
  Entity* e;
  Iterator* it = m->begin(dimension);
  while ((e = m->iterate(it))) {
    apf::CopyArray others;
    sh->getCopies(e, others);
    if (!others.getSize())
      continue;
    bool value = getFlag(a, e, flag);
    APF_ITERATE(apf::CopyArray, others, rit) {
      m->getPCU()->Pack(rit->peer, rit->entity);
      m->getPCU()->Pack(rit->peer, value);
    }
  }
  m->end(it);
  m->getPCU()->Send();
  bool ok = true;
  while (m->getPCU()->Receive()) {
    m->getPCU()->Unpack(e);
    bool value;
    m->getPCU()->Unpack(value);
    if(value != getFlag(a,e,flag))
      ok = false;
  }
  delete sh;
  return ok;
}

double getDistance(Adapt* a, Entity* v[2])
{
  return (getPosition(a->mesh,v[1]) - getPosition(a->mesh,v[0])).getLength();
}

double getDistance(Adapt* a, Entity* v0, Entity* v1)
{
  Entity* v[2] = {v0,v1};
  return getDistance(a,v);
}

int getClosestPair(Adapt* a, Entity* (*pairs)[2], int n)
{
  double min = DBL_MAX;
  double distance;
  int min_i = -1;
  for (int i=0; i < n; ++i)
    if ((distance = getDistance(a,pairs[i])) < min)
    {
      min = distance;
      min_i = i;
    }
  return min_i;
}

/* marks entities of a dimension for which the predicate
   returns true with the true flag, and uses the false
   flag to prevent duplicate checks of the same entity.
   Per the workings of an ma::Operator, it expects the
   true flag to be cleared from all entities, so it
   always re-evaluates those entities.

   returns the total global number of marked entities,
   counting shared entities once.
*/
long markEntities(
    Adapt* a,
    int dimension,
    Predicate& predicate,
    int trueFlag,
    int setFalseFlag,
    int allFalseFlags)
{
  if (!allFalseFlags) allFalseFlags = setFalseFlag;
  Entity* e;
  long count = 0;
  Mesh* m = a->mesh;
  Iterator* it = m->begin(dimension);
  while ((e = m->iterate(it)))
  {
    PCU_ALWAYS_ASSERT( ! getFlag(a,e,trueFlag));
    /* this skip conditional is powerful: it affords us a
       3X speedup of the entire adaptation in some cases */
    if (allFalseFlags & getFlags(a,e))
      continue;
    if (predicate(e))
    {
      setFlag(a,e,trueFlag);
      if (a->mesh->isOwned(e))
        ++count;
    }
    else
      setFlag(a,e,setFalseFlag);
  }
  m->end(it);
  return m->getPCU()->Add<long>(count);
}

void NewEntities::reset()
{
  entities.clear();
}

void NewEntities::addEntity(Entity* e)
{
  entities.push_back(e);
}

void NewEntities::call(Entity* e)
{
  addEntity(e);
}

void NewEntities::retrieve(EntityArray& a)
{
  a.setSize(entities.size());
  for (size_t i=0; i < entities.size(); ++i)
    a[i] = entities[i];
}

Cavity::Cavity()
{
  shouldTransfer = false;
  shouldFit = false;
  adapter = 0;
  solutionTransfer = 0;
  shape = 0;
  shouldTransferSizeField = false;
}

void Cavity::init(Adapt* a)
{
  adapter = a;
  solutionTransfer = a->solutionTransfer;
  shape = a->shape;
  Mesh* m = a->mesh;
  shouldTransfer = false;
  shouldFit = false;
  for (int d=1; d <= m->getDimension(); ++d) {
    if (solutionTransfer->hasNodesOn(d))
      shouldTransfer = true;
    if (shape->hasNodesOn(d))
      shouldFit = true;
    if (adapter->sizeField->hasNodesOn(d))
      shouldTransferSizeField = true;
  }
}

void Cavity::beforeBuilding()
{
  if (shouldTransfer || shouldFit)
  {
    newEntities.reset();
    setBuildCallback(adapter,&newEntities);
  }
}

void Cavity::afterBuilding()
{
  if (shouldTransfer || shouldFit)
    clearBuildCallback(adapter);
}

void Cavity::beforeTrying()
{
  if (shouldFit)
  {
    newEntities.reset();
    setBuildCallback(adapter,&newEntities);
  }
}

void Cavity::afterTrying()
{
  if (shouldFit)
    clearBuildCallback(adapter);
}

void Cavity::transfer(EntityArray& oldElements)
{
  if (shouldTransfer)
  {
    EntityArray a;
    newEntities.retrieve(a);
    solutionTransfer->onCavity(oldElements,a);
  }
  if (shouldTransferSizeField)
  {
    EntityArray a;
    newEntities.retrieve(a);
    adapter->sizeField->onCavity(oldElements,a);
  }
}

void Cavity::fit(EntityArray& oldElements)
{
  if (shouldFit)
  {
    EntityArray a;
    newEntities.retrieve(a);
    shape->onCavity(oldElements,a);
  }
}

Entity* buildVertex(
    Adapt* a,
    Model* c,
    Vector const& point,
    Vector const& param)
{
  Entity* v = a->mesh->createVertex(c,point,param);
  if (a->buildCallback)
    a->buildCallback->call(v);
  return v;
}

Entity* buildElement(
    Adapt* a,
    Model* c,
    int type,
    Entity** verts)
{
  return apf::buildElement(a->mesh,c,type,verts,a->buildCallback);
}

Entity* rebuildElement(
    Adapt* a,
    Entity* original,
    Entity* oldVert,
    Entity* newVert)
{
  return rebuildElement(a->mesh,original,oldVert,newVert,a->buildCallback);
}

void setBuildCallback(Adapt* a, apf::BuildCallback* cb)
{
  PCU_ALWAYS_ASSERT(a->buildCallback==0);
  a->buildCallback = cb;
}

void clearBuildCallback(Adapt* a)
{
  a->buildCallback = 0;
}

void print(pcu::PCU *PCUObj, const char* format, ...)
{
  if (PCUObj->Self())
    return;
  lion_oprint(1,"\nMeshAdapt: ");
  va_list ap;
  va_start(ap,format);
  lion_voprint(1,format,ap);
  va_end(ap);
  lion_oprint(1,"\n");
}

void setFlagOnClosure(Adapt* a, Entity* element, int flag)
{
  Mesh* m = a->mesh;
  int D = getDimension(m,element);
  for (int d=0; d <= D; ++d)
  {
    Downward down;
    int nd = m->getDownward(element,d,down);
    for (int i=0; i < nd; ++i)
      setFlag(a,down[i],flag);
  }
}

void syncFlag(Adapt* a, int dimension, int flag)
{
  Mesh* m = a->mesh;
  apf::Sharing* sh = apf::getSharing(m);
  m->getPCU()->Begin();
  Entity* e;
  Iterator* it = m->begin(dimension);
  while ((e = m->iterate(it))) {
    if (!getFlag(a, e, flag))
      continue;
    apf::CopyArray others;
    sh->getCopies(e, others);
    APF_ITERATE(apf::CopyArray, others, rit)
      m->getPCU()->Pack(rit->peer, rit->entity);
  }
  m->end(it);
  m->getPCU()->Send();
  while (m->getPCU()->Receive()) {
    m->getPCU()->Unpack(e);
    setFlag(a,e,flag);
  }
  delete sh;
}

HasTag::HasTag(Mesh* m, Tag* t)
{
  mesh = m;
  tag = t;
}

bool HasTag::operator()(Entity* e)
{
  return mesh->hasTag(e, tag);
}

HasFlag::HasFlag(Adapt* a, int f)
{
  adapter = a;
  flag = f;
}

bool HasFlag::operator()(Entity* e)
{
  return getFlag(adapter, e, flag);
}

}
