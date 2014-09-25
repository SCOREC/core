/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <PCU.h>
#include "apfCavityOp.h"
#include "apf.h"
#include "apfMesh2.h"

namespace apf {

CavityOp::CavityOp(Mesh* m, bool cm):
  mesh(m),
  isRequesting(false),
  canModify(cm)
{
}

/* these functions are over in a corner because they
   require Mesh2 functionality. This is more or
   less ok because deletion is also a Mesh2-only
   feature. */
void CavityOp::applyLocallyWithModification(int d)
{
  Mesh2* mesh2 = static_cast<Mesh2*>(mesh);
  MeshEntity* e;
  isRequesting = false;
  for (this->iterator = mesh2->begin(d);
       ! mesh2->isDone(this->iterator);)
  {
    e = mesh2->deref(this->iterator);
    if (mesh2->isOwned(e))
    {
      Outcome o = setEntity(e);
      if (o == OK)
      { 
        movedByDeletion = false;
        apply();
        if (movedByDeletion)
          continue;
      }
    }
    mesh2->increment(this->iterator);
  }
  mesh2->end(this->iterator);
  /* request any non-local cavities.
     note: it is critical that this loop
     be separated from the one above for
     mesh-modifying operators, since an apply()
     call could change other overlapping cavities */
  isRequesting = true;
  this->iterator = mesh2->begin(d);
  while ((e = mesh2->iterate(this->iterator)))
  {
    if ( ! mesh2->isOwned(e)) continue;
    setEntity(e);
  }
  mesh2->end(this->iterator);
}

void CavityOp::preDeletion(MeshEntity* e)
{
  Mesh2* mesh2 = static_cast<Mesh2*>(mesh);
  if (( ! mesh2->isDone(this->iterator))&&
      (e == mesh2->deref(this->iterator)))
  {
    /* save the iterator from invalidation */
    mesh2->increment(this->iterator);
    movedByDeletion = true;
  }
}

void CavityOp::applyLocallyWithoutModification(int d)
{
  /* if the mesh is static then computation and
     cavity requests can happen in the same loop */
  MeshIterator* entities = mesh->begin(d);
  MeshEntity* e;
  isRequesting = true;
  while ((e = mesh->iterate(entities)))
  {
    if ( ! mesh->isOwned(e))
      continue;
    Outcome o = setEntity(e);
    if (o == OK)
      apply();
  }
  mesh->end(entities);
}

void CavityOp::applyToDimension(int d)
{
  /* the iteration count of this loop is hard to predict,
   * but typical cavity definitions should cause a small
   * constant number of iterations that does not grow
   * with parallelism 
   */
  do {
    /* apply the operator to all local cavities
       and request missing cavity elements */
    if (this->canModify)
      this->applyLocallyWithModification(d);
    else
      this->applyLocallyWithoutModification(d);
    /* this is the exit of the loop:
       tryToPull will return false when no requests
       were made by any process, which should imply
       that all mesh entities that needed to be operated
       on have been. */
  } while (tryToPull());
}

bool CavityOp::requestLocality(MeshEntity** entities, int count)
{
  bool areLocal = true;
  for (int i=0; i < count; ++i)
    if (mesh->isShared(entities[i]))
      areLocal = false;
  if (isRequesting && ( ! areLocal))
    for (int i=0; i < count; ++i)
      requests.insert(requests.end(),entities,entities+count);
  return areLocal;
}

bool CavityOp::sendPullRequests(std::vector<PullRequest>& received)
{
  int done = requests.empty();
  PCU_Min_Ints(&done,1);
  if (done) return false;
  /* throw in the local pull requests */
  int self = PCU_Comm_Self();
  received.reserve(requests.size());
  APF_ITERATE(Requests,requests,it)
  {
    PullRequest request;
    request.to = self;
    request.e = *it;
    received.push_back(request);
  }
  /* now communicate the rest */
  PCU_Comm_Begin();
  APF_ITERATE(Requests,requests,it)
  {
    Copies remotes;
    mesh->getRemotes(*it,remotes);
    APF_ITERATE(Copies,remotes,rit)
    {
      int remotePart = rit->first;
      MeshEntity* remoteEntity = rit->second;
      PCU_COMM_PACK(remotePart,remoteEntity);
    }
  }
  requests.clear();
  PCU_Comm_Send();
  while (PCU_Comm_Listen())
  {
    PullRequest request;
    request.to = PCU_Comm_Sender();
    while ( ! PCU_Comm_Unpacked())
    {
      PCU_COMM_UNPACK(request.e);
      received.push_back(request);
    }
  }
  return true;
}

/* in the case of competing requests for an
   element, the winner is the maximum part ID.
   This has subtle consequences.
   The whole pulling algorithm favors the highest
   part ID, which could in theory pull the whole
   mesh to it if everything it gets is designated
   as "owned".
   FMDB's ownership rule is poor-to-rich, which will
   naturally counteract this effect and thus is a good rule.
   Other databases may use MIN part ID as their ownership rule,
   which is why this one is MAX (just in case) */

static int ownership(int a, int b)
{
  return std::max(a,b);
}

static void markElement(
    Migration* plan,
    MeshEntity* element,
    int requester)
{
  if (plan->has(element))
  {
    int current = plan->sending(element);
    plan->send(element,ownership(current,requester));
  }
  else
    plan->send(element,requester);
}

static void markElements(
    Migration* plan,
    MeshEntity* e,
    int requester)
{
  Mesh* m = plan->getMesh();
  Adjacent a;
  m->getAdjacent(e,m->getDimension(),a);
  for (size_t i=0; i < a.getSize(); ++i)
    markElement(plan,a[i],requester);
}

bool CavityOp::tryToPull()
{
  std::vector<PullRequest> pulls;
  if ( ! sendPullRequests(pulls))
    return false;
  Migration* plan = new Migration(mesh);
  for (std::size_t i=0; i < pulls.size(); ++i)
    markElements(plan,pulls[i].e,pulls[i].to);
  mesh->migrate(plan); //plan deleted here
  return true;
}

} //namespace apf
