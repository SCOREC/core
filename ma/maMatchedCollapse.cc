/****************************************************************************** 

  Copyright 2015 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "maMatchedCollapse.h"
#include "maAdapt.h"
#include <pcu_util.h>

#include <algorithm>

namespace ma {

Rebuild::Rebuild(Entity* a, Entity* b):
  e(a),original(b)
{
}

bool Rebuild::operator<(Rebuild const& other) const
{
  if (original != other.original)
    return original < other.original;
  if (e != other.e)
    return e < other.e;
  return false;
}

bool Rebuild::operator==(Rebuild const& other) const
{
  return original == other.original &&
         e == other.e;
}

Rebuilds::Rebuilds(Mesh* m):
  mesh(m)
{
}

void Rebuilds::rebuilt(Entity* e, Entity* original)
{
  int dim = apf::getDimension(mesh, e);
  if (0 < dim && dim < mesh->getDimension())
    v.push_back(Rebuild(e, original));
}

void Rebuilds::reset()
{
  v.clear();
}

struct IsFalseRebuild {
  bool operator()(Rebuild const& r) const
  {
    return r.original == r.e;
  }
};

void Rebuilds::match(apf::Sharing* sh)
{
  /* the ma::rebuildElement call will produce more logs than we want:
     1) entities which don't really change show up as rebuilds of themselves,
        i.e. edges along the boundary of the cavity. TODO: This seems
        unnecessary. Double check before removing it.
     2) rebuilds required my multiple elements show up multiple times
     below we filter out each of these cases to obtain the real rebuilds
     we're interested in */
  v.erase(std::remove_if(v.begin(), v.end(), IsFalseRebuild()), v.end());
  std::sort(v.begin(), v.end());
  v.erase(std::unique(v.begin(), v.end()), v.end());
  for (unsigned i = 0; i < v.size(); ++i) {
    Entity* orig = v[i].original;
    Entity* gen = v[i].e;
    PCU_ALWAYS_ASSERT(mesh->getType(orig) == mesh->getType(gen));
    apf::CopyArray orig_matches;
    sh->getCopies(orig, orig_matches);
    for (unsigned j = 0; j < orig_matches.getSize(); ++j) {
      PCU_ALWAYS_ASSERT(orig_matches[j].peer == mesh->getPCU()->Self());
      Entity* gen_match_j = 0;
      for (unsigned k = 0; k < v.size(); ++k)
        if (v[k].original == orig_matches[j].entity)
          gen_match_j = v[k].e;
      PCU_ALWAYS_ASSERT(gen_match_j);
      mesh->addMatch(gen, mesh->getPCU()->Self(), gen_match_j);
    }
  }
}

MatchedCollapse::MatchedCollapse(Adapt* a):
  adapt(a),
  rebuilds(a->mesh)
{
  mesh = a->mesh;
  sharing = apf::getSharing(mesh);
}

MatchedCollapse::~MatchedCollapse()
{
  delete sharing;
}

void MatchedCollapse::setEdge(Entity* e)
{
  collapses.setSize(0);
  collapses.setSize(1);
  collapses[0].Init(adapt);
  bool ok = collapses[0].setEdge(e);
  PCU_ALWAYS_ASSERT(ok);
}

bool MatchedCollapse::requestLocality(apf::CavityOp* o)
{
  return collapses[0].requestLocality(o);
}

void MatchedCollapse::setEdges()
{
  Entity* e = collapses[0].edge;
  collapses.setSize(0);
  apf::CopyArray copies;
  sharing->getCopies(e, copies);
  APF_ITERATE(apf::CopyArray, copies, it) {
    PCU_ALWAYS_ASSERT(it->peer == mesh->getPCU()->Self());
    PCU_ALWAYS_ASSERT(it->entity != e);
  }
  collapses.setSize(copies.getSize() + 1);
  for (unsigned i = 0; i < collapses.getSize(); ++i)
    collapses[i].Init(adapt);
  bool ok = collapses[0].setEdge(e);
  PCU_ALWAYS_ASSERT(ok);
  for (unsigned i = 0; i < copies.getSize(); ++i) {
    ok = collapses[i + 1].setEdge(copies[i].entity);
    PCU_ALWAYS_ASSERT(ok);
  }
}

bool MatchedCollapse::checkTopo2()
{
  for (unsigned i = 0; i < collapses.getSize(); ++i)
    if (!checkEdgeCollapseTopology(adapt, collapses[i].edge))
      return false;
  for (unsigned i = 0; i < collapses.getSize(); ++i)
    collapses[i].setVerts();
  /* assert that the vertices being collapsed are matched copies
     of one another, and that they have no off-processor
     copies of any kind */
  apf::CopyArray copies;
  sharing->getCopies(collapses[0].vertToCollapse, copies);
  for (unsigned i = 1; i < collapses.getSize(); ++i) {
    unsigned j;
    bool ok = false;
    for (j = 0; j < copies.getSize(); ++j) {
      PCU_ALWAYS_ASSERT(copies[j].peer == mesh->getPCU()->Self());
      if (copies[j].entity == collapses[i].vertToCollapse) {
        ok = true;
        break;
      }
    }
    PCU_ALWAYS_ASSERT(ok);
  }
  /* things go quite badly if the sub-collapse cavities overlap.
     one cheap way to check is to see if there exist mesh
     edges between the vertices to collapse */
  for (unsigned i = 1; i < collapses.getSize(); ++i)
    if (findEdge(mesh,
          collapses[0].vertToCollapse,
          collapses[i].vertToCollapse))
      return false;
  return true;
}

bool MatchedCollapse::checkTopo()
{
  if (!checkTopo2()) {
    unmark();
    return false;
  }
  return true;
}

void MatchedCollapse::unmark()
{
  // TODO: Something seems to be wrong with this.
  // -- in the second for loop the condition doesn't
  // seem to be the write one [look, for example, to
  // the corresponding function in maCollapse.cc]
  for (unsigned i = 0; i < collapses.getSize(); ++i)
    clearFlag(adapt, collapses[i].edge, COLLAPSE);
  bool required[2] = {false,false};
  for (unsigned i = 0; i < collapses.getSize(); ++i) {
    Entity* v[2];
    mesh->getDownward(collapses[i].edge, 0, v);
    for (unsigned j = 0; j < 2; ++j)
      if (isRequiredForAnEdgeCollapse(adapt, v[j]))
        required[j] = true;
  }
  for (unsigned i = 0; i < collapses.getSize(); ++i) {
    Entity* v[2];
    mesh->getDownward(collapses[i].edge, 0, v);
    for (unsigned j = 0; j < 2; ++j)
      if (required[j])
        clearFlag(adapt, v[j], COLLAPSE);
  }
}

void MatchedCollapse::cancel()
{
  for (unsigned i = 0; i < collapses.getSize(); ++i)
    collapses[i].destroyNewElements();
  unmark();
}

bool MatchedCollapse::overlapsSelf()
{
  for (unsigned i = 0; i < collapses.getSize(); ++i)
  for (unsigned j = 0; j < collapses.getSize(); ++j) {
    if (i == j)
      continue;
    EntityArray old[2];
    collapses[i].getOldElements(old[0]);
    collapses[j].getOldElements(old[1]);
    for (size_t k = 0; k < old[0].getSize(); ++k)
    for (size_t l = 0; l < old[1].getSize(); ++l)
      if (old[0][k] == old[1][l])
        return true;
  }
  return false;
}

bool MatchedCollapse::tryThisDirection(double qualityToBeat)
{
  for (unsigned i = 0; i < collapses.getSize(); ++i)
    collapses[i].computeElementSets();
  if (overlapsSelf()) {
    cancel();
    return false;
  }
  rebuilds.reset();
  for (unsigned i = 0; i < collapses.getSize(); ++i)
    collapses[i].rebuildCallback = &rebuilds;
  bool ok = true;
  for (unsigned i = 0; i < collapses.getSize(); ++i)
    if (!collapses[i].tryThisDirectionNoCancel(qualityToBeat))
      ok = false;
  if (ok)
    rebuilds.match(sharing);
  else
    cancel();
  return ok;
}

bool MatchedCollapse::tryBothDirections(double qualityToBeat)
{
  if (tryThisDirection(qualityToBeat))
    return true;
  for (unsigned i = 0; i < collapses.getSize(); ++i)
    if (!getFlag(adapt, collapses[i].vertToKeep, COLLAPSE))
      return false;
  for (unsigned i = 0; i < collapses.getSize(); ++i)
    std::swap(collapses[i].vertToKeep, collapses[i].vertToCollapse);
  return tryThisDirection(qualityToBeat);
}

void MatchedCollapse::destroyOldElements()
{
  for (unsigned i = 0; i < collapses.getSize(); ++i)
    collapses[i].destroyOldElements();
}

}
