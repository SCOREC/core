#include "phFilterMatching.h"
#include "phModelGeometry.h"
#include "phAxisymmetry.h"
#include <vector>
#include <map>
#include <set>
#include <apf.h>
#include <PCU.h>
#include <cassert>
#include <cstdlib>
#include <iostream>

namespace ph {

typedef std::vector<apf::Matches> SavedMatches;
typedef std::set<gmi_ent*> ModelSet;
typedef std::map<gmi_ent*, ModelSet> ModelMatching;

void saveMatches(apf::Mesh* m, int dim, SavedMatches& sm);
void restoreMatches(apf::Mesh2* m, int dim, SavedMatches& sm);
void getFullAttributeMatching(gmi_model* m, BCs& bcs, ModelMatching& mm);
void filterMatching(apf::Mesh2* m, ModelMatching& mm, int dim);

void saveMatches(apf::Mesh* m, int dim, SavedMatches& sm)
{
  sm.resize(m->count(dim));
  apf::MeshIterator* it = m->begin(dim);
  apf::MeshEntity* e;
  unsigned i = 0;
  while ((e = m->iterate(it))) {
    m->getMatches(e, sm[i]);
    ++i;
  }
  m->end(it);
}

void restoreMatches(apf::Mesh2* m, int dim, SavedMatches& sm)
{
  apf::MeshIterator* it = m->begin(dim);
  apf::MeshEntity* e;
  unsigned i = 0;
  while ((e = m->iterate(it))) {
    if (sm[i].getSize()) {
      m->clearMatches(e);
      APF_ITERATE(apf::Matches, sm[i], mit)
        m->addMatch(e, mit->peer, mit->entity);
    }
    ++i;
  }
  m->end(it);
}

static void completeEntMatching(gmi_ent* e, ModelMatching& mm, ModelSet& visited)
{
  if (visited.count(e))
    return;
  visited.insert(e);
  if (!mm.count(e))
    return;
  ModelSet& adj = mm[e];
  APF_ITERATE(ModelSet, adj, it)
    completeEntMatching(*it, mm, visited);
}

static void completeMatching(ModelMatching& mm)
{
  ModelMatching mm2;
  APF_ITERATE(ModelMatching, mm, it) {
    ModelSet reachable;
    completeEntMatching(it->first, mm, reachable);
    reachable.erase(it->first);
    mm2[it->first] = reachable;
  }
  mm = mm2;
}

static void addMatch(gmi_ent* a, gmi_ent* b, ModelMatching& mm)
{
  /* the hinge edge in axisymmetry can show up matching itself */
  if (a == b)
    return;
  mm[a].insert(b);
  mm[b].insert(a);
}

static void getAttributeMatching(gmi_model* gm, BCs& bcs, ModelMatching& mm)
{
  std::string name = "periodic slave";
  if (!haveBC(bcs, name))
    return;
  FieldBCs& fbcs = bcs.fields[name];
  APF_ITERATE(FieldBCs::Set, fbcs.bcs, it) {
    BC* bc = *it;
    gmi_ent* e = gmi_find(gm, bc->dim, bc->tag);
    double* val = bc->eval(apf::Vector3(0,0,0));
    int otherTag = *val;
    gmi_ent* oe = gmi_find(gm, bc->dim, otherTag);
    addMatch(e, oe, mm);
  }
}

static double getDistanceWithFrame(gmi_model* gm, gmi_ent* e, gmi_ent* oe,
    apf::Frame const& frame)
{
  apf::Vector3 pt = getCenter(gm, e);
  apf::Vector3 opt = getCenter(gm, oe);
  return (frame * pt - opt).getLength();
}

static void closeFaceMatchingWithFrame(gmi_model* gm, gmi_ent* f, gmi_ent* of,
    ModelMatching& mm, apf::Frame const& frame)
{
  for (int dim = 0; dim < 2; ++dim) {
    gmi_set* s = gmi_adjacent(gm, f, dim);
    gmi_set* os = gmi_adjacent(gm, of, dim);
    /* these faces are matched, they should have the
       same layout of bounding entities with only geometric differences */
    assert(s->n == os->n);
    /* warning! this is an O(N^2) comparison.
       if it takes up a large part of your runtime,
       you should rethink your life. */
    for (int i = 0; i < s->n; ++i) {
      double minDist = getDistanceWithFrame(
          gm, s->e[i], os->e[0], frame);
      gmi_ent* closest = os->e[0];
      for (int j = 1; j < os->n; ++j) {
        double dist = getDistanceWithFrame(
            gm, s->e[i], os->e[j], frame);
        if (dist < minDist) {
          minDist = dist;
          closest = os->e[j];
        }
      }
      addMatch(s->e[i], closest, mm);
    }
    gmi_free_set(s);
    gmi_free_set(os);
  }
}

static void closeFaceMatching(gmi_model* gm, gmi_ent* f, ModelMatching& mm)
{
  assert(mm[f].size() == 1);
  gmi_ent* of = *(mm[f].begin());
  if (f < of)
    return;
  assert(f != of);
  /* the key assumptions are these:
     1) both faces are planar
     2) the periodic mapping between them is either:
        a) translation along the normal between their
           parallel planes
        b) rotation about the line of intersection of
           their planes
     3) model edges are well represented by the point
        at the center of their parametric range
     note that these are quite restrictive and hard to
     check for violations, so use this code with great care. */
  apf::Frame frame;
  apf::Line axis;
  double angle;
  if (getAxisymmetry(gm, f, of, axis, angle)) {
    apf::Frame toline = apf::Frame::forTranslation(axis.origin * -1);
    apf::Frame rotation = apf::Frame::forRotation(axis.direction, angle);
    apf::Frame fromline = apf::Frame::forTranslation(axis.origin);
    frame = fromline * rotation * toline;
  } else {
    apf::Plane p = getFacePlane(gm, f);
    apf::Plane op = getFacePlane(gm, of);
    assert( ! apf::areClose(p, op, ph::tolerance));
    apf::Vector3 o = p.normal * p.radius;
    apf::Vector3 oo = op.normal * op.radius;
    frame = apf::Frame::forTranslation(oo - o);
  }
  closeFaceMatchingWithFrame(gm, f, of, mm, frame);
}

static void closeAttributeMatching(gmi_model* m, ModelMatching& mm)
{
  APF_ITERATE(ModelMatching, mm, it) {
    if (gmi_dim(m, it->first) == 2)
      closeFaceMatching(m, it->first, mm);
  }
}

void getFullAttributeMatching(gmi_model* m, BCs& bcs, ModelMatching& mm)
{
  getAttributeMatching(m, bcs, mm);
  closeAttributeMatching(m, mm);
  completeMatching(mm);
}

static void checkFilteredMatching(apf::Mesh* m, ModelMatching& mm, int dim)
{
  apf::MeshIterator* it = m->begin(dim);
  apf::MeshEntity* e;
  while ((e = m->iterate(it))) {
    apf::Matches matches;
    m->getMatches(e, matches);
    gmi_ent* ge = (gmi_ent*) m->toModel(e);
    if (!mm.count(ge)) {
      assert(matches.getSize() == 0);
      continue;
    }
    if (matches.getSize() < mm[ge].size()) {
      std::cerr << "solution periodicity requested "
                << "where mesh periodicity does not exist.\n"
                << "rebuild mesh to match solution request\n";
      abort();
    }
  }
  m->end(it);
}

void filterMatching(apf::Mesh2* m, ModelMatching& mm, int dim)
{
  gmi_model* gm;
  gm = m->getModel();
  PCU_Comm_Begin();
  apf::MeshIterator* it = m->begin(dim);
  apf::MeshEntity* e;
  int gd, gt;
  while ((e = m->iterate(it))) {
    apf::Matches matches;
    m->getMatches(e, matches);
    if (!matches.getSize())
      continue;
    gmi_ent* ge = (gmi_ent*) m->toModel(e);
    gd = gmi_dim(gm, ge);
    gt = gmi_tag(gm, ge);
    APF_ITERATE(apf::Matches, matches, mit) {
      PCU_COMM_PACK(mit->peer, mit->entity);
      PCU_COMM_PACK(mit->peer, e);
      PCU_COMM_PACK(mit->peer, gd);
      PCU_COMM_PACK(mit->peer, gt);
    }
    m->clearMatches(e);
  }
  m->end(it);
  PCU_Comm_Send();
  while (PCU_Comm_Receive()) {
    PCU_COMM_UNPACK(e);
    apf::MeshEntity* oe;
    PCU_COMM_UNPACK(oe);
    PCU_COMM_UNPACK(gd);
    PCU_COMM_UNPACK(gt);
    gmi_ent* ge = (gmi_ent*) m->toModel(e);
    if (!mm.count(ge))
      continue;
    ModelSet& ms = mm[ge];
    gmi_ent* oge = gmi_find(gm, gd, gt);
    if (oge == ge || ms.count(oge))
      m->addMatch(e, PCU_Comm_Sender(), oe);
  }
  checkFilteredMatching(m, mm, dim);
}

static SavedMatches* savedVertexMatches = 0;
static SavedMatches* savedFaceMatches = 0;

void enterFilteredMatching(apf::Mesh2* m, Input& in, BCs& bcs)
{
  if (!in.filterMatches)
    return;
  assert(PCU_Thrd_Peers() == 1);
  savedVertexMatches = new SavedMatches();
  saveMatches(m, 0, *savedVertexMatches);
  if (in.formElementGraph) {
    savedFaceMatches = new SavedMatches();
    saveMatches(m, 2, *savedFaceMatches);
  }
  ModelMatching mm;
  getFullAttributeMatching(m->getModel(), bcs, mm);
  filterMatching(m, mm, 0);
  if (in.formElementGraph)
    filterMatching(m, mm, 2);
}

void exitFilteredMatching(apf::Mesh2* m)
{
  if (savedVertexMatches)
    restoreMatches(m, 0, *savedVertexMatches);
  if (savedFaceMatches)
    restoreMatches(m, 2, *savedFaceMatches);
  delete savedVertexMatches;
  delete savedFaceMatches;
  savedVertexMatches = 0;
  savedFaceMatches = 0;
}

}
