#include "phAxisymmetry.h"
#include <apf.h>
#include <PCU.h>
#include <cassert>
#include <iostream>

namespace ph {

bool getAxisymmetry(gmi_model* gm, gmi_ent* f, gmi_ent* of,
    apf::Line& axis, double& angle)
{
  apf::Plane p = getFacePlane(gm, f);
  apf::Plane op = getFacePlane(gm, of);
  if (apf::areParallel(p, op, ph::tolerance))
    return false;
  axis = apf::intersect(p, op);
  /* we need a couple points to know which of the four
     quadrants between the planes we are dealing with */
  apf::Vector3 pt = getAnyPointOnFace(gm, f);
  apf::Vector3 opt = getAnyPointOnFace(gm, of);
  apf::Vector3 v = apf::reject(pt - axis.origin, axis.direction);
  apf::Vector3 ov = apf::reject(opt - axis.origin, axis.direction);
  angle = apf::getAngle(v, ov);
  if (apf::cross(v, ov) * axis.direction < 0)
    angle = -angle;
  return true;
}

static void attachAngleBC(BCs& bcs, gmi_model* gm, gmi_ent* f, double angle)
{
  ConstantBC* bc = makeConstantBC(bcs, "ph::angle", gmi_dim(gm, f),
      gmi_tag(gm, f), 1);
  bc->value[0] = angle;
}

static void tryAttachingAngleBCs(BCs& bcs, gmi_model* gm, gmi_ent* f, gmi_ent* of)
{
  apf::Line axis;
  double angle;
  if (!getAxisymmetry(gm, f, of, axis, angle))
    return;
  /* PHASTA is constrained to axisymmetry around the Z axis */
  assert(apf::areClose(axis.origin, apf::Vector3(0,0,0), ph::tolerance));
  assert(apf::areParallel(axis.direction, apf::Vector3(0,0,1), ph::tolerance));
  if (axis.direction.z() < 0)
    angle = -angle;
  attachAngleBC(bcs, gm, f, angle);
  attachAngleBC(bcs, gm, of, -angle);
}

void attachAllAngleBCs(gmi_model* gm, BCs& bcs)
{
  std::string name = "periodic slave";
  if (!haveBC(bcs, name))
    return;
  FieldBCs& pfbcs = bcs.fields[name];
  FieldBCs& afbcs = bcs.fields["ph::angle"];
  APF_ITERATE(FieldBCs::Set, pfbcs.bcs, it) {
    BC* bc = *it;
    gmi_ent* e = gmi_find(gm, bc->dim, bc->tag);
    double* val = bc->eval(apf::Vector3(0,0,0));
    int otherTag = *val;
    gmi_ent* oe = gmi_find(gm, bc->dim, otherTag);
    if (getBCValue(gm, afbcs, e)) {
      assert(getBCValue(gm, afbcs, oe));
      continue;
    }
    tryAttachingAngleBCs(bcs, gm, e, oe);
  }
}

static double* getAngleBC(gmi_model* gm, BCs& bcs, gmi_ent* e)
{
  int d = gmi_dim(gm, e);
  if (d > 2)
    return 0;
  if (d == 2)
    return getBCValue(gm, bcs.fields["ph::angle"], e);
  gmi_set* s = gmi_adjacent(gm, e, d + 1);
  double* angle = 0;
  for (int i = 0; i < s->n; ++i) {
    angle = getAngleBC(gm, bcs, s->e[i]);
    if (angle)
      break;
  }
  gmi_free_set(s);
  return angle;
}

static bool requiresRotation(gmi_model* gm, BCs& bcs, gmi_ent* e, gmi_ent* oe,
    double& angle)
{
  double a;
  double oa;
  double* ap;
  ap = getAngleBC(gm, bcs, e);
  if (!ap)
    return false;
  a = *ap;
  ap = getAngleBC(gm, bcs, oe);
  if (!ap)
    return false;
  oa = *ap;
  if (a * oa > 0)
    return false;
  angle = a;
  return true;
}

apf::MeshTag* tagAngles(apf::Mesh* m, BCs& bcs, apf::MatchedSharing* ms)
{
  apf::MeshTag* tag = m->createDoubleTag("ph_angle", 1);
  gmi_model* gm = m->getModel();
  PCU_Comm_Begin();
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* v;
  while ((v = m->iterate(it))) {
    apf::Matches matches;
    m->getMatches(v, matches);
    if (!matches.getSize())
      continue;
    if (!ms->isOwned(v))
      continue;
    apf::ModelEntity* me = m->toModel(v);
    int mdim, mtag;
    mdim = m->getModelType(me);
    mtag = m->getModelTag(me);
    APF_ITERATE(apf::Matches, matches, mit) {
      PCU_COMM_PACK(mit->peer, mit->entity);
      PCU_COMM_PACK(mit->peer, mdim);
      PCU_COMM_PACK(mit->peer, mtag);
    }
  }
  m->end(it);
  PCU_Comm_Send();
  while (PCU_Comm_Receive()) {
    PCU_COMM_UNPACK(v);
    int mdim, mtag;
    PCU_COMM_UNPACK(mdim);
    PCU_COMM_UNPACK(mtag);
    gmi_ent* oge = gmi_find(gm, mdim, mtag);
    gmi_ent* ge = (gmi_ent*) m->toModel(v);
    double angle;
    if (requiresRotation(gm, bcs, ge, oge, angle))
      m->setDoubleTag(v, tag, &angle);
  }
  return tag;
}

}
