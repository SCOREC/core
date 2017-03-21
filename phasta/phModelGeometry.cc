#include "phModelGeometry.h"
#include <apf.h>
#include <pcu_util.h>

namespace ph {

double const tolerance = 1e-9;

static apf::Vector3 getVertexCenter(gmi_model* gm, gmi_ent* v)
{
  double p[2] = {0,0};
  apf::Vector3 center;
  gmi_eval(gm, v, p, &center[0]);
  return center;
}

static apf::Vector3 getEdgeCenter(gmi_model* gm, gmi_ent* e)
{
  double r[2];
  gmi_range(gm, e, 0, r);
  double p[2];
  p[0] = (r[1] + r[0]) / 2;
  apf::Vector3 center;
  gmi_eval(gm, e, p, &center[0]);
  return center;
}

apf::Vector3 getCenter(gmi_model* gm, gmi_ent* e)
{
  switch (gmi_dim(gm, e)) {
    case 0:
      return getVertexCenter(gm, e);
    case 1:
      return getEdgeCenter(gm, e);
  };
  apf::fail("ph::getCenter called on something not a vertex or edge");
}

apf::Plane getFacePlane(gmi_model* gm, gmi_ent* f)
{
  double r[2];
  double p[2];
  gmi_range(gm, f, 0, r);
  p[0] = r[0];
  gmi_range(gm, f, 1, r);
  p[1] = r[0];
  apf::Vector3 origin;
  gmi_eval(gm, f, p, &origin[0]);
  apf::Vector3 normal;
  gmi_normal(gm, f, p, &normal[0]);
  return apf::Plane(normal, normal * origin);
}

apf::Vector3 getAnyPointOnFace(gmi_model* gm, gmi_ent* f)
{
  gmi_set* s = gmi_adjacent(gm, f, 1);
  PCU_ALWAYS_ASSERT(s);
  PCU_ALWAYS_ASSERT(s->n >= 1);
  apf::Vector3 p = getEdgeCenter(gm, s->e[0]);
  gmi_free_set(s);
  return p;
}

}
