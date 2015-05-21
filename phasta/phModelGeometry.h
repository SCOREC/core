#ifndef PH_MODEL_GEOMETRY_H
#define PH_MODEL_GEOMETRY_H

#include <gmi.h>
#include <apfGeometry.h>

namespace ph {

extern double const tolerance;

apf::Vector3 getCenter(gmi_model* gm, gmi_ent* e);
apf::Plane getFacePlane(gmi_model* gm, gmi_ent* f);
apf::Vector3 getAnyPointOnFace(gmi_model* gm, gmi_ent* f);

}

#endif
