#ifndef PH_AXISYMMETRY_H
#define PH_AXISYMMETRY_H

#include "phModelGeometry.h"
#include "phBC.h"
#include <apfMesh.h>

namespace ph {

bool getAxisymmetry(gmi_model* gm, gmi_ent* f, gmi_ent* of,
    apf::Line& axis, double& angle);
void attachAllAngleBCs(gmi_model* gm, BCs& bcs);
apf::MeshTag* tagAngles(apf::Mesh* m, BCs& bcs, apf::MatchedSharing* ms);

}

#endif
