#ifndef PH_INTERFACE_CUTTER_H
#define PH_INTERFACE_CUTTER_H

#include "phBC.h"
#include <apfMesh2.h>

namespace ph {

bool isInterface(gmi_model* gm, gmi_ent* ge, FieldBCs& fbcs);

void cutInterface(apf::Mesh2* m, BCs& bcs);

bool migrateInterface(apf::Mesh2*& m, ph::BCs& bcs);

}

#endif
