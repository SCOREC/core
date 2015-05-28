#ifndef PH_INTERFACE_CUTTER_H
#define PH_INTERFACE_CUTTER_H

#include "phBC.h"
#include <apfMesh2.h>

namespace ph {

void cutInterface(apf::Mesh2* m, BCs& bcs);

}

#endif
