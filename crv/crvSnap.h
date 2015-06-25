/******************************************************************************

  Copyright 2013 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.

*******************************************************************************/

#ifndef CRVSNAP_H
#define CRVSNAP_H

#include "apfMesh.h"
#include "apfMesh2.h"
#include "apfShape.h"

namespace crv {

void transferParametricOnTriSplit(
    apf::Mesh2* m,
    apf::MeshEntity* e,
    apf::Vector3& t,
    apf::Vector3& p);
void transferParametricOnGeometricEdgeSplit(
    apf::Mesh2* m,
    apf::MeshEntity* e,
    double t,
    apf::Vector3& p);
void transferParametricOnGeometricTriSplit(
    apf::Mesh2* m,
    apf::MeshEntity* e,
    apf::Vector3& t,
    apf::Vector3& p);

}

#endif
