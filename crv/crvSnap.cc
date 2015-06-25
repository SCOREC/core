/******************************************************************************

  Copyright 2013 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.

 *******************************************************************************/

#include "crvSnap.h"

namespace crv {

void transferParametricOnGeometricEdgeSplit(
    apf::Mesh2* m,
    apf::MeshEntity* e,
    double t,
    apf::Vector3& p)
{
  apf::ModelEntity* g = m->toModel(e);
  int modelDimension = m->getModelType(g);
  if (modelDimension==m->getDimension()) return;
  apf::MeshEntity* ev[2];
  m->getDownward(e,0,ev);
  apf::Vector3 p0,p1,cpt;
  m->getPoint(ev[0],0,p0);
  m->getPoint(ev[1],0,p1);
  apf::Vector3 pt = p0*(1.-t)+p1*t;
  m->getClosestPoint(g,pt,cpt,p);
}

void transferParametricOnGeometricTriSplit(
    apf::Mesh2* m,
    apf::MeshEntity* e,
    apf::Vector3& t,
    apf::Vector3& p)
{
  apf::ModelEntity* g = m->toModel(e);
  int modelDimension = m->getModelType(g);
  if (modelDimension==m->getDimension()) return;
  apf::MeshEntity* ev[3];
  m->getDownward(e,0,ev);
  // split in physical space, project
  apf::Vector3 p0,p1,p2,cpt;
  m->getPoint(ev[0],0,p0);
  m->getPoint(ev[1],0,p1);
  m->getPoint(ev[2],0,p2);
  apf::Vector3 pt = p0*(1.-t[0]-t[1])+p1*t[0]+p2*t[1];
  m->getClosestPoint(g,pt,cpt,p);
}


} // namespace crv

