/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "maRegionCollapse.h"
#include "maAdapt.h"
#include "maShape.h"
#include <apfCavityOp.h>
#include <pcu_util.h>

namespace ma {

void RegionCollapse::Init(Adapt* a, double fa)
{
  adapt = a;
  region = 0;
  reclassifyEdge = 0;
  numBdryFaces = 0;
  faces[0] = faces[1] = faces[2] = faces[3] = 0;
  flatAngle = fa;
}

bool RegionCollapse::requestLocality(apf::CavityOp* o)
{
/* get vertices again since this is sometimes used
   before setVerts can work */
  Entity* v[4];
  Mesh* m = adapt->mesh;
  m->getDownward(region,0,v);
  return o->requestLocality(v,4);
}

bool RegionCollapse::setRegion(Entity* r)
{
  if (getFlag(adapt,r,DONT_COLLAPSE))
    return false;
  region = r;
  return true;
}

static int relativeDirection(Vector n1, Vector n2, double flatAngle)
{
  double cosFlatAngleSquare;
  double dot,ll1,ll2;

  cosFlatAngleSquare = cos(flatAngle*0.017453293);   //  PI/180
  cosFlatAngleSquare = cosFlatAngleSquare * cosFlatAngleSquare;

  dot = n1 * n2;

  ll1 = n1 * n1;
  ll2 = n2 * n2;

  if(dot*dot/(ll1*ll2) < cosFlatAngleSquare) return 0;
  if(dot > 0) return 1;
  return -1;
}

static Entity* getTetVertOppositeFace(Mesh* m, Entity* r, Entity* f)
{
  Entity* rfs[4];
  m->getDownward(r, 2, rfs);
  int index = findIn(rfs, 4, f);
  PCU_ALWAYS_ASSERT(index > -1);

  Entity* rvs[4];
  Entity* fvs[3];

  m->getDownward(r, 0, rvs);
  m->getDownward(f, 0, fvs);

  for (int i = 0; i < 4; i++) {
    int index = findIn(fvs, 3, rvs[i]);
    if (index == -1)
      return rvs[i];
  }
}

Vector faceNormal(Mesh* mesh, Entity* face, Entity* rgn)
{
  Vector xyz[4];
  Vector v01, v02, v03;
  Entity* fvs[3];
  mesh->getDownward(face, 0, fvs);
  for (int i = 0; i < 3; i++) {
    xyz[i] = getPosition(mesh, fvs[i]);
  }


  Entity* vertex = getTetVertOppositeFace(mesh, rgn, face);
  xyz[3] = getPosition(mesh, vertex);

  v01 = xyz[1] - xyz[0];
  v02 = xyz[2] - xyz[0];
  v03 = xyz[3] - xyz[0];

  Vector normal = cross(v01, v02);
  if( (normal * v03) > 0 )
    normal = normal * (-1);

  return normal;
}

bool RegionCollapse::checkGeom()
{
  Vector normal1, normal2;
  Mesh* mesh = adapt->mesh;
  switch (numBdryFaces) {
    case 1:
      normal1 = faceNormal(mesh, faces[0], region);
      for (int i = 1; i < 4; i++) {
      	normal2 = faceNormal(mesh, faces[i], region);
      	if (relativeDirection(normal1, normal2, flatAngle) != -1)
      	  return false;
      }
      break;
    case 2:
      normal1 = faceNormal(mesh, faces[0], region);
      normal2 = faceNormal(mesh, faces[1], region);
      if (relativeDirection(normal1, normal2, flatAngle) != 1)
      	return false;
      normal1 = faceNormal(mesh, faces[2], region);
      normal2 = faceNormal(mesh, faces[3], region);
      if (relativeDirection(normal1, normal2, flatAngle) != 1)
      	return false;
      break;
    case 3:
      normal1 = faceNormal(mesh, faces[3], region);
      for (int i = 0; i < 3; i++) {
	normal2 = faceNormal(mesh, faces[i], region);
	if (!relativeDirection(normal1, normal2, flatAngle))
	  return false;
      }
      break;
    default:
      return false;
  }
  return true;
}


bool RegionCollapse::checkTopo()
{

  Mesh* mesh = adapt->mesh;

  Model* modelFace = 0;
  Model* uniqueModelFace = 0;
  numBdryFaces = 0;

  Entity* fs[4];
  mesh->getDownward(region, 2, fs);
  int m = 0;
  for (int i = 0; i < 4; i++) {
    Entity* face = fs[i];

    int faceModelType = mesh->getModelType(mesh->toModel(face));
    if (faceModelType != 2) {
      faces[3-m] = face;
      ++m;
      continue;
    }

    if (mesh->countUpward(face) == 2) {
      faces[3-m] = face;
      ++m;
      continue;
    }

    faces[numBdryFaces] = face;
    ++numBdryFaces;

    modelFace = mesh->toModel(face);
    if (uniqueModelFace) {
      if (uniqueModelFace != modelFace)
      	return false;
    }
    else
      uniqueModelFace = modelFace;
  }

  Entity* oppositeEdge = 0;
  switch (numBdryFaces) {
    case 0:
      return false;
    case 1:
      {
	Entity* vert = getTetVertOppositeFace(mesh, region, faces[0]);
	if (mesh->getModelType(mesh->toModel(vert)) != 3)
	  return false;
	break;
      }
    case 2:
      {
	Entity* f2es[3];
	Entity* f3es[3];
	mesh->getDownward(faces[2], 1, f2es);
	mesh->getDownward(faces[3], 1, f3es);
	int j;
	for (j = 0; j < 3; j++) {
	  int index = findIn(f3es, 3, f2es[j]);
	  if (index > -1)
	    break;
	}
	reclassifyEdge = f2es[j]; // this is the edge being reclassified

	Entity* f0es[0];
	Entity* f1es[1];
	mesh->getDownward(faces[0], 1, f0es);
	mesh->getDownward(faces[1], 1, f1es);
	for (j = 0; j < 3; j++) {
	  int index = findIn(f1es, 3, f0es[j]);
	  if (index > -1)
	    break;
	}
	oppositeEdge = f0es[j]; // this is the edge being deleted

	if (mesh->getModelType(mesh->toModel(reclassifyEdge)) != 3)
	  return false;
	if (mesh->getModelType(mesh->toModel(oppositeEdge)) == 1)
	  return false;
	break;
      }
    case 3:
      {
	Entity* f3es[3];
	for (int j = 0; j < 2; j++) {
	  Entity* fes[3];
	  mesh->getDownward(faces[j], 1, fes);
	  for (int k = 0; k < 3; k++) {
	    int index = findIn(f3es, 3, fes[k]);
	    if (index > -1) continue;
	    if (mesh->getModelType(mesh->toModel(fes[k])) != 2)
	      return false;
	  }
	}
	Entity* vert = getTetVertOppositeFace(mesh, region, vert);
	if (mesh->getModelType(mesh->toModel(vert)) != 2)
	  return false;
	if (mesh->countUpward(vert) != 3)
	  return false;
	break;
      }
    default:
      return false;
  }

  return true;
}

static void processNewBdryVert(Mesh* m, Entity* v)
{
  Model* c = m->toModel(v);
  int modelType = m->getModelType(c);
  PCU_ALWAYS_ASSERT(modelType == 1 || modelType == 2);

  Vector coords = getPosition(m, v);
  Vector newCoords;
  Vector newParams;
  m->getClosestPoint(c, coords, newCoords, newParams);
  m->setParam(v, newParams);
  // TODO: this must be done via snapping
  /* m->setPoint(v, 0, newCoords); */
}

void RegionCollapse::apply()
{
  Mesh* mesh = adapt->mesh;
  Model* modelFace = mesh->toModel(faces[0]);

  switch (numBdryFaces) {
    case 1:
      {
	// reclassify faces
	for (int i = 1; i < 4; i++) {
	  mesh->setModelEntity(faces[i], modelFace);
	  // TODO: do we need to change direction of faces?
	}
	// reclassify edges
	Entity* f0es[3];
	mesh->getDownward(faces[0], 1, f0es);
	Entity* fes[3];
	mesh->getDownward(faces[1], 1, fes);
	for (int i = 0; i < 3; i++) {
	  Entity* edge = fes[i];
	  int index = findIn(f0es, 3, edge);
	  if (index > -1) continue;
	  mesh->setModelEntity(edge, modelFace);
	}
	mesh->getDownward(faces[2], 1, fes);
	for (int i = 0; i < 3; i++) {
	  Entity* edge = fes[i];
	  int index = findIn(f0es, 3, edge);
	  if (index > -1) continue;
	  mesh->setModelEntity(edge, modelFace);
	}
	Entity* vert = getTetVertOppositeFace(mesh, region, faces[0]);
	mesh->setModelEntity(vert, modelFace);
	processNewBdryVert(mesh, vert);

	mesh->destroy(region);
	mesh->destroy(faces[0]);
	break;
      }
    case 2:
      {
	Entity* edgeToRemove;
	Entity* f0es[3];
	Entity* f1es[3];
	mesh->getDownward(faces[0], 1, f0es);
	mesh->getDownward(faces[1], 1, f1es);
	for (int i = 0; i < 3; i++) {
	  edgeToRemove = f0es[i];
	  int index = findIn(f1es, 3, edgeToRemove);
	  if (index > -1) break;
	}

	// reclassify faces
	for (int i = 2; i < 4; i++) {
	  mesh->setModelEntity(faces[i], modelFace);
	  // TODO: do we need to change direction of faces?
	}

	mesh->setModelEntity(reclassifyEdge, modelFace);

	mesh->destroy(region);
	mesh->destroy(faces[0]);
	mesh->destroy(faces[1]);
	mesh->destroy(edgeToRemove);
	break;
      }
    case 3:
      {
	mesh->setModelEntity(faces[3], modelFace);
	// TODO: do we need to change direction of faces?

	Entity* vert = getTetVertOppositeFace(mesh, region, faces[3]);
	int n = mesh->countUpward(vert);
	PCU_ALWAYS_ASSERT(n == 3);
	Entity* edges[3];
	for (int i = 0; i < n; i++) {
	  edges[i] = mesh->getUpward(vert, i);
	}


	mesh->destroy(region);
	mesh->destroy(faces[0]);
	mesh->destroy(faces[1]);
	mesh->destroy(faces[2]);
	mesh->destroy(edges[0]);
	mesh->destroy(edges[1]);
	mesh->destroy(edges[2]);
	mesh->destroy(vert);
	break;
      }
  }
}

void RegionCollapse::unmark()
{
  clearFlag(adapt,region,COLLAPSE);
}


bool setupRegionCollapse(RegionCollapse& rcollapse, Entity* region)
{
  Adapt* adapter = rcollapse.adapt;
  PCU_ALWAYS_ASSERT(adapter->mesh->getType(region) == apf::Mesh::TET);
  if ( ! rcollapse.setRegion(region))
    return false;
  if ( ! rcollapse.checkGeom())
    return false;
  if ( ! rcollapse.checkTopo())
    return false;

  return true;
}

}
