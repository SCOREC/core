/******************************************************************************

  Copyright 2013 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.

*******************************************************************************/
#include "maMesh.h"
#include "maSnap.h"
#include "maFaceSplit.h"
#include "maSolutionTransfer.h"
#include "maShapeHandler.h"

namespace ma {

static int getTetVertIdOppositeTri(Mesh* m, Entity* tet, Entity* tri)
{
  Entity* tetv[4];
  m->getDownward(tet,0,tetv);
  Entity* triv[3];
  m->getDownward(tri,0,triv);
  for (int i=0; i < 4; ++i)
    if (-1==findIn(triv,3,tetv[i]))
      return i;
  return -1;
}

/** \brief Rotates the tetrahedron to a standard orientation for splitting
 * \details Specifically, such that 4th vertex is opposite the split face
 * \param a The Adapt instance
 * \param tet The tetrahedron
 * \param tri The face to be split
 * \param v A Downward
 */
static void rotateForFaceSplit(Adapt* a, Entity* tet, Entity* tri, Downward v)
{
  Mesh* m = a->mesh;
  int type = m->getType(tet);
  int rots[4] = {4, 1, 2, 0};
  Downward vi;
  m->getDownward(tet, 0, vi);
  rotateEntity(type, vi, rots[getTetVertIdOppositeTri(m, tet, tri)], v);
}

FaceSplit::FaceSplit(Adapt* a):
  adapter(a)
{
}

bool FaceSplit::setFace(Entity* face)
{
  PCU_ALWAYS_ASSERT_VERBOSE(adapter->mesh->getType(face) == apf::Mesh::TRIANGLE,
			    "Only simplicial meshes are supported.\n");
  if (getFlag(adapter, face, DONT_SPLIT))
    return false;
  for (int d=1; d <= 3; ++d)
    toSplit[d].setSize(0);
  // Add face to toSplit.
  // No flags are set because we know exactly what is to be split
  toSplit[2].setSize(1);
  toSplit[2][0] = face;
  // Add up regions to toSplit.
  // No flags are set because we know exactly what is to be split
  apf::Up rs;
  adapter->mesh->getUp(face, rs);
  toSplit[3].setSize(rs.n);
  for (int i = 0; i < rs.n; ++i) {
    toSplit[3][i] = rs.e[i];
  }
  return true;
}

void FaceSplit::makeNewElements()
{
  Mesh* m = adapter->mesh;
  NewEntities cb;
  PCU_ALWAYS_ASSERT(toSplit[2].getSize() == 1);
  // Doing the BuildCallback routine (similar to ma::splitElements)
  // only for the split face
  newEntities[2].setSize(toSplit[2].getSize());
  setBuildCallback(adapter, &cb);
  cb.reset();
  Entity* tri = toSplit[2][0];
  // Split vertex. Also makes the new faces.
  // TODO: Add all these to newEntities, or do something similar
  Entity* sv = splitTri0(adapter, tri);
  cb.retrieve(newEntities[2][0]);
  clearBuildCallback(adapter);
  Downward tetv, ntetv;
  for (size_t i = 0; i < toSplit[3].getSize(); ++i) {
    Entity* tet = toSplit[3][i];
    // Orient Tets
    rotateForFaceSplit(adapter, tet, tri, tetv);
    // Create splits
    ntetv[0] = tetv[0]; ntetv[1] = tetv[1]; ntetv[2] = sv; ntetv[3] = tetv[3];
    buildElement(adapter, m->toModel(tet), apf::Mesh::TET, ntetv);
    ntetv[0] = tetv[1]; ntetv[1] = tetv[2]; ntetv[2] = sv; ntetv[3] = tetv[3];
    buildElement(adapter, m->toModel(tet), apf::Mesh::TET, ntetv);
    ntetv[0] = tetv[2]; ntetv[1] = tetv[0]; ntetv[2] = sv; ntetv[3] = tetv[3];
    buildElement(adapter, m->toModel(tet), apf::Mesh::TET, ntetv);
  }
}

void FaceSplit::cancel()
{
  Mesh* m = adapter->mesh;
  // TODO Find Split vert and delete all adjacent regions
  Entity* sv = getSplitVert();
  Upward deleting; // We'll delete all regions adjacent to sv
  m->getAdjacent(sv, m->getDimension(), deleting);
  APF_ITERATE(Upward, deleting, it)
    destroyElement(adapter, *it);
  // TODO Any cleanup?
  // TODO This code is based on the Splits operator, and there,
  // some more work is done, relating to removal of tags. Should investigate
  // if that is needed here too.
}

void FaceSplit::transfer()
{
  Mesh* m = adapter->mesh;
  // TODO
  SolutionTransfer* st = adapter->solutionTransfer;
  int td = st->getTransferDimension();
  for (int d = td; d <= m->getDimension(); ++d)
    for (size_t i = 0; i < toSplit[d].getSize(); ++i)
      st->onRefine(toSplit[d][i], newEntities[d][i]);

  td = adapter->shape->getTransferDimension();
  for (int d = td; d <= m->getDimension(); ++d)
    for (size_t i = 0; i < toSplit[d].getSize(); ++i)
      adapter->shape->onRefine(toSplit[d][i], newEntities[d][i]);
}

void FaceSplit::destroyOldElements()
{
  Mesh* m = adapter->mesh;
  // Destroy split elements
  int D = m->getDimension();
  for (size_t i = 0; i < toSplit[D].getSize(); ++i)
    destroyElement(adapter, toSplit[D][i]);
  // Clear everything
  for (int d = 2; d <= D; ++d){
    toSplit[d].setSize(0);
    newEntities[d].setSize(0);
  }
}

Entity* FaceSplit::getSplitVert()
{
  EntityArray& a = newEntities[2][0];
  for (size_t i=0; i < a.getSize(); ++i) {
    if (adapter->mesh->getType(a[i]) == apf::Mesh::VERTEX)
      return a[i];
  }
  return 0;
}

Entity* makeSplitVertOnFace(Adapt* a, Entity* face)
{
  Mesh* m = a->mesh;
  Model* c = m->toModel(face);
  SizeField* sf = a->sizeField;
  SolutionTransfer* st = a->solutionTransfer;
  // midpoint of face (parametric space)
  Vector xi(1.0/3.0, 1.0/3.0, 1.0/3.0);
  apf::MeshElement* me = apf::createMeshElement(m, face);
  Vector point;
  apf::mapLocalToGlobal(me,xi,point);
  Vector param(0,0,0);
  if (a->input->shouldTransferParametric)
    transferParametricOnTriSplit(m, face, xi, param);
  Entity* vert = buildVertex(a, c, point, param);
  st->onVertex(me, xi, vert);
  sf->interpolate(me, xi, vert);
  apf::destroyMeshElement(me);
  return vert;
}

Entity* splitTri0(Adapt* a, Entity* parent)
{
  Entity* sv = makeSplitVertOnFace(a, parent);
  Downward v;
  a->mesh->getDownward(parent, 0, v);
  Entity* tv[3];
  tv[0] = v[0]; tv[1] = v[1]; tv[2] = sv;
  buildElement(a, a->mesh->toModel(parent), apf::Mesh::TRIANGLE, tv);
  tv[0] = v[1]; tv[1] = v[2]; tv[2] = sv;
  buildElement(a, a->mesh->toModel(parent), apf::Mesh::TRIANGLE, tv);
  tv[0] = v[2]; tv[0] = v[0]; tv[2] = sv;
  buildElement(a, a->mesh->toModel(parent), apf::Mesh::TRIANGLE, tv);
  return sv;
}

}
