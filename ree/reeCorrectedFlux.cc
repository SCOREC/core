/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */
#include <apfCavityOp.h>
#include "apfElement.h"
#include "crv.h"
#include "crvShape.h"
#include "ree.h"

namespace ree {

enum {VISITED};

/* overall information useful during computation of corrected flux */
struct CorrectFlux
{
  apf::Mesh* mesh;
  /* mesh dimension, so far handling 3 only */
  int dim;
  /* polynomial order of nedelec space */
  int order;
  /* polynomial order of p+1 nedelec space */
  int orderp1;
  /* input nedelec field containing scalar dofs for electric field */
  apf::Field* ef;
  /* input theta field for flux correction
   * (3 scalar dofs on each face, 1 per edge) */
  apf::Field* theta;
  /* tags each face once it has been visited */
  apf::MeshTag* tag;
  /* output per-element field containing correctedFlux vectors,
   * which are computed using theta and electric fields */
  apf::Field* correctedFlux;
};

static void setupCorrectFlux(
    CorrectFlux* cf,
    apf::Field* f,
    apf::Field* theta,
    apf::Field* correctedFlux)
{
  cf->mesh = apf::getMesh(f);
  cf->dim = cf->mesh->getDimension();
  cf->order = f->getShape()->getOrder();
  cf->orderp1 = cf->order+1;
  cf->ef = f;
  cf->theta = theta;
  cf->tag = cf->mesh->createIntTag("isVisited", 1);

  cf->correctedFlux = correctedFlux;
  int nc = apf::countComponents(cf->correctedFlux);
  double zeros[nc];
  for (int i = 0; i < nc; i++)
    zeros[i] = 0.;
  apf::MeshEntity* tet;
  apf::MeshIterator* it = cf->mesh->begin(3);
  while ((tet = cf->mesh->iterate(it))) {
    apf::MeshElement* me = apf::createMeshElement(cf->mesh, tet);
    int orderp1 = cf->order+1;
    int np = apf::countIntPoints(me, 2*orderp1-1);

    for (int i = 0; i < np; i++) {
      apf::setComponents(cf->correctedFlux, tet, i, zeros);
    }
    apf::destroyMeshElement(me);
  }
  cf->mesh->end(it);
}

typedef std::vector<apf::MeshEntity*> EntityVector;
struct FaceCavity
{
  apf::Mesh* mesh;
  apf::MeshEntity* entity;
  CorrectFlux* correctflux;
  EntityVector tets;
};

static void setupFaceCavity(FaceCavity* fc, CorrectFlux* cf)
{
  fc->mesh = cf->mesh;
  fc->correctflux = cf;
  fc->entity = 0;
}

static void startFaceCavity(FaceCavity* fc, apf::MeshEntity* f)
{
  fc->entity = f;
  fc->tets.clear();
}

static void addEntityToCavity(FaceCavity* fc, apf::MeshEntity* e)
{
  PCU_ALWAYS_ASSERT(fc->mesh->getType(e) == apf::Mesh::TET);
  fc->tets.push_back(e);
}

static void addEntitiesToCavity(
    FaceCavity* fc, apf::DynamicArray<apf::MeshEntity*>& es)
{
  for (std::size_t i=0; i < es.getSize(); ++i)
    addEntityToCavity(fc, es[i]);
}

static bool getInitialFaceCavity(FaceCavity* fc, apf::CavityOp* o)
{
  if (! o->requestLocality(&fc->entity, 1))
    return false;

  apf::DynamicArray<apf::MeshEntity*> adjacent;
  fc->mesh->getAdjacent(fc->entity, 3, adjacent);
  addEntitiesToCavity(fc, adjacent);
  return true;
}

static bool buildFaceCavity(FaceCavity* fc, apf::CavityOp* o)
{
  if (!getInitialFaceCavity(fc, o)) return false;
  return true;
}


static void computeCorrectedFlux(FaceCavity* fc)
{
  apf::MeshEntity* face = fc->entity;

  // 1. get upward tets of the face
  apf::Up up;
  fc->mesh->getUp(face, up);
  if (crv::isBoundaryEntity(fc->mesh, face))
    PCU_ALWAYS_ASSERT(up.n == 1);
  else
    PCU_ALWAYS_ASSERT(up.n == 2);

  apf::MeshEntity* firstTet = up.e[0];
  apf::MeshEntity* secondTet = nullptr;
  if (up.n == 2)
    secondTet  = up.e[1];

  // 2. get positions of the face in downward faces of upward tets
  apf::Downward tet1_faces, tet2_faces;
  int nf = fc->mesh->getDownward(firstTet, 2, tet1_faces);
  if (up.n == 2)
    fc->mesh->getDownward(secondTet, 2, tet2_faces);
  int tet1_pos = -1;
  int tet2_pos = -1;
  tet1_pos = apf::findIn(tet1_faces, nf, face);
  if (up.n == 2)
    tet2_pos = apf::findIn(tet2_faces, nf, face);

  // 3. get downward edges of the face
  apf::Downward edges;
  int nedges = fc->mesh->getDownward(face, 1, edges);

  // 4. get theta coeffs on the face
  double components[3];
  apf::getComponents(fc->correctflux->theta, face, 0, components);
  mth::Vector<double> theta_coeffs(3);
  theta_coeffs(0) = components[0];
  theta_coeffs(1) = components[1];
  theta_coeffs(2) = components[2];


  // 5. Evaluate and save corrected flux vector in an auxiliary field
  int ftype = fc->mesh->getType(face);
  PCU_ALWAYS_ASSERT(ftype == apf::Mesh::TRIANGLE);
  int nfdofs = apf::countElementNodes(fc->correctflux->ef->getShape(), ftype);
  apf::NewArray<apf::Vector3> vectorshapes(nfdofs);

  apf::MeshElement* fme = apf::createMeshElement(fc->mesh, face);
  apf::Element* fel = apf::createElement(fc->correctflux->ef, fme);
  int int_order = 2*fc->correctflux->orderp1-1;
  int np = apf::countIntPoints(fme, int_order);
  int nc = apf::countComponents(fc->correctflux->correctedFlux);

  apf::Vector3 p, tet1xi, tet2xi, curl1, curl2, curl,
      fnormal1, fnormal2, tk, tk1, tk2, vshape;
  for (int n = 0; n < np; n++) {

    apf::getIntPoint(fme, int_order, n, p);

    // evaluate theta vector using theta coeffs
    apf::Vector3 theta_vector, theta_vector1, theta_vector2;
    theta_vector.zero(); theta_vector1.zero(); theta_vector2.zero();
    apf::NewArray<apf::Vector3> triVectorShapes (nfdofs);
    apf::getVectorShapeValues(fel, p, triVectorShapes);
    for (int i = 0; i < nedges; i++) {
      apf::Vector3 v = triVectorShapes[i];
      v = v * theta_coeffs[i];
      theta_vector += v;
    }

    // orient face outward normals wrt tets and theta vectors
    fnormal1 = computeFaceOutwardNormal(fc->mesh, firstTet, face, p);
    fnormal2 = apf::Vector3(0.,0.,0.);
    theta_vector1 = theta_vector;
    if (up.n == 2) {
      fnormal2 = computeFaceOutwardNormal(fc->mesh, secondTet, face, p);
      theta_vector2 = theta_vector * -1.;
    }

    curl.zero();
    // compute curl1
    tet1xi = apf::boundaryToElementXi(fc->mesh, face, firstTet, p);
    apf::MeshElement* me1 = apf::createMeshElement(fc->mesh, firstTet);
    apf::Element* el1 = apf::createElement(fc->correctflux->ef, me1);
    apf::getCurl(el1, tet1xi, curl1);
    apf::Vector3 temp1 = apf::cross(fnormal1, curl1);
    curl += temp1;
    apf::destroyElement(el1);
    apf::destroyMeshElement(me1);

    // compute curl2
    if (up.n == 2) {
      tet2xi = apf::boundaryToElementXi(fc->mesh, face, secondTet, p);
      apf::MeshElement* me2 = apf::createMeshElement(fc->mesh, secondTet);
      apf::Element* el2 = apf::createElement(fc->correctflux->ef, me2);
      apf::getCurl(el2, tet2xi, curl2);
      apf::Vector3 temp2 = apf::cross(fnormal2, curl2);
      curl += (temp2 * -1.);
      curl = curl * 1./2.;
      apf::destroyElement(el2);
      apf::destroyMeshElement(me2);
    }

    tk = curl;
    tk1 = tk;
    apf::Vector3 theta_plus_tk1 = theta_vector1 + tk1;

    // get and set components in the auxiliary field
    double comp1[nc];
    apf::getComponents(fc->correctflux->correctedFlux, firstTet, n, comp1);
    int id1 = tet1_pos * 3;
    comp1[id1] = theta_plus_tk1[0];
    comp1[id1+1] = theta_plus_tk1[1];
    comp1[id1+2] = theta_plus_tk1[2];
    apf::setComponents(fc->correctflux->correctedFlux, firstTet, n, comp1);

    if (up.n == 2) {
      tk2 = tk * -1.;
      apf::Vector3 theta_plus_tk2 = theta_vector2 + tk2;
      double comp2[nc];
      apf::getComponents(fc->correctflux->correctedFlux, secondTet, n, comp2);
      int id2 = tet2_pos * 3;
      comp2[id2] = theta_plus_tk2[0];
      comp2[id2+1] = theta_plus_tk2[1];
      comp2[id2+2] = theta_plus_tk2[2];
      apf::setComponents(fc->correctflux->correctedFlux, secondTet, n, comp2);
    }
  }
  apf::destroyElement(fel);
  apf::destroyMeshElement(fme);
}

class FaceCavityOp : public apf::CavityOp
{
public:
  FaceCavityOp(CorrectFlux* cf):
    apf::CavityOp(cf->mesh)
  {
    setupFaceCavity(&face_cavity, cf);
  }
  virtual Outcome setEntity(apf::MeshEntity* e)
  {
    if (face_cavity.mesh->hasTag(e, face_cavity.correctflux->tag))
      return SKIP;
    startFaceCavity(&face_cavity, e);
    if ( ! buildFaceCavity(&face_cavity, this)) {
      return REQUEST;
    }
    return OK;
  }
  virtual void apply()
  {
    computeCorrectedFlux(&face_cavity);
    int n = VISITED;
    face_cavity.mesh->setIntTag(
        face_cavity.entity, face_cavity.correctflux->tag, &n);
  }
  FaceCavity face_cavity;
};

apf::Field* computeCorrectedFlux(apf::Field* ef, apf::Field* theta)
{
  int dim = apf::getMesh(ef)->getDimension();
  int order = ef->getShape()->getOrder() + 1; // local BVPs require p+1
  int int_order = 2*order-1;
  int nc = 4*3; // 1 flux vector per face
  apf::Field* correctedFlux = createPackedField(
      apf::getMesh(ef), "correctedFlux", nc, apf::getIPShape(dim, int_order));

  CorrectFlux correctflux;
  setupCorrectFlux(&correctflux, ef, theta, correctedFlux);
  FaceCavityOp op (&correctflux);
  op.applyToDimension(2);
  return correctedFlux;

}
}
