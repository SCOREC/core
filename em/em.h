/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef EM_H
#define EM_H


/** \file em.h
 *  \brief The Elegtromagnetics Equilibrated Residual error estimator inteface
 */
#include "apf.h"
#include <apfMesh.h>

#include <mthQR.h>
#include <mth.h>
#include <mth_def.h>

namespace em {

// Takes the solution electric field and corrected flux field and solves
// local element level BVPs to estimate the error.
// Returns a per-element scalar error field.
apf::Field* emEstimateError(apf::Field* ef, apf::Field* correctedFlux);

// Takes the solution electric field and equilibrated field (of face vectors)
// and computes the 'correction' to the flux vectors on each face.
apf::Field* computeFluxCorrection(apf::Field* ef, apf::Field* g);

// Takes the solution electric field and computes edge equilibrations.
apf::Field* equilibrateResiduals(apf::Field* f);


void assembleCurlCurlElementMatrix(apf::Mesh* mesh, apf::MeshEntity* e,
    apf::Field* f, mth::Matrix<double>& elmat);

void assembleVectorMassElementMatrix(apf::Mesh* mesh, apf::MeshEntity* e,
    apf::Field* f, mth::Matrix<double>& elmat);

void assembleElementMatrix(apf::Mesh* mesh, apf::MeshEntity*e,
    apf::Field* f, mth::Matrix<double>& elmat);

void assembleDomainLFElementVector(apf::Mesh* mesh, apf::MeshEntity* e,
    apf::Field* f, mth::Vector<double>& elvect);

// TODO QUESTION redo this to allow user access from outside
void pumiUserFunction(const apf::Vector3& x, mth::Vector<double>& f,
    apf::MeshEntity* e, apf::Mesh* mesh)
{
  double freq = 1.;
  double kappa = freq * M_PI;
  int dim = apf::getDimension(mesh, e);
  if (dim == 3) {
      f(0) = (1. + kappa * kappa) * sin(kappa * x[1]);
      f(1) = (1. + kappa * kappa) * sin(kappa * x[2]);
      f(2) = (1. + kappa * kappa) * sin(kappa * x[0]);
  }
  else {
     f(0) = (1. + kappa * kappa) * sin(kappa * x[1]);
     f(1) = (1. + kappa * kappa) * sin(kappa * x[0]);
     f(2) = 0.0;
  }
}

apf::Vector3 computeFaceOutwardNormal(apf::Mesh* m,
    apf::MeshEntity* t, apf::MeshEntity* f, apf::Vector3 const& p);


}











#endif
