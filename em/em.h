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

#include <PCU.h>
#include <lionPrint.h>

#include "apfShape.h"
#include "apfField.h"

#include <mthQR.h>
#include <mth.h>
#include <mth_def.h>

namespace em {

/*
 * Computes nodal size field // TODO currently only element size field
 */
apf::Field* getTargetEMSizeField(
    apf::Field* ef,
    apf::Field* error_field,
    double alpha,
    double beta);
/*
 * Takes the solution electric field and computes edge equilibrations.
 */
apf::Field* equilibrateResiduals(apf::Field* f);

/*
 * Takes the solution electric field and equilibrated field (of face vectors)
 * and computes the 'correction' to the flux vectors on each face.
 */
apf::Field* computeFluxCorrection(apf::Field* ef, apf::Field* g);

/* Takes the solution electric field and corrected flux field and solves
 * local element level BVPs to estimate the error.
 * Returns a per-element scalar error field.
 */
apf::Field* computeErrorField(apf::Field* ef, apf::Field* THETA_Field);

apf::Field* estimateError(apf::Field* f);


void assembleCurlCurlElementMatrix(apf::Mesh* mesh, apf::MeshEntity* e,
    apf::Field* f, mth::Matrix<double>& elmat);

void assembleVectorMassElementMatrix(apf::Mesh* mesh, apf::MeshEntity* e,
    apf::Field* f, mth::Matrix<double>& elmat);

void assembleElementMatrix(apf::Mesh* mesh, apf::MeshEntity*e,
    apf::Field* f, mth::Matrix<double>& elmat);

void assembleDomainLFElementVector(apf::Mesh* mesh, apf::MeshEntity* e,
    apf::Field* f, mth::Vector<double>& elvect);

apf::Vector3 computeFaceOutwardNormal(apf::Mesh* m,
    apf::MeshEntity* t, apf::MeshEntity* f, apf::Vector3 const& p);


}











#endif
