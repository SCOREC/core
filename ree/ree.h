/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef REE_H
#define REE_H


/** \file ree.h
 *  \brief The Residual Error Estimator (REE) interface
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

/** \namespace ree
  * \brief All Residual Error Estimator functions
  */
namespace ree {

/** @brief Computes a nodal size field from the element error field
  *        obtained after running the residual error estimator.
  * @param ef (In) nedelec electric field
  * @param error_field (In) per-element residual error field
  * @param n a parameter to prescribe allowable error
  * @param alpha floor on the size field; alpha < h_new/h_old < beta
  * @param beta ceiling on the size field; alpha < h_new/h_old < beta
  * @returns a scalar mesh size field at mesh vertices
  */
apf::Field* getTargetEMSizeField(
    apf::Field* ef,
    apf::Field* error_field,
    int n,
    double alpha = 0.25,
    double beta = 2.0);

/** @brief Computes equilibrated residuals using the fem nedelec electric field.
  * @param f (In) nedelec electric field
  */
apf::Field* equilibrateResiduals(apf::Field* f);

/** @brief Uses the fem nedelec electric field and the equilibrated residuals to
  *        compute the 'correction' to the flux vectors on each face.
  * @param ef (In) nedelec electric field
  * @param g (In) equilibrated residuals field
  */
apf::Field* computeFluxCorrection(apf::Field* ef, apf::Field* g);

/** @brief Uses the fem nedelec electric field and the 'correction' to the flux
  *        vectors to compute the 'corrected' flux vectors on each face.
  * @param ef (In) nedelec electric field
  * @param theta (In) correction to the flux vectors
  */
apf::Field* computeCorrectedFlux(apf::Field* ef, apf::Field* theta);

/** @brief Solves additional local BVPs on a one order higher nedelec field on
  *        each element to estimate the dsicretization error.
  * @param ef (In) nedelec electric field
  * @param correctedFlux (In) flux field which provides Neumann BCs for each
  *        local BVP.
  * @returns a per-element scalar residual error field.
  */
apf::Field* computeErrorField(apf::Field* ef, apf::Field* correctedFlux);

/** @brief run the residual error estimator.
  * @param f the fem nedelec electric field
  *        scales the output size field.
  * @returns a per-element scalar residual error field.
  */
apf::Field* estimateError(apf::Field* f);

// helper functions
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

bool isOnDomainBoundary(apf::Mesh* m, apf::MeshEntity* e);

}

#endif
