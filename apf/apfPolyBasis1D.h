/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, incensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directroy.
 */

#ifndef APFPOLYBASES1D_H
#define APFPOLYBASES1D_H

#include "apfNew.h"

namespace apf {

/** \brief Gauss point types */
enum Gauss{
  GAUSS_LEGENDRE,
  GAUSS_LOBATTO,
  GAUSS
};


/** \brief Gauss Legendre points */
void getGaussLegendrePoints(int np, double* pts);

/** \brief Gauss Lobatto points */
void getGaussLobattoPoints(int, double*);

/** \brief get open points for a given order and Gauss point type */
void getOpenPoints(
    int order,
    NewArray<double>& op,
    int type = GAUSS_LEGENDRE);

/** \brief get closed points for a given order and Gauss point type */
void getClosedPoints(
    int order,
    NewArray<double>& cp,
    int type = GAUSS_LOBATTO);

/** \brief Chebyshev polynomials of the first kind */
void getChebyshevT(int order, double xi, double* u);

/** \brief Chebyshev polynomials of the first kind and 1st derivative */
void getChebyshevT(int order, double xi, double* u, double* d);

/** \brief Chebyshev polynomials of the first kind and 1st and 2nd derivative */
void getChebyshevT(int order, double xi, double* u, double* d, double* dd);

}

#endif
