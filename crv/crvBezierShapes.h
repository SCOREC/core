/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef CRVBEZIERSHAPES_H
#define CRVBEZIERSHAPES_H

namespace crv {

/** \brief polynomial part of bernstein polynomial */
double Bij(const int i, const int j,const double u, const double v);
double Bijk(const int i, const int j, const int k, const double u,
    const double v, const double w);
double Bijkl(const int i, const int j, const int k, const int l,
    const double u, const double v, const double w, const double t);

/** \brief a different form of the above */
double Bij(const int ij[], const double xi[]);
double Bijk(const int ijk[], const double xi[]);
double Bijkl(const int ijkl[], const double xi[]);

typedef void (*bezierShape)(int P,
    apf::Vector3 const& xi,
    apf::NewArray<double>& values);

typedef void (*bezierShapeGrads)(int P,
    apf::Vector3 const& xi,
    apf::NewArray<apf::Vector3>& grads);

extern const bezierShape bezier[apf::Mesh::TYPES];

extern const bezierShapeGrads bezierGrads[apf::Mesh::TYPES];

}

#endif
