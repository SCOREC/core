/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "spr.h"
#include <PCU.h>
#include <cmath>

namespace spr {

void backSub(apf::DynamicMatrix& R,
             apf::DynamicVector& b,
             apf::DynamicVector& x)
{
  int n = R.getColumns();
  x.setSize(n);
  for (int i=n-1; i>=0; i--) {
    x(i) = b(i);
    for (int j=i+1; j<n; j++) {
      x(i) -= R(i,j)*x(j);
    }
    x(i) = x(i)/R(i,i);
  }
}

double sign(double x) {
  double r = (x >= 0.0) ? 1.0 : -1.0;
  return r;
}

void solveQR(apf::DynamicMatrix& A,
             apf::DynamicVector& x,
             apf::DynamicVector& b)
{
  int m = A.getRows();
  int n = A.getColumns();
  apf::DynamicVector vk;
  apf::DynamicMatrix T;
  double Tb;
  for (int k=0; k<n; k++) {
    vk.setSize(m-k);
    for (int j=0; j<(m-k); j++) {
      vk(j) = A(k+j,k);
    }
    double v = sign(vk(0))*vk.getLength();
    vk(0) += v;
    vk /= vk.getLength();
    T.setSize(m-k,n-k);
    Tb = 0;
    for (int i=0; i<m-k; i++) {
      Tb += vk(i)*b(i+k);
      for (int j=0; j<n-k; j++) {
        T(i,j) = 0.0;
        for (int l=0; l<m-k; l++) {
          T(i,j) += 2.0 * vk(i) * vk(l) * A(l+k,j+k);
        }
      }
    }
    for (int i=0; i<m-k; i++) {
      b(i+k) -= 2*vk(i)*Tb;
      for (int j=0; j<n-k; j++) {
        A(i+k,j+k) -= T(i,j);
      }
    }
  }
  backSub(A,b,x);
}
}
