/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <PCU.h>
#include "spr.h"
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

/* V gets the Householder vectors as
   its rows upon return.
   this enables getImplicitB later on */
bool decompQR(apf::DynamicMatrix& A,
              apf::DynamicMatrix& V,
              apf::DynamicMatrix& R)
{
  int m = A.getRows();
  int n = A.getColumns();
  R = A;
  apf::DynamicVector vk;
  apf::DynamicMatrix T;
  V.setSize(n, m);
  T.setSize(m, n);
  vk.setSize(m);
  for (int k=0; k<n; k++) {
    for (int j=0; j < (m-k); j++) {
      vk(j) = R(k+j,k);
    }
    for (int j=(m-k); j < m; j++) {
      vk(j) = 0;
    }
    double lvk = vk.getLength();
    if (fabs(lvk) < 1e-10) /* tolerance check */
      return false;
    double v = sign(vk(0))*lvk;
    vk(0) += v;
    vk /= vk.getLength();
    for (int i=0; i<m-k; i++) {
      for (int j=0; j<n-k; j++) {
        T(i,j) = 0.0;
        for (int l=0; l<m-k; l++) {
          T(i,j) += 2.0 * vk(i) * vk(l) * R(l+k,j+k);
        }
      }
    }
    for (int i=0; i<m-k; i++) {
      for (int j=0; j<n-k; j++) {
        R(i+k,j+k) -= T(i,j);
      }
    }
    for (int j=0; j < m; j++) {
      V(k,j) = vk(j);
    }
  }
  return true;
}

/* computes Q'*b using the Householder
   vectors V(k,:) */
void getImplicitB(apf::DynamicMatrix& V,
                  apf::DynamicVector& b,
                  apf::DynamicVector& Qtb)
{
  int n = V.getRows();
  int m = V.getColumns();
  Qtb = b;
  for (int k=0; k<n; k++) {
    double Tb = 0;
    for (int i=0; i<m-k; i++) {
      Tb += V(k,i) * Qtb(i+k);
    }
    for (int i=0; i<m-k; i++) {
      Qtb(i+k) -= 2 * V(k,i) * Tb;
    }
  }
}

void solveFromQR(apf::DynamicMatrix& V,
                 apf::DynamicMatrix& R,
                 apf::DynamicVector& b,
                 apf::DynamicVector& x)
{
  apf::DynamicVector Qtb;
  getImplicitB(V, b, Qtb);
  backSub(R, Qtb, x);
}

}
