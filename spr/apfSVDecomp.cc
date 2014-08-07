/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfSVDecomp.h"

#include <stdio.h>
#include <cassert>
#include <math.h>

namespace apf {


/* returns absolute value of a times the sign of b */
double svdSign(double a, double b) 
{
  double result = ( b >= 0.0 ) ? fabs(a) : -fabs(a);
  return result;
}


/* returns c for a^2 + b^2 = c^2 */
double svdPythag(double a, double b)
{
  double at = fabs(a), bt = fabs(b), ct, result;
  if (at > bt)       { ct = bt / at; result = at * sqrt(1.0 + ct * ct); }
  else if (bt > 0.0) { ct = at / bt; result = bt * sqrt(1.0 + ct * ct); }
  else result = 0.0;
  return(result);
}


/* returns max of a and b */
double svdMax(double a, double b)
{
  double result = ( a > b ) ? a : b;
  return result;
}

void svdSolve(DynamicMatrix& A,
              DynamicVector& b,
              int maxIters,
              double tol,
              DynamicVector& x)
{
  DynamicVector w;
  DynamicMatrix V;
  svdDecompose(A,w,V,maxIters);

  /* invert W and zero out singular values */
  for (size_t i=0; i<w.getSize(); i++)
  {
    if ( w(i) < tol )
      w(i) = 0.0;
    else 
      w(i) = 1.0/w[i];
  }
  
  /* multiply V * inverse(w) (store value in V) */
  for (size_t i=0; i<V.getRows(); i++) 
  {
    for (size_t j=0; j<V.getColumns(); j++)
      V(i,j) *= w(j);
  }

  DynamicMatrix invU;
  transpose(A, invU);
  DynamicMatrix tmp;
  multiply(V, invU, tmp);
  multiply(tmp,b,x);
}

void svdDecompose(DynamicMatrix& A,
                  DynamicVector& w,
                  DynamicMatrix& V,
                  int maxIters) 
{
  int m, n;
  m = A.getRows();
  n = A.getColumns();

  if (m < n) 
    std::cerr << "num rows must be > num cols \n";

  w.setSize(n);
  V.setSize(n,n);

  int flag, i, its, j, jj, k, l, nm;
  double c, f, h, s, x, y, z;
  double anorm = 0.0, g = 0.0, scale = 0.0;
  DynamicVector rv1(n);

  l = 0; /* silence warning on less-able compilers */

  /* Householder reduction to bidiagonal form */
  for (i = 0; i < n; i++)
  {
    /* left-hand reduction */
    l = i + 1;
    rv1(i) = scale * g;
    g = s = scale = 0.0;
    if (i < m)
    {
      for (k = i; k < m; k++) 
        scale += fabs(A(k,i));
      if (scale) 
      {
        for (k = i; k < m; k++)
        {
          A(k,i) = (A(k,i)/scale);
          s += (A(k,i) * A(k,i));
        }
        f = A(i,i);
        g = -svdSign(sqrt(s), f);
        h = f * g - s;
        A(i,i) = (f - g);
        if (i != n - 1)
        {
          for (j = l; j < n; j++)
          {
            for (s = 0.0, k = i; k < m; k++) 
              s += (A(k,i) * A(k,j));
            f = s / h;
            for (k = i; k < m; k++) 
              A(k,j) += (f * A(k,i));
          }
        }
        for (k = i; k < m; k++) 
          A(k,i) = (A(k,i)*scale);
      }
    }
    w(i) = (scale * g);

    /* right-hand reduction */
    g = s = scale = 0.0;
    if (i < m && i != n - 1)
    {
      for (k = l; k < n; k++) 
        scale += fabs(A(i,k));
      if (scale)
      {
        for (k = l; k < n; k++)
        {
          A(i,k) = (A(i,k)/scale);
          s += (A(i,k) * A(i,k));
        }
        f = A(i,l);
        g = -svdSign(sqrt(s), f);
        h = f * g - s;
        A(i,l) = (f - g);
        for (k = l; k < n; k++) 
          rv1(k) = A(i,k) / h;
        if (i != m - 1) 
        {
          for (j = l; j < m; j++)
          {
            for (s = 0.0, k = l; k < n; k++) 
              s += (A(j,k) * A(i,k));
            for (k = l; k < n; k++) 
              A(j,k) += (s * rv1(k));
          }
        }
        for (k = l; k < n; k++) 
          A(i,k) = (A(i,k)*scale);
      }
    }
    anorm = svdMax(anorm, (fabs(w(i)) + fabs(rv1(i))));
  }

  /* accumulate the right-hand transformation */
  for (i = n - 1; i >= 0; i--)
  {
    if (i < n - 1)
    {
      if (g) 
      {
        for (j = l; j < n; j++)
          V(j,i) = ((A(i,j) / A(i,l)) / g);
        for (j = l; j < n; j++)
        {
          for (s = 0.0, k = l; k < n; k++) 
            s += (A(i,k) * V(k,j));
          for (k = l; k < n; k++) 
            V(k,j) += (s * V(k,i));
        }
      }
      for (j = l; j < n; j++) 
        V(i,j) = V(j,i) = 0.0;
    }
    V(i,i) = 1.0;
    g = rv1(i);
    l = i;
  }

  /* accumulate the left-hand transformation */
  for (i = n - 1; i >= 0; i--)
  {
    l = i + 1;
    g = w(i);
    if (i < n - 1) 
      for (j = l; j < n; j++) 
        A(i,j) = 0.0;
    if (g) 
    {
      g = 1.0 / g;
      if (i != n - 1) {
        for (j = l; j < n; j++)
        {
          for (s = 0.0, k = l; k < m; k++) 
            s += (A(k,i) * A(k,j));
          f = (s / A(i,i)) * g;
          for (k = i; k < m; k++) 
            A(k,j) += (f * A(k,i));
        }
      }
      for (j = i; j < m; j++) 
        A(j,i) = (A(j,i)*g);
    }
    else {
      for (j = i; j < m; j++) 
        A(j,i) = 0.0;
    }
    ++A(i,i);
  }

  /* diagonalize the bidiagonal form */
  for (k = n - 1; k >= 0; k--) 
  {
    for (its = 0; its < maxIters; its++) 
    {
      flag = 1;
      for (l = k; l >= 0; l--) 
      {
        nm = l - 1;
        if (fabs(rv1(l)) + anorm == anorm) 
        {
          flag = 0;
          break;
        }
        if (fabs(w(nm)) + anorm == anorm) 
          break;
      }
      if (flag) 
      {
        c = 0.0;
        s = 1.0;
        for (i = l; i <= k; i++) 
        {
          f = s * rv1(i);
          if (fabs(f) + anorm != anorm) 
          {
            g = w(i);
            h = svdPythag(f, g);
            w(i) = h; 
            h = 1.0 / h;
            c = g * h;
            s = (- f * h);
            for (j = 0; j < m; j++) 
            {
              y = A(j,nm);
              z = A(j,i);
              A(j,nm) = (y * c + z * s);
              A(j,i) = (z * c - y * s);
            }
          }
        }
      }
      z = w(k);

      /* convergence */
      if (l == k) 
      {
        if (z < 0.0)
        {
          w(k) = (-z);
          for (j = 0; j < n; j++) 
            V(j,k) = (-V(j,k));
        }
        break;
      }
      if (its >= maxIters) 
        std::cerr << "No convergence after " << maxIters << " iterations \n";

      /* shift from bottom 2x2 minor */
      x = w(l);
      nm = k - 1;
      y = w(nm);
      g = rv1(nm);
      h = rv1(k);
      f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
      g = svdPythag(f, 1.0);
      f = ((x - z) * (x + z) + h * ((y / (f + svdSign(g, f))) - h)) / x;

      /* next QR transformation */
      c = s = 1.0;
      for (j = l; j <= nm; j++) 
      {
        i = j + 1;
        g = rv1(i);
        y = w(i);
        h = s * g;
        g = c * g;
        z = svdPythag(f, h);
        rv1(j) = z;
        c = f / z;
        s = h / z;
        f = x * c + g * s;
        g = g * c - x * s;
        h = y * s;
        y = y * c;
        for (jj = 0; jj < n; jj++)
        {
          x = V(jj,j);
          z = V(jj,i);
          V(jj,j) = (x * c + z * s);
          V(jj,i) = (z * c - x * s);
        }
        z = svdPythag(f, h);
        w(j) = z;
        if (z)
        {
          z = 1.0 / z;
          c = f * z;
          s = h * z;
        }
        f = (c * g) + (s * y);
        x = (c * y) - (s * g);
        for (jj = 0; jj < m; jj++)
        {
          y = A(jj,j);
          z = A(jj,i);
          A(jj,j) = (y * c + z * s);
          A(jj,i) = (z * c - y * s);
        }
      }
      rv1(l) = 0.0;
      rv1(k) = f;
      w(k) = x;
    }
  }
}

}
