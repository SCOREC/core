/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFMATRIX_H
#define APFMATRIX_H

#include "apfVector.h"

namespace apf {

template <std::size_t M, std::size_t N>
//matrix is an array of row vectors
class Matrix : public Array<Vector<N>,M>
{
  public:
    Matrix<M,N> operator+(Matrix<M,N> const& b) const
    {
      Matrix<M,N> r;
      for (std::size_t i=0; i < M; ++i)
        r.elements[i] = this->elements[i] + b.elements[i];
      return r;
    }
    Matrix<M,N> operator-(Matrix<M,N> const& b) const
    {
      Matrix<M,N> r;
      for (std::size_t i=0; i < M; ++i)
        r.elements[i] = this->elements[i] - b.elements[i];
      return r;
    }
    Matrix<M,N> operator*(double s) const
    {
      Matrix<M,N> r;
      for (std::size_t i=0; i < M; ++i)
        r.elements[i] = this->elements[i] * s;
      return r;
    }
    Matrix<M,N> operator/(double s) const
    {
      Matrix<M,N> r;
      for (std::size_t i=0; i < M; ++i)
        r.elements[i] = this->elements[i] / s;
      return r;
    }
    Vector<M> operator*(Vector<N> const& b) const
    {
      Vector<M> r;
      for (std::size_t i=0; i < M; ++i)
        r[i] = this->elements[i] * b;
      return r;
    }
    template <std::size_t O>
    Matrix<M,O> operator*(Matrix<N,O> const& b) const
    {
      Matrix<M,O> r;
      for (std::size_t i=0; i < M; ++i)
      for (std::size_t j=0; j < O; ++j)
      {
        r[i][j] = this->elements[i][0]*b[0][j];
        for (std::size_t k=1; k < N; ++k)
          r[i][j] += this->elements[i][k]*b[k][j];
      }
      return r;
    }
};

template <std::size_t M, std::size_t N>
Matrix<N,M> transpose(Matrix<M,N> const& m)
{
  Matrix<N,M> r;
  for (std::size_t i=0; i < M; ++i)
  for (std::size_t j=0; j < N; ++j)
    r[j][i] = m[i][j];
  return r;
}

template <std::size_t M, std::size_t N>
Matrix<M,N> tensorProduct(Vector<M> const& a, Vector<N> const& b)
{
  Matrix<M,N> r;
  for (std::size_t i=0; i < M; ++i)
    r[i] = b * a[i];
  return r;
}

inline double det(Matrix<2,2> const& m)
{
  return m[0][0] * m[1][1] - m[1][0] * m[0][1];
}

inline double det(Matrix<3,3> const& m)
{
  return m[0]*cross(m[1],m[2]);
}

inline Matrix<3,3> cofactor(Matrix<3,3> const &m)
{
  Matrix<3,3> r;

  r[0][0] = m[1][1] * m[2][2] - m[2][1] * m[1][2];
  r[0][1] = m[2][1] * m[0][2] - m[2][2] * m[0][1];
  r[0][2] = m[0][1] * m[1][2] - m[1][1] * m[0][2];

  r[1][0] = m[1][2] * m[2][0] - m[1][0] * m[2][2];
  r[1][1] = m[0][0] * m[2][2] - m[0][2] * m[2][0];
  r[1][2] = m[1][0] * m[0][2] - m[0][0] * m[1][2];
  
  r[2][0] = m[1][0] * m[2][1] - m[2][0] * m[1][1];
  r[2][1] = m[2][0] * m[0][1] - m[0][0] * m[2][1];
  r[2][2] = m[0][0] * m[1][1] - m[0][1] * m[1][0];

  return r;
}

inline Matrix<3,3> invert(Matrix<3,3> const& m)
{
  Matrix<3,3> x = transpose(m);
  Matrix<3,3> r;
  r[0] = cross(x[1],x[2]);
  r[1] = cross(x[2],x[0]);
  r[2] = cross(x[0],x[1]);
  return r/det(m);
}

class Matrix3x3 : public Matrix<3,3>
{
  public:
    Matrix3x3() {}
    Matrix3x3(double a11, double a12, double a13,
              double a21, double a22, double a23,
              double a31, double a32, double a33)
    {
      this->elements[0] = Vector3(a11,a12,a13);
      this->elements[1] = Vector3(a21,a22,a23);
      this->elements[2] = Vector3(a31,a32,a33);
    }
    Matrix3x3(Matrix<3,3> const& other):
      Matrix<3,3>(other)
    {}
    Matrix3x3(double (*array)[3])
    {
      for (std::size_t i=0; i < 3; ++i)
      for (std::size_t j=0; j < 3; ++j)
        this->elements[i][j] = array[i][j];
    }
    void toArray(double (*array)[3]) const
    {
      for (std::size_t i=0; i < 3; ++i)
      for (std::size_t j=0; j < 3; ++j)
        array[i][j] = this->elements[i][j];
    }
};

template <std::size_t M, std::size_t N>
double getInnerProduct(Matrix<M,N> const& a,
                       Matrix<M,N> const& b)
{
  double r = a[0]*b[0];
  for (std::size_t i=1; i < M; ++i)
    r += a[i]*b[i];
  return r;
}

Matrix3x3 cross(Vector3 const& u);
Matrix3x3 rotate(Vector3 const& u, double a);
Matrix3x3 getFrame(Vector3 const& v);

int eigen(Matrix3x3 const& A,
          Vector3* eigenVectors,
          double* eigenValues);

}//namespace apf

std::ostream& operator<<(std::ostream& s, apf::Matrix3x3 const& v);

#endif
