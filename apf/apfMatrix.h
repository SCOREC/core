/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFMATRIX_H
#define APFMATRIX_H

/** \file apfMatrix.h
  \brief The APF linear algebra matrix interface */

#include "apfVector.h"

namespace apf {

/** \brief template-generic matrix of M by N doubles
  \details see apf::Vector for the rationale on templating.
  In short, this class is designed to handle matrices
  whose small sizes are known at compile time.
  They should be used in a functional programming style

  For matrices sized at runtime, see apf::DynamicMatrix.
  For sparse structures or parallel matrices, look
  outside of APF.
 
  for those interested in software design, notice how
  Array and Vector come together to form Matrix */
template <std::size_t M, std::size_t N>
class Matrix : public Array<Vector<N>,M>
{
  public:
    /** \brief mandatory */
    Matrix() {}
    /** \brief construct from an array */
    Matrix(double const (*array)[N])
    {
      for (std::size_t i=0; i < M; ++i)
      for (std::size_t j=0; j < N; ++j)
        this->elements[i][j] = array[i][j];
    }
    /** \brief add two matrices */
    Matrix<M,N> operator+(Matrix<M,N> const& b) const
    {
      Matrix<M,N> r;
      for (std::size_t i=0; i < M; ++i)
        r.elements[i] = this->elements[i] + b.elements[i];
      return r;
    }
    /** \brief subtract two matrices */
    Matrix<M,N> operator-(Matrix<M,N> const& b) const
    {
      Matrix<M,N> r;
      for (std::size_t i=0; i < M; ++i)
        r.elements[i] = this->elements[i] - b.elements[i];
      return r;
    }
    /** \brief multiply a matrix by a scalar */
    Matrix<M,N> operator*(double s) const
    {
      Matrix<M,N> r;
      for (std::size_t i=0; i < M; ++i)
        r.elements[i] = this->elements[i] * s;
      return r;
    }
    /** \brief divide a matrix by a scalar */
    Matrix<M,N> operator/(double s) const
    {
      Matrix<M,N> r;
      for (std::size_t i=0; i < M; ++i)
        r.elements[i] = this->elements[i] / s;
      return r;
    }
    /** \brief multiply a matrix by a vector */
    Vector<M> operator*(Vector<N> const& b) const
    {
      Vector<M> r;
      for (std::size_t i=0; i < M; ++i)
        r[i] = this->elements[i] * b;
      return r;
    }
    /** \brief multiply two matrices
      \details the extra template parameter generates
      code for all possible combinations of matrix
      sizes */
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

/** \brief transpose a matrix */
template <std::size_t M, std::size_t N>
Matrix<N,M> transpose(Matrix<M,N> const& m)
{
  Matrix<N,M> r;
  for (std::size_t i=0; i < M; ++i)
  for (std::size_t j=0; j < N; ++j)
    r[j][i] = m[i][j];
  return r;
}

/** \brief tensor product of two vectors */
template <std::size_t M, std::size_t N>
Matrix<M,N> tensorProduct(Vector<M> const& a, Vector<N> const& b)
{
  Matrix<M,N> r;
  for (std::size_t i=0; i < M; ++i)
    r[i] = b * a[i];
  return r;
}

/** \brief get the minor matrix associated with entry (i,j) of matrix A
 \details this is only instantiated for square matrices up to 4 by 4 */
template <std::size_t M, std::size_t N>
Matrix<M-1,N-1> getMinor(Matrix<M,N> const& A, std::size_t i, std::size_t j);

/** \brief get the cofactor associated with entry (i,j) of matrix A
 \details this is only instantiated for square matrices up to 4 by 4 */
template <std::size_t M, std::size_t N>
double getCofactor(Matrix<M,N> const& A, std::size_t i, std::size_t j);

/** \brief get the determinant of a matrix A
 \details this is only instantiated for square matrices up to 4 by 4 */
template <std::size_t M, std::size_t N>
double getDeterminant(Matrix<M,N> const& A);

/** \brief get the matrix of cofactors for a given matrix */
inline Matrix<3,3> cofactor(Matrix<3,3> const &m)
{
  Matrix<3,3> r;
  for (std::size_t i = 0; i < 3; ++i)
  for (std::size_t j = 0; j < 3; ++j)
    r[i][j] = getCofactor(m, i, j);
  return r;
}

/** \brief invert a 2 by 2 matrix */
inline Matrix<2,2> invert(Matrix<2,2> const& m)
{
  Matrix<2,2> a;
  a[0][0] =  m[1][1]; a[0][1] = -m[0][1];
  a[1][0] = -m[1][0]; a[1][1] =  m[0][0];
  return a / getDeterminant(m);
}

/** \brief invert a 3 by 3 matrix */
inline Matrix<3,3> invert(Matrix<3,3> const& m)
{
  Matrix<3,3> x = transpose(m);
  Matrix<3,3> r;
  r[0] = cross(x[1],x[2]);
  r[1] = cross(x[2],x[0]);
  r[2] = cross(x[0],x[1]);
  return r / getDeterminant(m);
}

/** \brief convenience wrapper over apf::Matrix<3,3>
  \details like apf::Vector3, this provides component-wise
  initialization */
class Matrix3x3 : public Matrix<3,3>
{
  public:
    /** \brief required default constructor */
    Matrix3x3() {}
    /** \brief component-wise constructor
     \details this is useful for hardcoded matrices */
    Matrix3x3(double a11, double a12, double a13,
              double a21, double a22, double a23,
              double a31, double a32, double a33)
    {
      this->elements[0] = Vector3(a11,a12,a13);
      this->elements[1] = Vector3(a21,a22,a23);
      this->elements[2] = Vector3(a31,a32,a33);
    }
    /** \brief constructor from base type */
    Matrix3x3(Matrix<3,3> const& other):
      Matrix<3,3>(other)
    {}
    /** \brief write matrix to an array
      \todo this could be templated */
    void toArray(double (*array)[3]) const
    {
      for (std::size_t i=0; i < 3; ++i)
      for (std::size_t j=0; j < 3; ++j)
        array[i][j] = this->elements[i][j];
    }
};

/** \brief get the component-wise inner product of two matrices */
template <std::size_t M, std::size_t N>
double getInnerProduct(Matrix<M,N> const& a,
                       Matrix<M,N> const& b)
{
  double r = a[0]*b[0];
  for (std::size_t i=1; i < M; ++i)
    r += a[i]*b[i];
  return r;
}

/** \brief get the skew-symmetric cross product matrix of a vector */
Matrix3x3 cross(Vector3 const& u);
/** \brief get the rotation matrix around an axis
  \param u the axis to rotate around
  \param a the amount of rotation in radians */
Matrix3x3 rotate(Vector3 const& u, double a);
/** \brief derive an orthogonal frame whose x axis is the given vector
  \details this is a robust algorithm for choosing an arbitrary
  frame that is aligned with a given vector while avoiding the
  numerical issues that arise when the vector is near global axes.

  Note, however, that the resulting frame is not normalized.
  some uses of this function require that the x basis vector be
  exactly the input, so that is the default behavior.
  It is the user's responsibility to normalize the result if desired */
Matrix3x3 getFrame(Vector3 const& v);

/** \brief get the eigenvectors and eigenvalues of a 3 by 3 matrix */
int eigen(Matrix3x3 const& A,
          Vector<3>* eigenVectors,
          double* eigenValues);

}//namespace apf

/** \brief write a 1 by 1 matrix to a C++ stream */
std::ostream& operator<<(std::ostream& s, apf::Matrix<1,1> const& A);
/** \brief write a 2 by 2 matrix to a C++ stream */
std::ostream& operator<<(std::ostream& s, apf::Matrix<2,2> const& A);
/** \brief write a 3 by 3 matrix to a C++ stream */
std::ostream& operator<<(std::ostream& s, apf::Matrix<3,3> const& A);
/** \brief write a 4 by 4 matrix to a C++ stream */
std::ostream& operator<<(std::ostream& s, apf::Matrix<4,4> const& A);

#endif
