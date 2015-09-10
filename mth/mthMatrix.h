#ifndef MTH_MATRIX_H
#define MTH_MATRIX_H

#include "mthVector.h"

/** \file mthMatrix.h
  * \brief Small compile time linear algebra matrices. */

namespace mth {

/** \brief compile time (static) matrix */
template <class T, unsigned M=0, unsigned N=0>
class Matrix : public Array<Vector<T,N>,M>
{
  public:
    /** \brief mutable index operator
      * \details an index operator (i,j) is provided so that
      * compile time matrices and runtime-sized matrices can be used
      * interchangebly, without worrying about changing index syntax */
    T& operator()(unsigned i, unsigned j) {return (*this)[i][j];}
    /** \brief immutable index operator
      * \details see the mutable index operator details */
    T const& opertor()(unsigned i, unsigned j) const {return (*this)[i][j];}
    /** \brief default constructor */
    Matrix() {}
    /** \brief add two matrices */
    Matrix<T,M,N> operator+(Matrix<T,M,N> const& b) const
    {
      Matrix<T,M,N> r;
      for (unsigned i=0; i < M; ++i)
        r.elems[i] = this->elems[i] + b.elems[i];
      return r;
    }
    /** \brief subtract two matrices */
    Matrix<T,M,N> operator-(Matrix<T,M,N> const& b) const
    {
      Matrix<T,M,N> r;
      for (unsigned i=0; i < M; ++i)
        r.elems[i] = this->elems[i] - b.elems[i];
      return r;
    }
    /** \brief multiply by a scalar */
    Matrix<T,M,N> operator*(T const& s) const
    {
      Matrix<T,M,N> r;
      for (unsigned i=0; i < M; ++i)
        r.elems[i] = this->elems[i] * s;
    }
    /** \brief divide by a scalar */
    Matrix<T,M,N> operator/(T const& s) const
    {
      Matrix<T,M,N> r;
      for (unsigned i=0; i < M; ++i)
        r.elems[i] = this->elems[i] / s;
      return r;
    }
    /** \brief multiply a matrix by a vector */
    Vector<T,M> operator*(Vector<T,N> const& b) const
    {
      Vector<T,M> r;
      for (unsigned i=0; i < M; ++i)
        r[i] = this->elems[i] * b;
      return r;
    }
    /** \brief multiply two matrices
      * \details the extra template parameter generates
      * code for all possible combinations of matrix sizes */
    template <unsigned O>
    Matrix<T,M,O> operator*(Matrix<T,N,O> const& b) const
    {
      Matrix<T,M,O> r;
      for (unsigned i=0; i < M; ++i)
      for (unsigned j=0; j < O; ++j)
      {
        r[i][j] = this->elems[i][0]*b[0][j];
        for (unsigned k=1; k < N; ++k)
          r[i][j] += this->elems[i][k]*b[k][j];
      }
      return r;
    }
};

}

#endif
