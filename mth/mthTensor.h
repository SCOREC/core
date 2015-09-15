/******************************************************************************

  Copyright 2015 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/

#ifndef MTH_TENSOR_H
#define MTH_TENSOR_H

#include "mthMatrix.h"

/** \file mthTensor.h
  * \brief Small run-time tensor. */

namespace mth {

/** \brief run-time (dynamic) tensor */
template <class T>
class Tensor : public Matrix<T>
{
  public:
    /** \brief default constructor */
    Tensor() : Matrix<T>() {}
    /** \brief construct with dimension d */
    Tensor(unsigned d) : mth::Matrix<T>(d,d)
    {
      assert((d==2) || (d==3));
    }
    /** \brief copy constructor */
    Tensor(Tensor<T> const& b) : mth::Matrix<T>() {copy(b);}
    /** \brief copy from a matrix */
    Tensor(Matrix<T> const& b) : mth::Matrix<T>() {copy(b);}
    /** \brief assignment operator */
    Tensor<T> operator=(Tensor<T> const& b) {copy(b);}
    /** \brief assignent to a matrix */
    Tensor<T> operator=(Matrix<T> const& b) {copy(b);}
    /** \brief add a tensor to a tensor */
    Tensor<T> operator+(Tensor<T> const& b)
    {
      Tensor<T> r(dim());
      for (unsigned i=0; i < dim(); ++i)
      for (unsigned j=0; j < dim(); ++j)
        r(i,j) = (*this)(i,j) + b(i,j);
      return r;
    }
    /** \brief subtract a tensor from a tensor */
    Tensor<T> operator-(Tensor<T> const& b)
    {
      Tensor<T> r(dim());
      for (unsigned i=0; i < dim(); ++i)
      for (unsigned j=0; j < dim(); ++j)
        r(i,j) = (*this)(i,j) - b(i,j);
      return r;
    }
    /** \brief multiply a tensor by a tensor */
    Tensor<T> operator*(Tensor<T> const& b)
    {
      Tensor<T> r(dim());
      for (unsigned i=0; i < dim(); ++i)
      for (unsigned j=0; j < dim(); ++j)
      {
        r(i,j) = (*this)(i,0) * b(0,j);
        for (unsigned k=1; k < dim(); ++k)
          r(i,j) += (*this)(i,k) * b(k,j);
      }
      return r;
    }
    /** \brief multiply a tensor by a scalar */
    Tensor<T> operator*(T const& s)
    {
      Tensor<T> r(dim());
      for (unsigned i=0; i < dim(); ++i)
      for (unsigned j=0; j < dim(); ++j)
        r(i,j) = (*this)(i,j) * s;
      return r;
    }
    /** \brief divide a tensor by a scalar */
    Tensor<T> operator/(T const& s)
    {
      Tensor<T> r(dim());
      for (unsigned i=0; i < dim(); ++i)
      for (unsigned j=0; j < dim(); ++j)
        r(i,j) = (*this)(i,j) / s;
      return r;
    }
    /** \brief resize this tensor to dimension d
      * \details d must be 2 or 3 */
    void resize(unsigned d)
    {
      assert((d==2) || (d==3));
      this->columns = d;
      this->elems.resize(d*d);
    }
    /** \brief resize this tensor
      * \details this method overrides Matrix.resize(m,n) to
      * ensure that a Tensor is not improperly resized */
    void resize(unsigned m, unsigned n)
    {
      assert(m==n);
      assert((m==2) || (m==3));
      this->columns = m;
      this->elems.resize(m*m);
    }
    /** \brief get the dimension of this tensor */
    unsigned dim() const {return this->columns;}
  private:
    void copy(Tensor<T> const& b)
    {
      unsigned d = b.dim();
      resize(d);
      for (unsigned i=0; i < d; ++i)
      for (unsigned j=0; j < d; ++j)
        (*this)(i,j) = b(i,j);
    }
    void copy(Matrix<T> const& b)
    {
      unsigned m = b.rows();
      unsigned n = b.cols();
      resize(m,n);
      for (unsigned i=0; i < m; ++i)
      for (unsigned j=0; j < m; ++j)
        (*this)(i,j) = b(i,j);
    }
};

}

#endif
