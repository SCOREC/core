/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFDYNAMICMATRIX_H
#define APFDYNAMICMATRIX_H

/** \file apfDynamicMatrix.h
  \brief Small runtime-sized matrices */

#include "apfDynamicVector.h"
#include "apfMatrix.h"
#include "iostream"

namespace apf {

/** \brief A runtime-sized dense matrix
  \details see apf::DynamicVector for some general
  guidance on apf::DynamicMatrix as opposed to apf::Matrix.
  This class is meant to be used for small matrices whose
  size is not known at compile time.
  For big, sparse, or parallel matrices, look outside of APF. */
class DynamicMatrix
{
  public:
    /** \brief defautl constructor, no allocation */
    DynamicMatrix() {}
    /** \brief construct with size m by n */
    DynamicMatrix(std::size_t m, std::size_t n):
      columns(n),
      values(m*n)
    {}
    /** \brief get the number of rows (first index) */
    std::size_t getRows() const {return values.getSize()/columns;}
    /** \brief get the number of columns (second index) */
    std::size_t getColumns() const {return columns;}
    /** \brief resize to m by n */
    void setSize(std::size_t m, std::size_t n)
    {
      columns = n;
      values.setSize(m*n);
    }
    /** \brief immutable index operator */
    double operator()(std::size_t i, std::size_t j) const
    {
      return values[i*columns + j];
    }
    /** \brief mutable index operator */
    double& operator()(std::size_t i, std::size_t j)
    {
      return values[i*columns + j];
    }
    /** \brief add a matrix to this matrix */
    DynamicMatrix& operator+=(DynamicMatrix const& b)
    {
      for (std::size_t i=0; i < this->values.getSize(); ++i)
        this->values[i] += b.values[i];
      return *this;
    }
    /** \brief subtract a matrix from this matrix */
    DynamicMatrix& operator-=(DynamicMatrix const& b)
    {
      for (std::size_t i=0; i < this->values.getSize(); ++i)
        this->values[i] += b.values[i];
      return *this;
    }
    /** \brief multiply this matrix by a scalar */
    DynamicMatrix& operator*=(double s)
    {
      for (std::size_t i=0; i < this->values.getSize(); ++i)
        this->values[i] *= s;
      return *this;
    }
    /** \brief divide this matrix by a scalar */
    DynamicMatrix& operator/=(double s)
    {
      for (std::size_t i=0; i < this->values.getSize(); ++i)
        this->values[i] *= s;
      return *this;
    }
    /** \brief copy row data into a DynamicVector */
    void getRow(std::size_t i, DynamicVector& r) const
    {
      r.setSize(columns);
      for (std::size_t j=0; j < columns; ++j)
        r(j) = (*this)(i,j);
    }
    /** \brief copy column data into a DynamicVector */
    void getColumn(std::size_t j, DynamicVector& r) const
    {
      std::size_t rows = getRows();
      r.setSize(rows);
      for (std::size_t i=0; i < rows; ++i)
        r(i) = (*this)(i,j);
    }
    /** \brief copy row data from a DynamicVector */
    void setRow(std::size_t i, DynamicVector const& r)
    {
      for (std::size_t j=0; j < columns; ++j)
        (*this)(i,j) = r(j);
    }
    /** \brief copy column data from a DynamicVector */
    void setColumn(std::size_t j, DynamicVector const& r)
    {
      std::size_t rows = getRows();
      for (std::size_t i=0; i < rows; ++i)
        (*this)(i,j) = r(i);
    }
  protected:
    std::size_t columns;
    DynamicArray<double> values;
};

/** \brief multiply a DynamicMatrix by a DynamicVector */
inline void multiply(DynamicMatrix const& a,
                     DynamicVector const& b,
                     DynamicVector& r)
{
  std::size_t rows = a.getRows();
  std::size_t columns = a.getColumns();
  r.setSize(rows);
  for (std::size_t i=0; i < rows; ++i)
  {
    r[i] = a(i,0)*b[0];
    for (std::size_t j=1; j < columns; ++j)
      r[i] += a(i,j)*b[j];
  }
}

/** \brief multiply a DynamicVector by a DynamicMatrix */
inline void multiply(DynamicVector const& b,
                     DynamicMatrix const& a,
                     DynamicVector& r)
{
  std::size_t rows = a.getRows();
  std::size_t columns = a.getColumns();
  r.setSize(columns);
  for (std::size_t j=0; j < columns; ++j)
  {
    r[j] = a(0,j)*b[0];
    for (std::size_t i=1; i < rows; ++i)
      r[j] += b[i]*a(i,j);
  }
}

/** \brief multiply two DynamicMatrix objects */
inline void multiply(DynamicMatrix const& a,
                     DynamicMatrix const& b,
                     DynamicMatrix& r)
{
  std::size_t rows = a.getRows();
  std::size_t columns = b.getColumns();
  std::size_t depth = b.getRows();
  r.setSize(rows,columns);
  for (std::size_t i=0; i < rows; ++i)
  for (std::size_t j=0; j < columns; ++j)
  {
    r(i,j) = a(i,0)*b(0,j);
    for (std::size_t k=1; k < depth; ++k)
      r(i,j) += a(i,k)*b(k,j);
  }
}

/** \brief get the transpose of a DynamicMatrix */
inline void transpose(DynamicMatrix const& a,
                      DynamicMatrix& r)
{
  std::size_t rows = a.getRows();
  std::size_t columns = a.getColumns();
  r.setSize(columns,rows);
  for (std::size_t i=0; i < rows; ++i)
  for (std::size_t j=0; j < columns; ++j)
    r(j,i) = a(i,j);
}

/** \brief convert an apf::Matrix into an apf::DynamicMatrix */
template <std::size_t N, std::size_t M>
inline DynamicMatrix fromMatrix(Matrix<N,M> other)
{
  DynamicMatrix result(N,M);
  for(std::size_t ii = 0; ii < N; ii++)
    for(std::size_t jj = 0; jj < M; jj++)
      result(ii,jj) = other[ii][jj];
  return result;
}

}//namespace apf

/** \brief write an apf::DynamicMatrix to a C++ stream */
std::ostream& operator<<(std::ostream& s, apf::DynamicMatrix const& A);

#endif
