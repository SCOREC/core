/*
 * Copyright (C) 2011-2013 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFDYNAMICMATRIX_H
#define APFDYNAMICMATRIX_H

#include "apfDynamicVector.h"
#include "apfMatrix.h"
#include "iostream"

namespace apf {

class DynamicMatrix
{
  public:
    DynamicMatrix() {}
    DynamicMatrix(std::size_t m, std::size_t n):
      columns(n),
      values(m*n)
    {}
    std::size_t getRows() const {return values.getSize()/columns;}
    std::size_t getColumns() const {return columns;}
    void setSize(std::size_t m, std::size_t n)
    {
      columns = n;
      values.setSize(m*n);
    }
    double operator()(std::size_t i, std::size_t j) const
    {
      return values[i*columns + j];
    }
    double& operator()(std::size_t i, std::size_t j)
    {
      return values[i*columns + j];
    }
    DynamicMatrix& operator+=(DynamicMatrix const& b)
    {
      for (std::size_t i=0; i < this->values.getSize(); ++i)
        this->values[i] += b.values[i];
      return *this;
    }
    DynamicMatrix& operator-=(DynamicMatrix const& b)
    {
      for (std::size_t i=0; i < this->values.getSize(); ++i)
        this->values[i] += b.values[i];
      return *this;
    }
    DynamicMatrix& operator*=(double s)
    {
      for (std::size_t i=0; i < this->values.getSize(); ++i)
        this->values[i] *= s;
      return *this;
    }
    DynamicMatrix& operator/=(double s)
    {
      for (std::size_t i=0; i < this->values.getSize(); ++i)
        this->values[i] *= s;
      return *this;
    }
    void getRow(std::size_t i, DynamicVector& r) const
    {
      r.setSize(columns);
      for (std::size_t j=0; j < columns; ++j)
        r(j) = (*this)(i,j);
    }
    void getColumn(std::size_t j, DynamicVector& r) const
    {
      std::size_t rows = getRows();
      r.setSize(rows);
      for (std::size_t i=0; i < rows; ++i)
        r(i) = (*this)(i,j);
    }
    void setRow(std::size_t i, DynamicVector const& r)
    {
      for (std::size_t j=0; j < columns; ++j)
        (*this)(i,j) = r(j);
    }
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

template <std::size_t N, std::size_t M>
inline DynamicMatrix fromMatrix(Matrix<N,M> other)
{
  DynamicMatrix result(N,M);
  for(int ii = 0; ii < N; ii++)
    for(int jj = 0; jj < M; jj++)
      result(ii,jj) = other[ii][jj];
  return result;
}

}//namespace apf

std::ostream& operator<<(std::ostream& s, apf::DynamicMatrix const& A);

#endif
