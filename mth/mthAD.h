#ifndef MTH_AD_H
#define MTH_AD_H

#include <cmath>

#include "../can/canArray.h"
namespace mth {

/** \brief forward automatic differentiation variable */

template <unsigned int N=0>
class AD
{
  protected:
    /** \brief the variable value */
    double x_;
    /** \brief the derivative array */
    double dx_[N];
  public:
    /** \brief the number of derivatives */
    enum { degree = N };
    /** \brief default constructor */
    AD():x_(0.) {zero();}
    /** \brief construct from a double */
    AD(double x):x_(x) {zero();}
    /** \brief copy constructor */
    AD(AD<N> const& other) {copy(other);}
    /** \brief get the size of the derivative array */
    unsigned int size() const {return N;}
    /** \brief set as the ith variable of N */
    void diff(unsigned int i)
    {
      zero();
      dx_[i] = 1.;
    }
    /** \brief get the value of the variable (mutable) */
    double& val() {return x_;}
    /** \brief get the value of the variable (immutable) */
    const double& val() const {return x_;}
    /** \brief get the ith derivative value (mutable) */
    double& dx(unsigned int i) {return dx_[i];}
    /** \brief get the ith derivative value (immutable) */
    const double& dx(unsigned int i) const {return dx_[i];}
    /** \brief assignment to a double */
    AD<N>& operator=(double other)
    {
      x_ = other;
      zero();
      return *this;
    }
    /** \brief assignment to another AD variable */
    AD<N>& operator=(AD<N> const& other)
    {
      copy(other);
      return *this;
    }
    /** \brief addition assignment with a double */
    AD<N>& operator+=(double other)
    {
      x_ += other;
      return *this;
    }
    /** \brief addition assignment with another AD variable */
    AD<N>& operator+=(AD<N> const& other)
    {
      x_ += other.x_;
      for (unsigned int i=0; i < N; ++i)
        dx_[i] += other.dx_[i];
      return *this;
    }
    /** \brief subtraction assignment with a double */
    AD<N>& operator-=(double other)
    {
      x_ -= other;
      return *this;
    }
    /** \brief subtraction assignment with another AD variable */
    AD<N>& operator-=(AD<N> const& other)
    {
      x_ -= other.x_;
      for (unsigned int i=0; i < N; ++i)
        dx_[i] -= other.dx_[i];
      return *this;
    }
    /** \brief multiplication assignment with a double */
    AD<N>& operator*=(double other)
    {
      x_ *= other;
      for (unsigned int i=0; i < N; ++i)
        dx_[i] *= other;
      return *this;
    }
    /** \brief multiplication assignment with another AD variable */
    AD<N>& operator*=(AD<N> const& other)
    {
      x_ *= other.x_;
      for (unsigned int i=0; i < N; ++i)
        dx_[i] = dx_[i]*other.x_ + x_*other.dx_[i];
    }
    /** \brief division assignment with a double */
    AD<N>& operator/=(double other)
    {
      x_ /= other;
      for (unsigned int i=0; i < N; ++i)
        dx_[i] /= other;
      return *this;
    }
    /** \brief division assignment with another AD variable */
    AD<N>& operator/=(AD<N> const& other)
    {
      x_ /= other.x_;
      for (unsigned int i=0; i < N; ++i)
        dx_[i] = (dx_[i]*other.x_ - x_*other.dx_[i]) / (other.x_*other.x_);
    }
  private:
    void zero()
    {
      for (unsigned int i=0; i < N; ++i)
        dx_[i] = 0.;
    }
    void copy(AD<N> const& other)
    {
      x_ = other.x_;
      for (unsigned int i=0; i < N; ++i)
        dx_[i] = other.dx_[i];
    }
};
/** \brief forward automatic differentiation variable with dynamic variable array */
template<>
class AD<0>
{
  protected:
    /** \brief the variable value */
    double x_;
    /** \brief the dynamic derivative array*/
    can::Array<double> dx_;
  public:
    const static double _zero_;
    /** \brief default constructor */
    AD():x_(0), dx_() {zero();}
    /** \brief constructs from a double*/
    AD(double val):x_(val), dx_() {zero();}
    /** \brief constructs from other AD*/
    AD(AD<0> const& other) {copy(other);}
    unsigned size() { return dx_.size();}
    const unsigned size() const { return dx_.size();}
    void diff(unsigned int i)
    {
      if(dx_.size() < i + 1)
        dx_.resize_copy(i + 1);
      zero();
      dx_[i] = 1;
    }
    /** \brief get the value of the variable (mutable) */
    double& val() {return x_;}
    /** \brief get the value of the variable (immutable) */
    const double& val() const {return x_;}
    /** \brief get the ith derivative value (mutable) */
    double& dx(unsigned int i)
    {
      if(i >= size())
      {
        dx_.resize_copy(i);
        for(unsigned int x = size(); x < i; x++)
          dx_[x] = 0.;
      }
      return dx_[i];
    }
    /** \brief get the ith deriative value (immutable) */
    const double& dx(unsigned int i) const {
      if(i >= size())
        return _zero_;
      return dx_[i];
    }
    /** \brief assignment to a double */
    AD<0>& operator=(double other)
    {
      x_ = other;
      zero();
      return *this;
    }
    /** \brief assignment to another AD variable */
    AD<0>& operator=(AD const& other)
    {
      resize(other.size());
      copy(other);
      return *this;
    }
    /** \brief addition assignment with a double */
    AD<0>& operator+=(double other)
    {
      x_ += other;
      return *this;
    }
    /** \brief addition assignment with another AD variable */
    AD<0>& operator+=(AD<0> const& other)
    {
      resize(other.size());
      x_ += other.x_;
      for (unsigned int i=0; i < size(); ++i)
        dx_[i] += other.dx_[i];
      return *this;
    }
    /** \brief subtraction assignment with a double */
    AD<0>& operator-=(double other)
    {
      x_ -= other;
      return *this;
    }
    /** \brief subtraction assignment with another AD variable */
    AD<0>& operator-=(AD<0> const& other)
    {
      resize(other.size());
      x_ -= other.x_;
      for (unsigned int i=0; i < size(); ++i)
        dx_[i] -= other.dx_[i];
      return *this;
    }
    /** \brief multiplication assignment with a double */
    AD<0>& operator*=(double other)
    {
      x_ *= other;
      for (unsigned int i=0; i < size(); ++i)
        dx_[i] *= other;
      return *this;
    }
    /** \brief multiplication assignment with another AD variable */
    AD<0>& operator*=(AD<0> const& other)
    {
      resize(other.size());
      x_ *= other.x_;
      for (unsigned int i=0; i < size(); ++i)
        dx_[i] = dx_[i]*other.x_ + x_*other.dx_[i];
      return *this;
    }
    /** \brief division assignment with a double */
    AD<0>& operator/=(double other)
    {
      x_ /= other;
      for (unsigned int i=0; i < size(); ++i)
        dx_[i] /= other;
      return *this;
    }
    /** \brief division assignment with another AD variable */
    AD<0>& operator/=(AD<0> const& other)
    {
      resize(other.size());
      x_ /= other.x_;
      for (unsigned int i=0; i < size(); ++i)
        dx_[i] = (dx_[i]*other.x_ - x_*other.dx_[i]) / (other.x_*other.x_);
      return *this;
    }
    void resize(unsigned int n)
    {
      dx_.resize(n);
    }
  private:
    void zero()
    {
      for (unsigned int i=0; i < size(); ++i)
        dx_[i] = 0.;
    }
    void copy(AD<0> const& other)
    {
      if(size() != other.size())
        dx_.resize(other.size());
      x_ = other.x_;
      for (unsigned int i=0; i < size(); ++i)
        dx_[i] = other.dx_[i];
    }
};
const double mth::AD<0>::_zero_ = 0;
/**********************
  * UNARY OPERATIONS *
***********************/

/** \brief unary subtraction */
template <unsigned int N>
AD<N> operator-(AD<N> const& A)
{
  AD<N> tmp;
  tmp -= A;
  return tmp;
}

/** \brief unary subtraction for dynamic AD*/
AD<0> operator-(AD<0> const& A)
{
  AD<0> tmp;
  tmp.resize(A.size());
  tmp -= A;
  return tmp;
}

/***********************
  * BINARY OPERATIONS *
************************/

/** \brief binary addition between a double and an AD variable */
template <unsigned int N>
AD<N> operator+(double L, AD<N> const& R)
{
  AD<N> tmp;
  for (unsigned int i=0; i < N; ++i)
    tmp.dx(i) = R.dx(i);
  tmp.val() = L + R.val();
  return tmp;
}
/** \brief binary addition between a double and a dynamic AD variable */
AD<0> operator+(double L, AD<0>& R)
{
  AD<0> tmp;
  tmp.resize(R.size());
  for (unsigned int i=0; i < tmp.size(); ++i)
    tmp.dx(i) = R.dx(i);
  tmp.val() = L + R.val();
  return tmp;
}

/** \brief binary addition between an AD variable and a double */
template <unsigned int N>
AD<N> operator+(AD<N> const& L, double R)
{
  AD<N> tmp;
  for (unsigned int i=0; i < N; ++i)
    tmp.dx(i) = L.dx(i);
  tmp.val() = L.val() + R;
  return tmp;
}

/** \brief binary addition between an AD variable and a double */
AD<0> operator+(AD<0> const& L, double R)
{
  AD<0> tmp;
  tmp.resize(L.size());
  for (unsigned int i=0; i < tmp.size(); ++i)
    tmp.dx(i) = L.dx(i);
  tmp.val() = L.val() + R;
  return tmp;
}
/** \brief binary addition between two AD variables */

template <unsigned int N>
AD<N> operator+(AD<N> const& L, AD<N> const& R)
{
  AD<N> tmp;
  for (unsigned int i=0; i < N; ++i)
    tmp.dx(i) = L.dx(i) + R.dx(i);
  tmp.val() = L.val() + R.val();
  return tmp;
}
/** \brief binary addition between two dynamic AD variables */
AD<0> operator+(AD<0> const& L, AD<0> const& R)
{
  AD<0> tmp;
  unsigned int max = L.size() > R.size() ? L.size() : R.size();
  tmp.resize(max);
  for (unsigned int i=0; i < tmp.size(); ++i)
    tmp.dx(i) = L.dx(i) + R.dx(i);
  tmp.val() = L.val() + R.val();
  return tmp;
}

/** \brief binary subtraction between a double and an AD variable */
template <unsigned int N>
AD<N> operator-(double L, AD<N> const& R)
{
  AD<N> tmp;
  for (unsigned int i=0; i < N; ++i)
    tmp.dx(i) = -R.dx(i);
  tmp.val() = L - R.val();
  return tmp;
}

/** \brief binary subtraction between a double and a dynamic AD variable */
AD<0> operator-(double L, AD<0> const& R)
{
  AD<0> tmp;
  tmp.resize(R.size());
  for (unsigned int i=0; i < tmp.size(); ++i)
    tmp.dx(i) = -R.dx(i);
  tmp.val() = L - R.val();
  return tmp;
}

/** \brief binary subtraction between an AD variable and a double */
template <unsigned int N>
AD<N> operator-(AD<N> const& L, double R)
{
  AD<N> tmp;
  for (unsigned int i=0; i < N; ++i)
    tmp.dx(i) = L.dx(i);
  tmp.val() = L.val() - R;
  return tmp;
}

/** \brief binary subtraction between a dynamic AD variable and a double */
AD<0> operator-(AD<0> const& L, double R)
{
  AD<0> tmp;
  tmp.resize(L.size());
  for (unsigned int i=0; i < tmp.size(); ++i)
    tmp.dx(i) = L.dx(i);
  tmp.val() = L.val() - R;
  return tmp;
}

/** \brief binary subtraction between two AD variables */
template <unsigned int N>
AD<N> operator-(AD<N> const& L, AD<N> const& R)
{
  AD<N> tmp;
  for (unsigned int i=0; i < N; ++i)
    tmp.dx(i) = L.dx(i) - R.dx(i);
  tmp.val() = L.val() - R.val();
  return tmp;
}
/**brief binary subtraction between two dynamic AD variables */
AD<0> operator-(AD<0> const& L, AD<0> const& R)
{
  AD<0> tmp;
  unsigned int max = L.size() > R.size() ? L.size() : R.size();
  tmp.resize(max);
  for (unsigned int i=0; i < tmp.size(); ++i)
    tmp.dx(i) = L.dx(i) - R.dx(i);
  tmp.val() = L.val() - R.val();
  return tmp;
}

/** \brief binary multiplication between a double and an AD variable */
template <unsigned int N>
AD<N> operator*(double L, AD<N> const& R)
{
  AD<N> tmp;
  for (unsigned int i=0; i < N; ++i)
    tmp.dx(i) = L * R.dx(i);
  tmp.val() = L * R.val();
  return tmp;
}

/** \brief binary multiplication between a double and a dynamic AD variable */
AD<0> operator*(double L, AD<0> const& R)
{
  AD<0> tmp;
  tmp.resize(R.size());
  for (unsigned int i=0; i < tmp.size(); ++i)
    tmp.dx(i) = L * R.dx(i);
  tmp.val() = L * R.val();
  return tmp;
}

/** \brief binary multiplication between an AD variable and a double */
template <unsigned int N>
AD<N> operator*(AD<N> const& L, double R)
{
  AD<N> tmp;
  for (unsigned int i=0; i < N; ++i)
    tmp.dx(i) = L.dx(i) * R;
  tmp.val() = L.val() * R;
  return tmp;
}

/** \brief binary multiplication between an AD variable and a double */
AD<0> operator*(AD<0> const& L, double R)
{
  AD<0> tmp;
  tmp.resize(L.size());
  for (unsigned int i=0; i < tmp.size(); ++i)
    tmp.dx(i) = L.dx(i) * R;
  tmp.val() = L.val() * R;
  return tmp;
}

/** \brief binary multiplication between two AD variables */
template <unsigned int N>
AD<N> operator*(AD<N> const& L, AD<N> const& R)
{
  AD<N> tmp;
  for (unsigned int i=0; i < N; ++i)
    tmp.dx(i) = L.dx(i) * R.val() + L.val() * R.dx(i);
  tmp.val() = L.val() * R.val();
  return tmp;
}

/** \brief binary multiplication between two dynamic AD variables */
AD<0> operator*(AD<0> const& L, AD<0> const& R)
{
  AD<0> tmp;
  unsigned int max = L.size() > R.size() ? L.size() : R.size();
  tmp.resize(max);
  for (unsigned int i=0; i < tmp.size(); ++i)
    tmp.dx(i) = L.dx(i) * R.val() + L.val() * R.dx(i);
  tmp.val() = L.val() * R.val();
  return tmp;
}

/** \brief binary division between a double and an AD variable */
template <unsigned int N>
AD<N> operator/(double L, AD<N> const& R)
{
  AD<N> tmp;
  for (unsigned int i=0; i < N; ++i)
    tmp.dx(i) = ( -L*R.dx(i) ) / (R.val()*R.val());
  tmp.val() = L / R.val();
  return tmp;
}

/** \brief binary division between a double and a dynamic AD variable */
AD<0> operator/(double L, AD<0> const& R)
{
  AD<0> tmp;
  tmp.resize(R.size());
  for (unsigned int i=0; i < tmp.size(); ++i)
    tmp.dx(i) = ( -L*R.dx(i) ) / (R.val()*R.val());
  tmp.val() = L / R.val();
  return tmp;
}

/** \brief binary division between an AD variable and a double */
template <unsigned int N>
AD<N> operator/(AD<N> const& L, double R)
{
  AD<N> tmp;
  for (unsigned int i=0; i < N; ++i)
    tmp.dx(i) = L.dx(i) / R;
  tmp.val() = L.val() / R;
  return tmp;
}

/** \brief binary division between a dynamic AD variable and a double */
AD<0> operator/(AD<0> const& L, double R)
{
  AD<0> tmp;
  tmp.resize(L.size());
  for (unsigned int i=0; i < L.size(); ++i)
    tmp.dx(i) = L.dx(i) / R;
  tmp.val() = L.val() / R;
  return tmp;
}
/** \brief binary division between two AD variables */
template <unsigned int N>
AD<N> operator/(AD<N> const& L, AD<N> const& R)
{
  AD<N> tmp;
  for (unsigned int i=0; i < N; ++i)
    tmp.dx(i) = (L.dx(i) * R.val() - L.val() * R.dx(i) ) / (R.val() * R.val());
  tmp.val() = L.val() / R.val();
  return tmp;
}

/** \brief binary division between two  dynamic AD variables */
AD<0> operator/(AD<0> const& L, AD<0> const& R)
{
  AD<0> tmp;
  unsigned int max = L.size() > R.size() ? L.size() : R.size();
  tmp.resize(max);
  for (unsigned int i=0; i < tmp.size(); ++i)
    tmp.dx(i) = (L.dx(i) * R.val() - L.val() * R.dx(i) ) / (R.val() * R.val());
  tmp.val() = L.val() / R.val();
  return tmp;
}

/********************
  * FANCY FUNCIONS *
*********************/

/** \brief exponent of an AD variable */
template <unsigned int N>
AD<N> exp(AD<N> const& A)
{
  AD<N> tmp(std::exp(A.val()));
  for (unsigned int i=0; i < N; ++i)
    tmp.dx(i) = A.dx(i) * std::exp(A.val());
  return tmp;
}
/** \brief exponent of a dynamic AD variable */
AD<0> exp(AD<0> const& A)
{
  AD<0> tmp;
  tmp.resize(A.size());
  tmp.val() = std::exp(A.val());
  for (unsigned int i=0; i < A.size(); ++i)
    tmp.dx(i) = A.dx(i) * std::exp(A.val());
  return tmp;
}

/** \brief logarithm of an AD variable */
template <unsigned int N>
AD<N> log(AD<N> const& A)
{
  AD<N> tmp(std::log(A.val()));
  for (unsigned int i=0; i < N; ++i)
    tmp.dx(i) = A.dx(i) / A.val();
  return tmp;
}

/** \brief logarithm of an AD variable */
template <unsigned int N>
AD<0> log(AD<0> const& A)
{
  AD<0> tmp;
  tmp.resize(A.size());
  tmp.val() = std::exp(A.val());
  for (unsigned int i=0; i < N; ++i)
    tmp.dx(i) = A.dx(i) / A.val();
  return tmp;
}

/** \brief AD variable raised to an integer power */
template <unsigned int N>
AD<N> pow(AD<N> const& A, int e)
{
  AD<N> tmp(std::pow(A.val(), (double)e));
  for (unsigned int i=0; i < N; ++i)
    tmp.dx(i) = e*A.dx(i)*std::pow(A.val(), (double)e-1.);
  return tmp;
}
/** \brief dynamic AD variable raised to a integer power */
AD<0> pow(AD<0> const& A, int e)
{
  AD<0> tmp;
  tmp.resize(A.size());
  tmp.val() = std::pow(A.val(), e);
  for (unsigned int i=0; i < tmp.size(); ++i)
    tmp.dx(i) = e*A.dx(i)*std::pow(A.val(), (double)e-1.);
  return tmp;
}

/** \brief AD variable raised to a double power */
template <unsigned int N>
AD<N> pow(AD<N> const& A, double e)
{
  AD<N> tmp(std::pow(A.val(), e));
  for (unsigned int i=0; i < N; ++i)
    tmp.dx(i) = e*A.dx(i)*std::pow(A.val(), e-1.);
  return tmp;
}

/** \brief dynamic AD variable raised to a double power */
AD<0> pow(AD<0> const& A, double e)
{
  AD<0> tmp;
  tmp.resize(A.size());
  tmp.val() = std::pow(A.val(), e);
  for (unsigned int i=0; i < tmp.size(); ++i)
    tmp.dx(i) = e*A.dx(i)*std::pow(A.val(), e-1.);
  return tmp;
}
/** \brief AD variable raised to an AD variable power */
template <unsigned int N>
AD<N> pow(AD<N> const& A, AD<N> const& e)
{
  AD<N> tmp(std::pow(A.val(), e.val()));
  for (unsigned int i=0; i < N; ++i)
    tmp.dx(i) = e.dx(i) * std::log(A.val()) * std::pow(A.val(), e.val()) +
      e.val() * A.dx(i) * std::pow(A.val(), e.val()-1.);
  return tmp;
}

/** \brief AD variable raised to an AD variable power */
AD<0> pow(AD<0> const& A, AD<0> const& e)
{
  AD<0> tmp;
  unsigned int max = A.size() > e.size() ? A.size() : e.size();
  tmp.resize(max);
  tmp.val() = std::pow(A.val(), e.val());
  for (unsigned int i=0; i < tmp.size(); ++i)
    tmp.dx(i) = e.dx(i) * std::log(A.val()) * std::pow(A.val(), e.val()) +
      e.val() * A.dx(i) * std::pow(A.val(), e.val()-1.);
  return tmp;
}
/** \brief square root of an AD variable */
template <unsigned int N>
AD<N> sqrt(AD<N> const& A)
{
  AD<N> tmp(std::sqrt(A.val()));
  for (unsigned int i=0; i < N; ++i)
    tmp.dx(i) = A.dx(i) / (2. * std::sqrt(A.val()));
  return tmp;
}
/** \brief square root of an AD variable */
AD<0> sqrt(AD<0> const& A)
{
  AD<0> tmp;
  tmp.resize(A.size());
  tmp.val() = std::sqrt(A.val());
  for (unsigned int i=0; i < tmp.size(); ++i)
    tmp.dx(i) = A.dx(i) / (2. * std::sqrt(A.val()));
  return tmp;
}


/** \brief absolute value of an AD variable */
template <unsigned int N>
AD<N> abs(AD<N> const& A)
{
  int sign = A.val() > 0 ? 1 : 0;
  if (sign) return A;
  else return (-A);
}

/**************************
  * COMPARISON OPERATORS *
***************************/

/** \brief double less than an AD variable */
template <unsigned int N>
bool operator<(double L, AD<N> const& R)
{
  return L < R.val();
}

/** \brief AD variable less than a double */
template <unsigned int N>
bool operator<(AD<N> const& R, double L)
{
  return R.val() < L;
}

/** \brief AD variable less than an AD variable */
template <unsigned int N>
bool operator<(AD<N> const& R, AD<N> const& L)
{
  return R.val() < L.val();
}

/** \brief double less than or equal to an AD variable */
template <unsigned int N>
bool operator<=(double L, AD<N> const& R)
{
  return L <= R.val();
}

/** \brief AD variable less than or equal to a double */
template <unsigned int N>
bool operator<=(AD<N> const& R, double L)
{
  return R.val() <= L;
}

/** \brief AD variable less than or equal to an AD variable */
template <unsigned int N>
bool operator<=(AD<N> const& R, AD<N> const& L)
{
  return R.val() <= L.val();
}

/** \brief double greater than an AD variable */
template <unsigned int N>
bool operator>(double L, AD<N> const& R)
{
  return L > R.val();
}

/** \brief AD variable greater than a double */
template <unsigned int N>
bool operator>(AD<N> const& R, double L)
{
  return R.val() > L;
}

/** \brief AD variable greater than an AD variable */
template <unsigned int N>
bool operator>(AD<N> const& R, AD<N> const& L)
{
  return R.val() > L.val();
}

/** \brief double greater than or equal to an AD variable */
template <unsigned int N>
bool operator>=(double L, AD<N> const& R)
{
  return L >= R.val();
}

/** \brief AD variable greater than or equal to a double */
template <unsigned int N>
bool operator>=(AD<N> const& R, double L)
{
  return R.val() >= L;
}

/** \brief AD variable greater than or equal to an AD variable */
template <unsigned int N>
bool operator>=(AD<N> const& R, AD<N> const& L)
{
  return R.val() >= L.val();
}
}
#endif