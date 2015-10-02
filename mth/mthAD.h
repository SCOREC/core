#ifndef MTH_AD_H
#define MTH_AD_H

#include <cmath>

namespace mth {

/** \brief forward automatic differentiation variable */
template <unsigned int N>
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
    /** \brief assigment to a double */
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
      for (int i=0; i < N; ++i)
        dx_[i] = 0.;
    }
    void copy(AD<N> const& other)
    {
      x_ = other.x_;
      for (int i=0; i < N; ++i)
        dx_[i] = other.dx_[i];
    }
};

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

/** \brief logarithm of an AD variable */
template <unsigned int N>
AD<N> log(AD<N> const& A)
{
  AD<N> tmp(std::log(A.val()));
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

/** \brief AD variable raised to a double power */
template <unsigned int N>
AD<N> pow(AD<N> const& A, double e)
{
  AD<N> tmp(std::pow(A.val(), e));
  for (unsigned int i=0; i < N; ++i)
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
