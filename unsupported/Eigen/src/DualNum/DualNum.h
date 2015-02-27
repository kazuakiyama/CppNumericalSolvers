// -*-c++-*-
// This file is part of Eigen, a lightweight C++ template library
// for linear algebra. 
//
// The MIT License (MIT)
// 
// Copyright (c) 2006 Jeffrey A. Fike
// Copyright (C) 2015 Michael Tesch tesch a tum de
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
//  of this software and associated documentation files (the "Software"), to deal
//  in the Software without restriction, including without limitation the rights
//  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//  copies of the Software, and to permit persons to whom the Software is
//  furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
//  all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
//  THE SOFTWARE.
// 
//
#ifndef EIGEN_DUALNUM_H
#define EIGEN_DUALNUM_H

#include <cmath>

namespace Eigen {

/* dual-numbers implementation for calculation of differentials */
template <class _Tp>
class DualNum {

private:
  _Tp _f0, _f1;
  
public:

  typedef _Tp value_type;
  // todo? - scalar type of value_type ?

  inline DualNum() : _f0(), _f1() {}
  inline DualNum(const _Tp & re) : _f0(re), _f1() {}
  inline DualNum(const _Tp & re, const _Tp & du) : _f0(re), _f1(du) {}

  inline _Tp &
  realpart() { return _f0; }
  
  inline const _Tp &
  realpart() const { return _f0; }
  
  inline void
  realpart(_Tp f0) { _f0 = f0; }
  
  inline _Tp &
  epart() { return _f1; }
  
  inline const _Tp &
  epart() const { return _f1; }
  
  inline void
  epart(_Tp f1) { _f1 = f1; }

  // basic ops
  DualNum<_Tp> operator+() const {
    return *this;
  }

  DualNum<_Tp> operator+(const DualNum<_Tp> rhs) const {
    return DualNum<_Tp>(_f0 + rhs._f0, _f1 + rhs._f1);
  }

  DualNum<_Tp> operator-(const DualNum<_Tp> rhs) const {
    return DualNum<_Tp>(_f0 - rhs._f0, _f1 - rhs._f1);
  }

  DualNum<_Tp> operator-() const {
    return DualNum<_Tp>(-_f0, -_f1);
  }

  DualNum<_Tp> operator*(const DualNum<_Tp> rhs) const {
    return DualNum<_Tp>(_f0 * rhs._f0, _f0 * rhs._f1 + _f1 * rhs._f0);
  }

  DualNum<_Tp> & operator+=(DualNum<_Tp> rhs) {
    _f0 += rhs._f0;
    _f1 += rhs._f1;
    return *this;
  }

  DualNum<_Tp> & operator-=(DualNum<_Tp> rhs) {
    _f0 -= rhs._f0;
    _f1 -= rhs._f1;
    return *this;
  }

  DualNum<_Tp> & operator*=(DualNum<_Tp> rhs) {
    _Tp tf0, tf1;
    tf0 = _f0;
    tf1 = _f1;
    _f0 = tf0 * rhs._f0;
    _f1 = tf0 * rhs._f1 + tf1 * rhs._f0;
    return *this;
  }

  DualNum<_Tp> & operator*=(_Tp rhs) {
    _f0 *= rhs;
    _f1 *= rhs;
    return *this;
  }

  DualNum<_Tp> & operator/=(_Tp rhs) {
    _f0 /= rhs;
    _f1 /= rhs;
  }

};

// function defs
template <class eeT>
inline eeT
realpart(const DualNum<eeT> & d)
{
  return d.realpart();
}

template <class eeT>
inline eeT
epart(const DualNum<eeT> & d)
{
  return d.epart();
}

// basic ops
template <class eeT>
inline DualNum<eeT>
operator+(const eeT lhs, const DualNum<eeT> rhs)
{
  return DualNum<eeT>(lhs) += rhs;
}

template <class eeT>
inline DualNum<eeT>
operator-(const eeT lhs, const DualNum<eeT> rhs)
{
  return DualNum<eeT>(lhs) -= rhs;
}

template <class eeT>
inline DualNum<eeT>
operator*(const eeT lhs, const DualNum<eeT> rhs)
{
  return DualNum<eeT>(lhs) *= rhs;
}

template <class eeT>
inline DualNum<eeT>
operator/(const DualNum<eeT> lhs, const DualNum<eeT> rhs)
{
  DualNum<eeT> temp = lhs;
  return lhs /= rhs;
}

template <class eeT>
inline DualNum<eeT>
operator/(const eeT lhs, const DualNum<eeT> rhs)
{
  return DualNum<eeT>(lhs) /= rhs;
}

template <class eeT>
inline DualNum<eeT>
operator/(const DualNum<eeT> lhs, const eeT rhs)
{
  DualNum<eeT> temp = lhs;
  return temp /= rhs;
}

// math.h functions
template <class _Tp>
DualNum<_Tp>
pow(DualNum<_Tp> x, _Tp a)
{
  DualNum<_Tp> temp;
  _Tp deriv, xval, tol;
  xval = x.realpart();
  tol = _Tp(1e-15); // TODO- should use type traits of _Tp instead of 1e-15
  if (std::abs(xval) < std::abs(tol)) {
    xval = x.realpart() / (std::abs(x.realpart()) / _Tp(1e-15));
    //if (xval >= 0)
    //  xval = tol;
    //if (xval < 0)
    //  xval = -tol;
  }
  deriv = a * std::pow(xval, (a - _Tp(1.0)));
  temp.realpart() = std::pow(x.realpart(), a);  //Use actual x value, only use tol for derivs
  temp.epart() = x.epart() * deriv;
  return temp;
}

//template <class eeT> DualNum<eeT> pow(DualNum<eeT> x, DualNum<eeT> a) {
//
//}

template <class eeT>
DualNum<eeT> exp(DualNum<eeT> x)
{
  eeT deriv;
  deriv = std::exp(x.realpart());
  return DualNum<eeT>(deriv, deriv * x.epart());
}

template <class eeT>
DualNum<eeT>
log(DualNum<eeT> x)
{
  eeT deriv1;
  deriv1 = x.epart() / x.realpart();
  return DualNum<eeT>(log(x.realpart()), deriv1);
}

template <class eeT>
DualNum<eeT>
sin(DualNum<eeT> x)
{
  DualNum<eeT> temp;
  eeT funval, deriv;
  funval = sin(x.realpart());
  deriv = cos(x.realpart());
  temp.realpart() = funval;
  temp.epart() = deriv * x.epart();
  return temp;
}

template <class eeT>
DualNum<eeT>
cos(DualNum<eeT> x)
{
  DualNum<eeT> temp;
  eeT funval, deriv;
  funval = cos(x.realpart());
  deriv = -sin(x.realpart());
  temp.realpart() = funval;
  temp.epart() = deriv * x.epart();
  return temp;
}

template <class eeT>
DualNum<eeT>
tan(DualNum<eeT> x)
{
  DualNum<eeT> temp;
  eeT funval, deriv;
  funval = tan(x.realpart());
  deriv  = funval*funval + 1.0;
  temp.realpart() = funval;
  temp.epart() = deriv*x.epart();
  return temp;
}

template <class eeT>
DualNum<eeT>
asin(DualNum<eeT> x)
{
  DualNum<eeT> temp;
  eeT funval, deriv1, deriv;
  funval = asin(x.realpart());
  deriv1 = 1.0-x.realpart()*x.realpart();
  deriv = 1.0/sqrt(deriv1);
  temp.realpart() = funval;
  temp.epart() = deriv*x.epart();
  return temp;
}

#if 0
template <class eeT>
DualNum<eeT>
acos(DualNum<eeT> x)
{
  DualNum<eeT> temp;
  eeT funval, deriv1, deriv;
  funval = acos(x.realpart());
  deriv1 = 1.0 - x.realpart() * x.realpart();
  deriv = 1.0 / sqrt(deriv1);
  temp.realpart() = funval;
  temp.epart() = deriv*x.epart();
  return temp;
}
#endif

template <class eeT>
DualNum<eeT>
atan(DualNum<eeT> x)
{
  DualNum<eeT> temp;
  eeT funval, deriv1, deriv;
  funval = atan(x.realpart());
  deriv1 = 1.0 + x.realpart() * x.realpart();
  deriv = 1.0 / deriv1;
  temp.realpart() = funval;
  temp.epart() = deriv * x.epart();
  return temp;
}

template <class eeT>
DualNum<eeT>
atan2(DualNum<eeT> y, DualNum<eeT> x)
{
  DualNum<eeT> temp;
  eeT funval, deriv1, deriv;
  funval = atan2(y.realpart(), x.realpart());
  // unsure from here on...
  deriv1 = 1.0 + x.realpart() * x.realpart();
  deriv = 1.0 / deriv1;
  temp.realpart() = funval;
  temp.epart() = deriv * x.epart();
  return temp;
}

template <class eeT>
DualNum<eeT>
sqrt(DualNum<eeT> x)
{
  return pow(x, (eeT)0.5);
}

template <class eeT>
DualNum<eeT>
max(DualNum<eeT> x1, DualNum<eeT> x2)
{
  return x1.realpart() >= x2.realpart() ? x1 : x2;
}

template <class eeT>
DualNum<eeT>
max(DualNum<eeT> x1, eeT x2)
{
  return x1.realpart() >= x2 ? x1 : DualNum<eeT>(x2);
}

template <class eeT>
DualNum<eeT>
max(eeT x1, DualNum<eeT> x2)
{
  return x1 >= x2.realpart() ? DualNum<eeT>(x1) : x2;
}

template <class eeT>
DualNum<eeT>
min(DualNum<eeT> x1, DualNum<eeT> x2)
{
  return x1.realpart() <= x2.realpart() ? x1 : x2;
}

template <class eeT>
DualNum<eeT>
min(DualNum<eeT> x1, eeT x2)
{
  return x1.realpart() <= x2 ? x1 : DualNum<eeT>(x2);
}

template <class eeT>
DualNum<eeT>
min(eeT x1, DualNum<eeT> x2)
{
  return x1 <= x2.realpart() ? DualNum<eeT>(x1) : x2;
}

template <class eeT>
DualNum<eeT>
conj(const DualNum<eeT> & x)
{
  return DualNum<eeT>(conj(x.realpart()), conj(x.epart()));
}

template <class eeT>
DualNum<eeT>
real(const DualNum<eeT> & x)
{
  // todo - dont just make things up
  return DualNum<eeT>(real(x.realpart()), real(x.epart()));
}

template <class eeT>
DualNum<eeT>
imag(const DualNum<eeT> & x)
{
  // todo - dont just make things up
  return DualNum<eeT>(imag(x.realpart()), imag(x.epart()));
}

template <class eeT>
inline DualNum<eeT>
abs(const DualNum<eeT> & x)
{
  return std::abs(x.realpart()) == x.realpart() ? x : -x;
}

template <class eeT>
DualNum<eeT>
abs2(const DualNum<eeT> & x)
{
  return x * x;
}

/// comparison
template <class eeT>
inline bool
operator>(DualNum<eeT> lhs, DualNum<eeT> rhs)
{
  return lhs.realpart() > rhs.realpart();
}

template <class eeT>
inline bool
operator>(eeT lhs, DualNum<eeT> rhs)
{
  return lhs > rhs.realpart();
}

template <class eeT>
inline bool
operator>(DualNum<eeT> lhs, eeT rhs)
{
  return lhs.realpart() > rhs;
}

template <class eeT>
inline bool
operator>=(DualNum<eeT> lhs, DualNum<eeT> rhs)
{
  return lhs.realpart() >= rhs.realpart();
}

template <class eeT>
inline bool
operator>=(eeT lhs, DualNum<eeT> rhs)
{
  return lhs >= rhs.realpart();
}

template <class eeT>
inline bool
operator>=(DualNum<eeT> lhs, eeT rhs)
{
  return lhs.realpart() >= rhs;
}

template <class eeT>
inline bool
operator<(DualNum<eeT> lhs, DualNum<eeT> rhs)
{
  return lhs.realpart() < rhs.realpart();
}

template <class eeT>
inline bool
operator<(eeT lhs, DualNum<eeT> rhs)
{
  return lhs < rhs.realpart();
}

template <class eeT>
inline bool
operator<(DualNum<eeT> lhs, eeT rhs)
{
  return lhs.realpart() < rhs;
}

template <class eeT>
inline bool
operator<=(DualNum<eeT> lhs, DualNum<eeT> rhs)
{
  return lhs.realpart() <= rhs.realpart();
}

template <class eeT>
inline bool
operator<=(eeT lhs, DualNum<eeT> rhs)
{
  return lhs <= rhs.realpart();
}

template <class eeT>
inline bool
operator<=(DualNum<eeT> lhs, eeT rhs)
{
  return lhs.realpart() <= rhs;
}

template <class eeT>
inline bool
operator==(DualNum<eeT> lhs, DualNum<eeT> rhs)
{
  return lhs.realpart() == rhs.realpart();
}

template <class eeT>
inline bool
operator==(eeT lhs, DualNum<eeT> rhs)
{
  return lhs == rhs.realpart();
}

template <class eeT>
inline bool
operator==(DualNum<eeT> lhs, eeT rhs)
{
  return lhs.realpart() == rhs;
}

template <class eeT>
inline bool
operator!=(DualNum<eeT> lhs, DualNum<eeT> rhs)
{
  return lhs.realpart() != rhs.realpart();
}

template <class eeT>
inline bool
operator!=(eeT lhs, DualNum<eeT> rhs)
{
  return lhs != rhs.realpart();
}

template <class eeT>
inline bool
operator!=(DualNum<eeT> lhs, eeT rhs)
{
  return lhs.realpart() != rhs;
}

} // namespace Eigen

#include <Eigen/Core>
namespace Eigen {
template<typename _Scalar>
struct NumTraits<DualNum<_Scalar> > : GenericNumTraits<_Scalar>
{
  typedef _Scalar Scalar;
  enum {
    IsInteger = std::numeric_limits<Scalar>::is_signed,
    IsSigned = std::numeric_limits<Scalar>::is_signed,
    IsComplex = NumTraits< _Scalar >::IsComplex,
    RequireInitialization = NumTraits<Scalar>::RequireInitialization,
    ReadCost = 2 * NumTraits<Scalar>::ReadCost,
    AddCost = 2 * NumTraits<Scalar>::AddCost,
    MulCost = 4 * NumTraits<Scalar>::MulCost + 2 * NumTraits<Scalar>::AddCost
  };

  static inline Scalar epsilon() { return NumTraits<Scalar>::epsilon(); }
  static inline Scalar dummy_precision() { return NumTraits<Scalar>::dummy_precision(); }
};

namespace internal {

template<typename T> struct scalar_product_traits<T, DualNum<T> > {
  enum {
    // Cost = 2*NumTraits<T>::MulCost,
    Defined = 1
  };
  typedef DualNum<T> ReturnType;
};

template<typename T> struct scalar_product_traits<DualNum<T>, T> {
  enum {
    // Cost = 2*NumTraits<T>::MulCost,
    Defined = 1
  };
  typedef DualNum<T> ReturnType;
};

} // namespace internal
} // namespace Eigen

#endif
