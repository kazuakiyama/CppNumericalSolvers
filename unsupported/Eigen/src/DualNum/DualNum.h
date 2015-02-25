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

#include <ostream>

namespace Eigen {
/* dual-numbers implementation for calculation of differentials */
template <class eT>
class DualNum {
private:
  eT _re, _du;
public:
  inline DualNum() : _re(), _du() {}
  inline DualNum(const eT & re, const eT & du = eT()) : _re(re), _du(du) {}

  inline eT real() const { return _re; }
  inline eT dual() const { return _du; }

  template <class eeT> friend std::ostream & operator<<(std::ostream & output, const DualNum<eeT> & rhs);
  template <class eeT> friend eeT real(const DualNum<eeT> & d) { return d.real(); }
  template <class eeT> friend eeT dual(const DualNum<eeT> & d) { return d.dual(); }

  // basic ops
  DualNum<eT> operator+() const { return *this; };
  DualNum<eT> operator+(const DualNum<eT> rhs) const { return DualNum<eT>(_re+rhs._re, _du+rhs._du); }
  template <class eeT>
  friend DualNum<eeT> operator+(const eeT lhs, const DualNum<eeT> rhs);
  DualNum<eT> operator-() const;
  DualNum<eT> operator-(const DualNum<eT> rhs) const;
  template <class eeT>
  friend DualNum<eeT> operator-(const eeT lhs, const DualNum<eeT> rhs);
  DualNum<eT> operator*(const DualNum<eT> rhs) const;
  template <class eeT>
  friend DualNum<eeT> operator*(const eeT lhs, const DualNum<eeT> rhs);
  template <class eeT>
  friend DualNum<eeT> operator/(const DualNum<eeT> lhs, const DualNum<eeT> rhs);
  template <class eeT>
  friend DualNum<eeT> operator/(const eeT lhs, const DualNum<eeT> rhs);
  template <class eeT>
  friend DualNum<eeT> operator/(const DualNum<eeT> lhs, const eeT rhs);
  DualNum<eT> & operator+=(DualNum<eT> rhs) {_re += rhs._re; _du += rhs._du; return *this; }
  DualNum<eT> & operator-=(DualNum<eT> rhs);
  DualNum<eT> & operator*=(DualNum<eT> rhs);
  DualNum<eT> & operator*=(eT rhs);
  DualNum<eT> & operator/=(eT rhs);

  // math.h functions
  template <class eeT> friend DualNum<eeT> pow(DualNum<eeT> x, eeT a);
  template <class eeT> friend DualNum<eeT> pow(DualNum<eeT> x, DualNum<eeT> a);
  template <class eeT> friend DualNum<eeT> exp(DualNum<eeT> x);
  template <class eeT> friend DualNum<eeT> log(DualNum<eeT> x);
  template <class eeT> friend DualNum<eeT> sin(DualNum<eeT> a);
  template <class eeT> friend DualNum<eeT> cos(DualNum<eeT> a);
  template <class eeT> friend DualNum<eeT> tan(DualNum<eeT> a);
  template <class eeT> friend DualNum<eeT> asin(DualNum<eeT> x);
  template <class eeT> friend DualNum<eeT> acos(DualNum<eeT> x);
  template <class eeT> friend DualNum<eeT> atan(DualNum<eeT> x);
  template <class eeT> friend DualNum<eeT> atan2(DualNum<eeT> y, DualNum<eeT> x);
  template <class eeT> friend DualNum<eeT> sqrt(DualNum<eeT> x);
  template <class eeT> friend DualNum<eeT> max(DualNum<eeT> x1, DualNum<eeT> x2);
  template <class eeT> friend DualNum<eeT> max(DualNum<eeT> x1, eeT x2);
  template <class eeT> friend DualNum<eeT> max(eeT x1, DualNum<eeT> x2);
  template <class eeT> friend DualNum<eeT> min(DualNum<eeT> x1, DualNum<eeT> x2);
  template <class eeT> friend DualNum<eeT> min(DualNum<eeT> x1, eeT x2);
  template <class eeT> friend DualNum<eeT> min(eeT x1, DualNum<eeT> x2);
  template <class eeT> friend DualNum<eeT> conj(const DualNum<eeT> & x) { return x; }
  //template <class eeT> friend DualNum<eeT> real(const DualNum<eeT> & x) { return x.real(); }
  template <class eeT> friend DualNum<eeT> imag(const DualNum<eeT> & x) { return DualNum<eeT>(0.0,0.0); }
  template <class eeT> friend DualNum<eeT> abs(const DualNum<eeT> & x) { return fabs(x._re); }
  template <class eeT> friend DualNum<eeT> abs2(const DualNum<eeT> & x) { return x * x; }

  // comparison
  template <class eeT> friend bool operator>(DualNum<eeT> lhs, DualNum<eeT> rhs);
  template <class eeT> friend bool operator>(eeT lhs, DualNum<eeT> rhs);
  template <class eeT> friend bool operator>(DualNum<eeT> lhs, eeT rhs);
  template <class eeT> friend bool operator>=(DualNum<eeT> lhs, DualNum<eeT> rhs);
  template <class eeT> friend bool operator>=(eeT lhs, DualNum<eeT> rhs);
  template <class eeT> friend bool operator>=(DualNum<eeT> lhs, eeT rhs);
  template <class eeT> friend bool operator<(DualNum<eeT> lhs, DualNum<eeT> rhs);
  template <class eeT> friend bool operator<(eeT lhs, DualNum<eeT> rhs);
  template <class eeT> friend bool operator<(DualNum<eeT> lhs, eeT rhs);
  template <class eeT> friend bool operator<=(DualNum<eeT> lhs, DualNum<eeT> rhs);
  template <class eeT> friend bool operator<=(eeT lhs, DualNum<eeT> rhs);
  template <class eeT> friend bool operator<=(DualNum<eeT> lhs, eeT rhs);
  template <class eeT> friend bool operator==(DualNum<eeT> lhs, DualNum<eeT> rhs);
  template <class eeT> friend bool operator==(eeT lhs, DualNum<eeT> rhs);
  template <class eeT> friend bool operator==(DualNum<eeT> lhs, eeT rhs);
  template <class eeT> friend bool operator!=(DualNum<eeT> lhs, DualNum<eeT> rhs);
  template <class eeT> friend bool operator!=(eeT lhs, DualNum<eeT> rhs);
  template <class eeT> friend bool operator!=(DualNum<eeT> lhs, eeT rhs);
};

}

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
  }
}

template <class eT>
std::ostream & operator<<(std::ostream & output, const DualNum<eT> & rhs)
{
  output << "(" << rhs._re << " + e*" << rhs._du << ")";
  return output;
}

#endif
