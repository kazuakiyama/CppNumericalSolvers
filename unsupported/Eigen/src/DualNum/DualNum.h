// -*-c++-*-
// This file is part of Eigen, a lightweight C++ template library
// for linear algebra. 
//
// Copyright (C) 2015 Michael Tesch tesch a tum de
// 
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
#ifndef EIGEN_DUALNUM_H
#define EIGEN_DUALNUM_H

#include <cmath>
#include "dual"

namespace Eigen {

using namespace cxxduals;
template <class eT> using DualNum = cxxduals::dual<eT>;

} // namespace Eigen

#include <Eigen/Core>
namespace Eigen {

template<typename _Scalar>
struct NumTraits<DualNum<_Scalar> > : GenericNumTraits<_Scalar>
{
  typedef DualNum<_Scalar> Real;
  typedef DualNum<typename NumTraits<_Scalar>::NonInteger> NonInteger;
  typedef DualNum<_Scalar> Nested;

  enum {
    IsInteger           =   NumTraits<_Scalar>::IsInteger,
    IsSigned            =   NumTraits<_Scalar>::IsSigned,
    IsComplex           =   NumTraits<_Scalar>::IsComplex,
    RequireInitialization = NumTraits<_Scalar>::RequireInitialization,
    ReadCost            = 2 * NumTraits<_Scalar>::ReadCost,
    AddCost             = 2 * NumTraits<_Scalar>::AddCost,
    MulCost             = 4 * NumTraits<_Scalar>::MulCost + 2 * NumTraits<_Scalar>::AddCost
  };

  static inline Real epsilon()          { return Real(NumTraits<_Scalar>::epsilon()); }
  static inline Real dummy_precision()  { return Real(NumTraits<_Scalar>::dummy_precision()); }
};

namespace internal {

template<typename T> struct scalar_product_traits<T, DualNum<T> > {
  enum {
    Cost = 2*NumTraits<T>::MulCost,
    Defined = 1
  };
  typedef DualNum<T> ReturnType;
};

template<typename T> struct scalar_product_traits<DualNum<T>, T> {
  enum {
    Cost = 2*NumTraits<T>::MulCost,
    Defined = 1
  };
  typedef DualNum<T> ReturnType;
};

template <typename _Tp>
struct conj_helper<DualNum<_Tp>, _Tp, false, false>
{
  EIGEN_STRONG_INLINE DualNum<_Tp> pmadd(const DualNum<_Tp> & x, const DualNum<_Tp> & y,
                                         const DualNum<_Tp>& c) const
  { return x * y + c; }

  EIGEN_STRONG_INLINE DualNum<_Tp> pmul(const DualNum<_Tp> & a, const DualNum<_Tp> & b) const
  { return a * b; }
};

// Eigen needs this to print
template <typename _Tp>
struct cast_impl<DualNum<_Tp>, int>
{
  static inline int run(const DualNum<_Tp> & x) {
    return int(x.part(0));
  }
};

} // namespace internal
} // namespace Eigen

#endif
