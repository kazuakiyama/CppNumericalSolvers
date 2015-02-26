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

using namespace duals;
template <class eT> using DualNum = duals::dual<eT>;

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
