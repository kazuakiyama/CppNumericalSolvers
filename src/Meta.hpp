/**
 * Copyright (c) 2014 Patrick Wieschollek
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */


#ifndef META_H_
#define META_H_

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <stdexcept>
#include <unsupported/Eigen/DualNum>

#define USE_COMPLEX_DUAL

namespace pwie
{

template <typename Func>
double checkGradient(const Func & func, const typename Func::InputType & x, 
                     const typename Func::JacobianType & grad, const double eps = 1e-8);
template <typename Func>
double checkGradientC(const Func & func, const typename Func::InputType & x,
                      const typename Func::JacobianType & grad, const double eps = 1e-10);
template <typename Func>
void computeGradient(const Func & func, const typename Func::InputType & x,
                     typename Func::JacobianType & grad, const double eps = 1e-8);
template <typename Func>
void computeGradientC(const Func & func, const typename Func::InputType & x,
                      typename Func::JacobianType & grad, const double eps = 1e-10);
template <typename Func>
void computeHessian(const Func & func, const typename Func::InputType & x, 
                    typename Func::HessianType & hessian, const double eps = 1e-8);

// Check for member function with given name and signature.
template<typename T, T>
struct sig_check : std::true_type {};

#define CREATE_MEMBER_FUNC_SIG_CHECK(func_name, func_sig, templ_postfix)    \
template<typename T, typename = std::true_type>                             \
struct has_member_func_##templ_postfix : std::false_type {};                \
template<typename T>                                                        \
struct has_member_func_##templ_postfix<                                     \
    T, std::integral_constant<                                              \
        bool                                                                \
        , sig_check<func_sig, &T::func_name>::value                         \
    >                                                                       \
> : std::true_type {}

template<typename Func>
class Functor : public Func {
public:
  typedef typename Func::Scalar Scalar;
#ifdef USE_COMPLEX_DUAL
  typedef std::complex<Scalar> DualScalar; // todo: use actual dual numbers
#else
  typedef Eigen::DualNum<Scalar> DualScalar;
#endif
  enum {
    InputDim = Func::InputsAtCompileTime,
    ValueDim = Func::ValuesAtCompileTime
  };
  typedef typename Func::InputType InputType;
  typedef Eigen::Matrix<DualScalar,InputDim,1> DualInputType;
  typedef typename Func::JacobianType JacobianType;
  typedef Eigen::Matrix<Scalar,InputDim,InputDim> HessianType; // tensor by ValueDim

  //
  // Constructors / Destructors
  //
  Functor() : Func() {
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(InputType);
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(JacobianType);
    EIGEN_STATIC_ASSERT_SAME_VECTOR_SIZE(InputType,JacobianType);
  }
  template <typename T0> Functor(const T0 & f) : Func(f) {}
  virtual ~Functor() {}

  // evaluate function f(x)->v
  //ValueType operator() (const InputType & x) const { return Func::template f <Scalar> (x); }

  // template magic to create gradient function if none is provided
  CREATE_MEMBER_FUNC_SIG_CHECK(gradient, void (T::*)(const InputType & x, JacobianType & grad) const, gradient);
  inline void gradient(const InputType & x, JacobianType & grad, std::true_type) const { Func::gradient(x, grad); }
  inline void gradient(const InputType & x, JacobianType & grad, std::false_type) const { gradientDual(x, grad); }
  // return the gradient/jacobian
  inline void gradient(const InputType & x, JacobianType & grad) const {
    gradient(x, grad, has_member_func_gradient<Func>());
  }

  // template magic to create hessian function if none is provided
  CREATE_MEMBER_FUNC_SIG_CHECK(hessian, void (T::*)(const InputType & x, HessianType & hes) const, hessian);
  inline void hessian(const InputType & x, HessianType & hes, std::true_type) const { Func::hessian(x, hes); }
  inline void hessian(const InputType & x, HessianType & hes, std::false_type) const { hessianFiniteDiff(x, hes); }
  // evaluate hessian
  inline void hessian(const InputType & x, HessianType & hes) const {
    hessian(x, hes, has_member_func_hessian<Func>());
  }

  // double-check the gradient computation (using either finite-differences or dual numbers)
  Scalar checkGradient(const InputType & x) const {
    const size_t DIM = x.rows();
    JacobianType grad(DIM), gradD(DIM), gradFD(DIM);
    Func::gradient(x, grad);
    gradientDual(x, gradD);
    size_t row, col;
    JacobianType diff = (gradD - grad);
    //Scalar error = diff.norm();
    Scalar error = diff.cwiseAbs().maxCoeff(&row, &col);
    if (error > sqrt(std::numeric_limits<Scalar>::epsilon())) {
      gradientFiniteDiff(x, gradFD);
      std::cerr << "(checkGradient error)=" << error
                << "(" << row << "," << col << ") "
        "g=" << grad(row,col) << " gD=" << gradD(row,col) << " gFD=" << gradFD(row,col)
                << "\n";
    }
    return error;
  }

  // calculate the gradient/jacobian using dual numbers
  void gradientDual(const InputType & x, JacobianType & grad) const {
    Scalar eps = sqrt(std::numeric_limits<Scalar>::epsilon());
    std::cerr << "gDU, eps=" << eps << "\n";
    const size_t DIM = x.rows();
    JacobianType gg(DIM);
#pragma omp parallel
    {
      DualInputType xx = x.template cast<DualScalar>();
#pragma omp for
      for (size_t i = 0; i < DIM; i++) {
#ifdef USE_COMPLEX_DUAL
        xx[i] += DualScalar(0, eps);
        gg[i] = imag(Func::f(xx)) / eps;
#else
        xx[i] += DualScalar(0, 1.0);
        gg[i] = Func::f(xx).ipart();
#endif
        xx[i] = x[i];
      }
    }
    grad = gg;
  }

  // calculate the gradient using finite differences
  void gradientFiniteDiff(const InputType & x, JacobianType & grad,
                          Scalar eps = sqrt(std::numeric_limits<Scalar>::epsilon())) const
  {
    eps = sqrt(eps);
    std::cerr << "gFD, eps=" << eps << "\n";
    const size_t DIM = x.rows();
    typename Func::JacobianType finite(DIM);
#pragma omp parallel
    {
      typename Func::InputType xx = x;
      typename Func::InputType xy = x;
#pragma omp for
      for(size_t i = 0; i < DIM; i++) {
        xx[i] += eps;
        xy[i] -= eps;
        finite[i] = (Func::f(xx) - Func::f(xy)) / (2.0 * eps);
        xx[i] = x[i];
        xy[i] = x[i];
      }
    }
    grad = finite;
  }

  // calculate the hessian using finite differences
  void hessianFiniteDiff(const InputType & x, HessianType & hes, Scalar eps = 1e-6) const {
    const size_t DIM = x.rows();
    typename Func::InputType xx = x;
    hes.resize(DIM,DIM);
    for (size_t i = 0; i < DIM; i++) {
      for (size_t j = 0; j < DIM; j++) {
        Scalar f4 = Func::f(xx);
        xx[i] += eps;
        xx[j] += eps;
        Scalar f1 = Func::f(xx);
        xx[j] -= eps;
        Scalar f2 = Func::f(xx);
        xx[j] += eps;
        xx[i] -= eps;
        Scalar f3 = Func::f(xx);
        hes(i, j) = (f1 - f2 - f3 + f4) / (eps * eps);
        xx[i] = x[i];
        xx[j] = x[j];
      }
    }
  }
  // void hessianHyperDual(const InputType & x, HessianType & hes) const {
  // }
};

template<typename T>
bool AssertSimiliar(T a, T b)
{
  return std::abs(a - b) <=  1e-2;
}
template<typename T>
bool AssertGreaterThan(T a, T b)
{
  return (a - b) > ((std::abs(a) < std::abs(b) ? std::abs(b) : std::abs(a)) * 1e-3);
}
template<typename T>
bool AssertLessThan(T a, T b)
{
  return (b - a) > ((std::abs(a) < std::abs(b) ? std::abs(b) : std::abs(a)) * 1e-3);
}
template<typename T>
bool AssertEqual(T a, T b)
{
  return (a == b);
}
}

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define Assert(x,m) if (!(x)) { throw (std::runtime_error(m)); }

#define UNUSED(arg) (void)arg;    //  trick of Qt

#define FAST

#ifdef FAST
#define Debug(x)
#else
#define Debug(x) if(true){std::cout << "DEBUG: "<< x;std::cout<< std::endl;}
#endif

#include "CppNumericalSolvers/src/Meta.cpp"

#endif /* META_H_ */
