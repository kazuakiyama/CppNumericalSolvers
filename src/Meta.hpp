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
class Functor {

public:

  enum {
    InputDim = Func::InputsAtCompileTime,
    ValueDim = Func::ValuesAtCompileTime
  };
  typedef typename Func::Scalar Scalar;
#ifdef USE_COMPLEX_DUAL
  typedef std::complex<Scalar> DualScalar; // todo: use actual dual numbers
#else
  typedef Eigen::DualNum<Scalar> DualScalar;
#endif
  typedef typename Func::InputType InputType;
  typedef Eigen::Matrix<DualScalar,InputDim,1> DualInputType;
  typedef typename Func::JacobianType JacobianType;
  typedef Eigen::Matrix<Scalar,InputDim,InputDim> HessianType; // tensor by ValueDim

  //
  // Constructors / Destructors
  //
private:
  const Func & _func;

public:

  Functor(const Func & func)
    : _func(func)
  {
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(InputType);
    EIGEN_STATIC_ASSERT_VECTOR_ONLY(JacobianType);
    EIGEN_STATIC_ASSERT_SAME_VECTOR_SIZE(InputType,JacobianType);
  }
  virtual ~Functor() {}

  // evaluate function f(x)->v
  inline
  Scalar f(const InputType & x) const {
    return _func.f(x);
  }

  //ValueType operator() (const InputType & x) const { return Func::template f <Scalar> (x); }

  // template magic to create gradient function if none is provided
private:
  CREATE_MEMBER_FUNC_SIG_CHECK(gradient, void (T::*)(const InputType & x, JacobianType & grad) const, gradient);
  inline void gradient(const InputType & x, JacobianType & grad, std::true_type) const { _func.gradient(x, grad); }
  inline void gradient(const InputType & x, JacobianType & grad, std::false_type) const { gradientDual(x, grad); }
public:
  // return the gradient/jacobian
  inline void gradient(const InputType & x, JacobianType & grad) const {
    gradient(x, grad, has_member_func_gradient<Func>());
  }

  // template magic to create hessian function if none is provided
private:
  CREATE_MEMBER_FUNC_SIG_CHECK(hessian, void (T::*)(const InputType & x, HessianType & hes) const, hessian);
  inline void hessian(const InputType & x, HessianType & hes, std::true_type) const { _func.hessian(x, hes); }
  inline void hessian(const InputType & x, HessianType & hes, std::false_type) const { hessianFiniteDiff(x, hes); }
public:
  // evaluate hessian
  inline void hessian(const InputType & x, HessianType & hes) const {
    hessian(x, hes, has_member_func_hessian<Func>());
  }

private:
  CREATE_MEMBER_FUNC_SIG_CHECK(getLowerBound, InputType (T::*)(int DIM) const, getLowerBound);
  inline InputType getLowerBound(int DIM, std::true_type) const {
    return _func.getLowerBound(DIM);
  }
  inline InputType getLowerBound(int DIM, std::false_type) const {
    return -std::numeric_limits<Scalar>::max() * InputType::Ones(DIM);
  }
public:
  virtual InputType getLowerBound(int DIM=Func::InputDim) const {
    return getLowerBound(DIM, has_member_func_getLowerBound<Func>());
  }

private:
  CREATE_MEMBER_FUNC_SIG_CHECK(getUpperBound, InputType (T::*)(int DIM) const, getUpperBound);
  inline InputType getUpperBound(int DIM, std::true_type) const {
    return _func.getUpperBound(DIM);
  }
  inline InputType getUpperBound(int DIM, std::false_type) const {
    return  std::numeric_limits<Scalar>::max() * InputType::Ones(DIM);
  }
public:
  virtual InputType getUpperBound(int DIM=Func::InputDim) const {
    return getUpperBound(DIM, has_member_func_getUpperBound<Func>());
  }

private:
  CREATE_MEMBER_FUNC_SIG_CHECK(constraintDim, int (T::*)() const, constraintDim);
  inline int constraintDim(std::true_type) const {
    return _func.constraintDim();
  }
  inline int constraintDim(std::false_type) const {
    return 0;
  }
public:
  virtual inline int constraintDim() const {
    return constraintDim(has_member_func_constraintDim<Func>());
  }
private:
  CREATE_MEMBER_FUNC_SIG_CHECK(inputDim, int (T::*)() const, inputDim);
  inline int inputDim(std::true_type) const {
    return _func.inputDim();
  }
  inline int inputDim(std::false_type) const {
    return 0;
  }
public:
  virtual inline int inputDim() const {
    return inputDim(has_member_func_inputDim<Func>());
  }

  // double-check the gradient computation (using either finite-differences or dual numbers)
  Scalar checkGradient(const InputType & x,
                       JacobianType * Pgrad = NULL,
                       JacobianType * PgradD = NULL) const
  {
    const size_t DIM = x.rows();
    JacobianType grad(DIM), gradD(DIM);
    _func.gradient(x, grad);
    gradientDual(x, gradD);
    JacobianType diff = (gradD - grad);
    //Scalar error = diff.norm();
    size_t row, col;
    Scalar error = diff.cwiseAbs().maxCoeff(&row, &col);
    if (error > sqrt(std::numeric_limits<Scalar>::epsilon())) {
      std::cerr << "(checkGradient error)=" << error
                << " ||error||=" << diff.norm()
                << " @(" << row << "," << col << ")\n"
                << " *** g=" << grad(row,col) << " gD=" << gradD(row,col)
                << "\n";
    }
    if (Pgrad)
      *Pgrad = grad;
    if (PgradD)
      *PgradD = gradD;
    return error;
  }

  // calculate the gradient/jacobian using dual numbers
  void gradientDual(const InputType & x, JacobianType & grad) const {
#ifdef USE_COMPLEX_DUAL
    DualScalar eps(0, sqrt(std::numeric_limits<Scalar>::epsilon()));
    //std::cerr << "gDU, eps=" << eps << "\n";
#else
    DualScalar eps(0, 1.0);
    //std::cerr << "gDU, eps=" << eps << "\n";
#endif
    const size_t DIM = x.rows();
    JacobianType gg(DIM);
#pragma omp parallel
    {
      DualInputType xx = x.template cast<DualScalar>();
#pragma omp for
      for (size_t i = 0; i < DIM; i++) {
        xx[i] += eps;
#ifdef USE_COMPLEX_DUAL
        gg[i] = imag(_func.f(xx)) / imag(eps);
#else
        gg[i] = _func.f(xx).epart();
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
        finite[i] = (_func.f(xx) - _func.f(xy)) / (2.0 * eps);
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
        Scalar f4 = _func.f(xx);
        xx[i] += eps;
        xx[j] += eps;
        Scalar f1 = _func.f(xx);
        xx[j] -= eps;
        Scalar f2 = _func.f(xx);
        xx[j] += eps;
        xx[i] -= eps;
        Scalar f3 = _func.f(xx);
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
