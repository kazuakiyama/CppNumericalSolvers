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

#ifndef ISOLVER_H_
#define ISOLVER_H_

#include <list>
#include "Meta.hpp"

namespace pwie
{

template <typename _Scalar>
struct Options
{
  typedef _Scalar Scalar;
  Scalar gradTol;       // minimum gradient norm
  Scalar tol;           // 
  Scalar rate;          // 
  size_t maxIter;       // maximum solver iterations
  size_t numIters;      // actual number of solver iterations
  size_t m;             // parameter for Lbfgs[b] - size of hessian estimate
  int verbosity;        // verbosity level, 0=quiet

  Options()
  {
    // stop when the ||gradient||_inf is below this.
    gradTol = 1e-9;
    // stop when successive 
    tol = 1e-9;
    rate = 0.00005;
    maxIter = 100000;
    numIters = 0;
    m = 10;
    verbosity = 0;
  }
};

#include "lbfgs.h"

template <typename Func>
class ISolver
{
protected:

  typedef typename Func::Scalar Scalar;
  typedef typename Func::InputType InputType;
  typedef typename Func::JacobianType JacobianType;
  typedef typename Functor<Func>::HessianType HessianType;
  typedef Functor<Func> FunctorType;

  FunctorType _functor;

private:
  lbfgs_parameter_t _param;
  line_search_proc _linesearch;

public:

  Options<Scalar> settings;

  ISolver(const Func & func,
          line_search_proc linesearch = line_search_morethuente)
    : _functor(func),
      _linesearch(linesearch),
      settings()
  {
    lbfgs_parameter_init(&_param);
  }
  virtual ~ISolver() {}

  void solve(InputType & x0);

protected:

  virtual void internalSolve(InputType & x0) = 0;
  Scalar linesearch(const InputType & x, const JacobianType & direction);
  Scalar linesearch(const InputType & x, const JacobianType & direction,
                    const Eigen::MatrixXd & hessian);
  Scalar linesearch2(const InputType & x, const JacobianType & direction);

  virtual int LineSearch(InputType & x,
                         const InputType & dx,
                         Scalar & f,
                         JacobianType & g,
                         Scalar & step);

private:
  CREATE_MEMBER_FUNC_SIG_CHECK(checkConverged,
                               bool (T::*)(int iter, InputType & x, Scalar & step,
                                           InputType & dir, JacobianType & grad) const,
                               checkConverged);
  inline bool checkConverged(int iter, InputType & x, Scalar & step, InputType & dir,
                             JacobianType & g, std::true_type) const
  {
    return _functor.checkConverged(x, step, dir, g);
  }
  inline bool checkConverged(size_t iter, InputType & x, Scalar & step, InputType & dir,
                             JacobianType & g, std::false_type) const
  {
    (void)x;
    (void)step;
    (void)dir;
    if (iter >= settings.maxIter)
      return true;
    return g.template lpNorm<Eigen::Infinity>() < settings.tol;
  }

protected:
  virtual bool checkConverged(int iter, InputType & x, Scalar & step,
                              InputType & dir, JacobianType & grad) const
  {
    return checkConverged(iter, x, step, dir, grad,
                          has_member_func_checkConverged<Func>());
  }

private:

  static lbfgsfloatval_t lbfgs_evaluate(void *instance,
                                        const lbfgsfloatval_t *x,
                                        lbfgsfloatval_t *g,
                                        const int n,
                                        const lbfgsfloatval_t step);
};

} /* namespace pwie */

#include "CppNumericalSolvers/src/ISolver.cpp"

#endif /* ISOLVER_H_ */
