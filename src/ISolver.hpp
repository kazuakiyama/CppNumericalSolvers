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
//#include <functional>
#include "Meta.hpp"

namespace pwie
{

typedef struct Options
{
  double gradTol;
  double tol;
  double rate;
  size_t maxIter;
  size_t m; // parameter for Lbfgs
  int verbosity;

  Options()
  {
    gradTol = 1e-9; // when the gradient is below this, goal reached
    tol = 1e-4;
    rate = 0.00005;
    maxIter = 100000;
    m = 10;
    verbosity = 0;
  }
} Options;

#include "lbfgs.h"

template <typename Func>
class ISolver : public Functor<Func>
{
private:

  typedef typename Func::InputType InputType;
  typedef typename Func::JacobianType JacobianType;

  lbfgs_parameter_t _param;
  std::list<InputType> _xHistory;

  inline InputType getLowerBound(int DIM, std::true_type) const {
    return Func::getLowerBound(DIM);
  }
  inline InputType getUpperBound(int DIM, std::true_type) const {
    return Func::getUpperBound(DIM);
  }
  inline InputType getLowerBound(int DIM, std::false_type) const {
    return  -INF * InputType::Ones(DIM);
  }
  inline InputType getUpperBound(int DIM, std::false_type) const {
    return  INF * InputType::Ones(DIM);
  }
  CREATE_MEMBER_FUNC_SIG_CHECK(getLowerBound, InputType (T::*)(int DIM) const, getLowerBound);
  CREATE_MEMBER_FUNC_SIG_CHECK(getUpperBound, InputType (T::*)(int DIM) const, getUpperBound);

public:

  ISolver() :
    Functor<Func>(),
    settings(),
    _linesearch(line_search_morethuente)
    //    _linesearch(line_search_backtracking),
    //    _linesearch(line_search_backtracking_owlqn),
  {
    lbfgs_parameter_init(&_param);
  }

  template <typename T0>
  ISolver(const T0 & f) :
    Functor<Func>(f),
    settings(),
    _linesearch(line_search_morethuente)
    //    _linesearch(line_search_backtracking),
    //    _linesearch(line_search_backtracking_owlqn),
  {
    lbfgs_parameter_init(&_param);
  }
  virtual ~ISolver() {}
  
  void solve(InputType & x0);
  std::list<InputType> & getXHistory() { return _xHistory; }

  virtual InputType getLowerBound(int DIM=Func::InputDim) const {
    return ISolver<Func>::getLowerBound(DIM, has_member_func_getLowerBound<Func>());
  }
  virtual InputType getUpperBound(int DIM=Func::InputDim) const {
    return ISolver<Func>::getUpperBound(DIM, has_member_func_getUpperBound<Func>());
  }

  Options settings;

protected:

  virtual void internalSolve(InputType & x0) = 0;
  double linesearch(const InputType & x, const JacobianType & direction);
  double linesearch(const InputType & x, const JacobianType & direction,
                    const Eigen::MatrixXd & hessian);
  double linesearch2(const InputType & x, const JacobianType & direction);

  int LineSearch(InputType & x,
                 const InputType & dx,
                 double & f,
                 JacobianType & g,
                 double & step);
private:
  line_search_proc _linesearch;
  static lbfgsfloatval_t lbfgs_evaluate(void *instance,
                                        const lbfgsfloatval_t *x,
                                        lbfgsfloatval_t *g,
                                        const int n,
                                        const lbfgsfloatval_t step);
};

} /* namespace pwie */

#include "CppNumericalSolvers/src/ISolver.cpp"

#endif /* ISOLVER_H_ */
