/**
 * Copyright (c) 2014-2015 Patrick Wieschollek
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

#include "ISolver.hpp"

namespace pwie
{

template <typename Func>
ISolver<Func>::ISolver()
  : settings(),
    _linesearch(line_search_morethuente)
    //    _linesearch(line_search_backtracking),
    //    _linesearch(line_search_backtracking_owlqn),
{
  lbfgs_parameter_init(&_param);
}

template <typename Func>
ISolver<Func>::~ISolver()
{
  // TODO Auto-generated destructor stub
}

template <typename Func>
void
ISolver<Func>::solve(InputType & x0)
{
  internalSolve(x0);
}

template <typename Func>
lbfgsfloatval_t ISolver<Func>::lbfgs_evaluate(void *instance,
                                              const lbfgsfloatval_t *x,
                                              lbfgsfloatval_t *g,
                                              const int n,
                                              const lbfgsfloatval_t step
                                              )
{
  (void)step;
  const Eigen::Map<const typename Func::InputType> xPROXY(x, n);
  const InputType xx = xPROXY;
  Eigen::Map<typename Func::InputType> gPROXY(g, n);
  JacobianType gg = gPROXY;
  ISolver<Func> * that = (ISolver<Func> *)instance;
  that->gradient(xx, gg);
  gPROXY = gg;
  return that->f(xx);
}

template <typename Func>
int
ISolver<Func>::LineSearch(InputType & x,
                          const InputType & dx,
                          double & f,
                          JacobianType & g,
                          double & step)
{
  InputType xp = x;
  f = Func::f(x);
  this->gradient(x, g);
  JacobianType gp = g;
  InputType wp(x.rows());

  /* Construct a callback data. */
  callback_data_t cd;
  cd.n = x.rows();
  cd.instance = (void *)this;
  cd.proc_evaluate = &ISolver<Func>::lbfgs_evaluate;
  cd.proc_progress = NULL;

  return _linesearch(cd.n,
                     x.data(),   /* current position */
                     &f,         /* current function value @x */
                     g.data(),   /* current function gradient @x */
                     dx.data(),  /* search direction */
                     &step,      /* step size - initial & final */
                     xp.data(),
                     gp.data(),
                     wp.data(),
                     &cd,
                     &_param);
}

template <typename Func>
double ISolver<Func>::linesearch(const InputType & xp, const JacobianType & direction)
{
  lbfgsfloatval_t step = 1.0;
  lbfgsfloatval_t f;
  InputType x = xp;
  JacobianType g(xp.rows());
  int status = LineSearch(x, direction, f, g, step);
  if (status < 0) {
    std::cerr << "_linesearch failed: " << status << "\n";
    throw std::exception();
  }
  return step;
}

template <typename Func>
double ISolver<Func>::linesearch2(const InputType & x, const JacobianType & direction)
{
  const double alpha = 0.001; // c1
  const double beta = 0.1;    // c2
  double t = 1.0;

  InputType xx = x + t * direction;
  double f = Func::f(xx);
  const double f_in = Func::f(x);
  JacobianType grad(x.rows());
  this->gradient(x, grad);
  const double Cache = alpha * grad.dot(direction);

  while (f > f_in + t * Cache) {
    t *= beta;
    xx = x + t * direction;
    f = Func::f(xx);
  }

  return t;
}

double ISolver::linesearch(const Vector & x, const Vector & direction,
                         const Eigen::MatrixXd & hessian,
                           const FunctionOracleType & FunctionValue,
                           const GradientOracleType & FunctionGradient)
{

    const double alpha = 0.2;
    const double beta = 0.9;
    double t = 1.0;

    double f = FunctionValue(x + t * direction);
    const double f_in = FunctionValue(x);
    Vector grad(x.rows());
    FunctionGradient(x, grad);
    const double Cache = alpha * grad.dot(direction) + 0.5*alpha*direction.transpose()*(hessian*direction);

    while(f > f_in + t * Cache)
    {
        t *= beta;
        f = FunctionValue(x + t * direction);
    }

    return t;

}

} /* namespace pwie */
