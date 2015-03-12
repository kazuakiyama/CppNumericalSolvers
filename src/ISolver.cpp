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
  const Eigen::Map<const Eigen::VectorXd> xPROXY(x, n);
  const InputType xx = xPROXY.cast<Scalar>();;
  Eigen::Map<Eigen::VectorXd> gPROXY(g, n);
  JacobianType gg = gPROXY.cast<Scalar>();
  ISolver<Func> * that = (ISolver<Func> *)instance;
  that->_functor.gradient(xx, gg);
  gPROXY = gg.template cast<double>();
  return (lbfgsfloatval_t)that->_functor.f(xx);
}

template <typename Func>
int
ISolver<Func>::LineSearch(InputType & x,
                          const InputType & dx,
                          Scalar & f,
                          JacobianType & g,
                          Scalar & step)
{
  Eigen::VectorXd xx = x.template cast<double>();
  Eigen::VectorXd xp = x.template cast<double>();
  Eigen::VectorXd dd = dx.template cast<double>();
  f = _functor.f(x);
  double ff = (double)f;
  double stepPROXY = (double)step;
  _functor.gradient(x, g);
  Eigen::VectorXd gg = g.template cast<double>();
  Eigen::VectorXd gp = g.template cast<double>();
  Eigen::VectorXd wp(x.rows());

  /* Construct a callback data. */
  callback_data_t cd;
  cd.n = x.rows();
  cd.instance = (void *)this;
  cd.proc_evaluate = &ISolver<Func>::lbfgs_evaluate;
  cd.proc_progress = NULL;

  int ret;
  ret = _linesearch(cd.n,
                     xx.data(),  /* current position */
                     &ff,        /* current function value @x */
                     gg.data(),   /* current function gradient @x */
                     dd.data(),  /* search direction (const) */
                     &stepPROXY, /* step size - initial & final */
                     xp.data(),
                     gp.data(),  /* (const) */
                     wp.data(),  /* (const) */
                     &cd,
                     &_param);   /* (const)  */
  x = xx.cast<Scalar>();
  f = ff;
  g = gg.cast<Scalar>();
  step = stepPROXY;
  return ret;
}

template <typename Func>
typename ISolver<Func>::Scalar
ISolver<Func>::linesearch(const InputType & xp, const JacobianType & direction)
{
  Scalar step = 1.0;
  Scalar f;
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
typename ISolver<Func>::Scalar
ISolver<Func>::linesearch2(const InputType & x, const JacobianType & direction)
{
  const Scalar alpha = 0.001; // c1
  const Scalar beta = 0.1;    // c2
  Scalar t = 1.0;

  InputType xx = x + t * direction;
  Scalar f = _functor.f(xx);
  const Scalar f_in = f;
  JacobianType grad(x.rows());
  _functor.gradient(x, grad);
  const Scalar Cache = alpha * grad.dot(direction);

  while (f > f_in + t * Cache) {
    t *= beta;
    xx = x + t * direction;
    f = _functor.f(xx);
  }

  return t;
}

template <typename Func>
typename ISolver<Func>::Scalar
ISolver<Func>::linesearch(const InputType & x, const JacobianType & direction,
                          const Eigen::MatrixXd & hessian)
{
    const Scalar alpha = 0.2;
    const Scalar beta = 0.9;
    Scalar t = 1.0;
    InputType xx = x + t * direction;
    Scalar f = _functor.f(xx);
    const Scalar f_in = f;
    JacobianType grad(x.rows());
    _functor.gradient(x, grad);
    const Scalar Cache = alpha * grad.dot(direction) + 0.5*alpha*direction.transpose()*(hessian*direction);

    while(f > f_in + t * Cache)
    {
        t *= beta;
        f = FunctionValue(x + t * direction);
    }

    return t;

}

} /* namespace pwie */
