/**
 * Copyright (c) 2014 Patrick Wieschollek
 * Copyright (c) 2015 Michael Tesch
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

#include "AsaCgSolver.hpp"
#include <iostream>
#include <list>
#include "stopwatch.hpp"

extern "C" {
#include "asa_user.h"
}

namespace pwie
{

template <typename Func>
AsaCgSolver<Func>::AsaCgSolver(const Func & func)
  : ISolver<Func>(func)
{
}

template <typename T>
class Global {
  typedef T * tptr;
public:
  static tptr & get() {
    static T * thing;
    return thing;
  }
private:
  Global() {}
  Global(const Global& rhs) { (void)rhs; }
  void operator=(const Global& rhs) { (void)rhs; }
};

template <typename FUNCTOR>
static double myvalue(asa_objective *asa)
{
  Eigen::Map<Eigen::VectorXd> x(asa->x, asa->n);
  typename FUNCTOR::InputType xx = x.cast<typename FUNCTOR::InputType::Scalar>();
  return (double)Global<FUNCTOR>::get()->f(xx);
}

template <typename FUNCTOR>
void mygrad(asa_objective *asa)
{
  Eigen::Map<Eigen::VectorXd> x(asa->x, asa->n);
  Eigen::Map<Eigen::VectorXd> g(asa->g, asa->n);
  typename FUNCTOR::InputType xx = x.cast<typename FUNCTOR::InputType::Scalar>();
  typename FUNCTOR::JacobianType grad(asa->n);
  Global<FUNCTOR>::get()->gradient(xx, grad);
  g = grad.template cast<double>();
}

template <typename FUNCTOR>
double myvalgrad(asa_objective *asa)
{
  Eigen::Map<Eigen::VectorXd> x(asa->x, asa->n);
  Eigen::Map<Eigen::VectorXd> g(asa->g, asa->n);
  typename FUNCTOR::InputType xx = x.cast<typename FUNCTOR::InputType::Scalar>();
  typename FUNCTOR::JacobianType grad(asa->n);
  Global<FUNCTOR>::get()->gradient(xx, grad);
  g = grad.template cast<double>();
  return (double)Global<FUNCTOR>::get()->f(xx);
}

template <typename Func>
void
AsaCgSolver<Func>::internalSolve(InputType & x0)
{
    int DIM = x0.rows();

    _lb = _functor.getLowerBound();
    _ub = _functor.getUpperBound();

    /* if you want to change parameter value, you need the following: */
    asacg_parm cgParm;
    asa_parm asaParm;
    asa_stat asaStat;
    
    /* allocate arrays for problem solution and bounds */
    Eigen::VectorXd x  = x0.template cast<double>();
    Eigen::VectorXd lo = _lb.template cast<double>();
    Eigen::VectorXd hi = _ub.template cast<double>();

    /* if you want to change parameter value, initialize strucs with default */
    asa_cg_default(&cgParm) ;
    asa_default(&asaParm) ;

    if (settings.verbosity) {
      std::cout << "internalSolve ||x||=" << x.norm() << " ||hi||="
                << hi.norm() << " ||lo||=" << lo.norm() << "\n";
    }

    /* if you want to change parameters, change them here: */
    cgParm.PrintLevel = settings.verbosity > 1 ? settings.verbosity - 1 : 0;
    cgParm.PrintParms = settings.verbosity > 0;
    cgParm.maxit = settings.maxIter;
    asaParm.PrintLevel = settings.verbosity > 1 ? settings.verbosity - 1 : 0;
    asaParm.PrintParms = settings.verbosity > 0;
    asaParm.PrintFinal = settings.verbosity > 0;
    asaParm.HardConstraint = TRUE;

    /* run the code */
    Global<FunctorGlobalType>::get() = &_functor;
    asa_cg(x.data(), lo.data(), hi.data(),
           DIM, &asaStat, &cgParm, &asaParm,
           (double)settings.gradTol, /* grad_tol */
           &myvalue<FunctorGlobalType>,
           &mygrad<FunctorGlobalType>,
           &myvalgrad<FunctorGlobalType>,
           NULL, NULL);

    x0 = x.cast<Scalar>();
    settings.numIters = asaStat.cgiter;
}
}

/* namespace pwie */
