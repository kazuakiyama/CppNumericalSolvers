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
AsaCgSolver<Func>::AsaCgSolver() : ISolver<Func>()
{
    // TODO Auto-generated constructor stub
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
  Eigen::Map<typename FUNCTOR::InputType> x(asa->x, asa->n);
  return Global<FUNCTOR>::get()->f(typename FUNCTOR::InputType(x));
}

template <typename FUNCTOR>
void mygrad(asa_objective *asa)
{
  Eigen::Map<typename FUNCTOR::InputType> x(asa->x, asa->n);
  Eigen::Map<typename FUNCTOR::InputType> g(asa->g, asa->n);
  typename FUNCTOR::JacobianType grad(asa->n);
  Global<FUNCTOR>::get()->gradient(typename FUNCTOR::InputType(x), grad);
  g = grad;
}

template <typename FUNCTOR>
double myvalgrad(asa_objective *asa)
{
  Eigen::Map<typename FUNCTOR::InputType> x(asa->x, asa->n);
  Eigen::Map<typename FUNCTOR::JacobianType> g(asa->g, asa->n);
  typename FUNCTOR::JacobianType grad(asa->n);
  Global<FUNCTOR>::get()->gradient(typename FUNCTOR::InputType(x), grad);
  g = grad;
  return Global<FUNCTOR>::get()->f(typename FUNCTOR::InputType(x));
}

template <typename Func>
void
AsaCgSolver<Func>::internalSolve(InputType & x0)
{
    int DIM = x0.rows();

    _lower = this->getLowerBound(DIM);
    _upper = this->getUpperBound(DIM);

    /* if you want to change parameter value, you need the following: */
    asacg_parm cgParm;
    asa_parm asaParm;

    /* allocate arrays for problem solution and bounds */
    double * x  = x0.data();
    double * lo = _lower.data();
    double * hi = _upper.data();

    /* if you want to change parameter value, initialize strucs with default */
    asa_cg_default (&cgParm) ;
    asa_default (&asaParm) ;

    /* if you want to change parameters, change them here: */
    //cgParm.PrintParms = TRUE ;
    cgParm.PrintLevel = 0 ;
    asaParm.HardConstraint = TRUE ;
    //asaParm.PrintParms = TRUE ;
    asaParm.PrintLevel = 0 ;

    /* run the code */
    Global<Functor<Func> >::get() = this;
    asa_cg(x, lo, hi, DIM, NULL, &cgParm, &asaParm, 1.e-9,
           &myvalue<Functor<Func> >, &mygrad<Functor<Func> >, &myvalgrad<Functor<Func> >,
           NULL, NULL);

    // x0 = x;
}
}

/* namespace pwie */
