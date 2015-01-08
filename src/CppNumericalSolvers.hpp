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
#ifndef CPPNUMERICAL_SOLVERS_H_
#define CPPNUMERICAL_SOLVERS_H_

/*
 * This file includes all of the Numerical Solvers and provides a generator
 * function for them all.
 */
#include "ConjugateGradientSolver.hpp"
#include "NewtonDescentSolver.hpp"
#include "LbfgsSolver.hpp"
#include "AsaCgSolver.hpp"
#include "GradientDescentSolver.hpp"
#include "BfgsSolver.hpp"
#include "LbfgsSolver.hpp"
#include "LbfgsbSolver.hpp"
#ifdef HAVE_IPOPT
#include "IpoptSolver.hpp"
#endif

namespace pwie {

typedef enum {
  SOLVER_GRADD,
  SOLVER_NEWTON,
  SOLVER_CG,
  SOLVER_ASA_CG,
  SOLVER_BFGS,
  SOLVER_LBFGS,
  SOLVER_LBFGSB,
#ifndef HAVE_IPOPT
#define SOLVER_COUNT (SOLVER_LBFGSB+1)
#else
  SOLVER_IPOPT,
#define SOLVER_COUNT (SOLVER_IPOPT+1)
#endif
} solver_id;

template <typename FUNCTOR>
static pwie::ISolver<FUNCTOR> * getSolver(solver_id type)
{
  switch (type) {
  case SOLVER_GRADD:    return new pwie::GradientDescentSolver<FUNCTOR>();
  case SOLVER_NEWTON:   return new pwie::NewtonDescentSolver<FUNCTOR>();
  case SOLVER_CG:       return new pwie::ConjugateGradientSolver<FUNCTOR>();
  case SOLVER_ASA_CG:   return new pwie::AsaCgSolver<FUNCTOR>();
  case SOLVER_BFGS:     return new pwie::BfgsSolver<FUNCTOR>();
  case SOLVER_LBFGS:    return new pwie::LbfgsSolver<FUNCTOR>();
  case SOLVER_LBFGSB:   return new pwie::LbfgsbSolver<FUNCTOR>();
#ifdef HAVE_IPOPT
  case SOLVER_IPOPT:    return new pwie::IpoptSolver<FUNCTOR>();
#endif
  }
  std::cerr << "ERROR: unknown solver type requested:" << type << "\n";
  throw std::exception();
}

}

#endif

