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

template <typename Func>
static ISolver<Func> * getSolver(const Func & f, solver_id type)
{
  switch (type) {
  case SOLVER_GRADD:    return new GradientDescentSolver<Func>(f);
  case SOLVER_NEWTON:   return new NewtonDescentSolver<Func>(f);
  case SOLVER_CG:       return new ConjugateGradientSolver<Func>(f);
  case SOLVER_ASA_CG:   return new AsaCgSolver<Func>(f);
  case SOLVER_BFGS:     return new BfgsSolver<Func>(f);
  case SOLVER_LBFGS:    return new LbfgsSolver<Func>(f);
  case SOLVER_LBFGSB:   return new LbfgsbSolver<Func>(f);
#ifdef HAVE_IPOPT
  case SOLVER_IPOPT:    return new IpoptSolver<Func>(f);
#endif
  }
  std::cerr << "ERROR: unknown solver type requested:" << type << "\n";
  throw std::exception();
}

}

#endif

