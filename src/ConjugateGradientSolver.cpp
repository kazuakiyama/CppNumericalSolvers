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

#include "ConjugateGradientSolver.hpp"
#include "stopwatch.hpp"
#include <iostream>
namespace pwie
{

template <typename Func>
ConjugateGradientSolver<Func>::ConjugateGradientSolver(const Func & func)
  : ISolver<Func>(func)
{
}

template <typename Func>
void
ConjugateGradientSolver<Func>::internalSolve(InputType & x)
{
  size_t iter = 0;

  JacobianType grad(x.rows());
  JacobianType grad_old(x.rows());
  JacobianType Si(x.rows());
  JacobianType Si_old(x.rows());
  Stopwatch<> stopwatch;

  _functor.gradient(x, grad);

  while ((grad.template lpNorm<Eigen::Infinity>() > settings.gradTol) &&
         (iter < settings.maxIter)) { 

    if (iter==0){
      Si = -grad;
    } else {
      Scalar beta;
      // fletcher-reeves
      beta = grad.dot(grad)/(grad_old.dot(grad_old));
      // polak ribiere
      //beta = grad.dot(grad - grad_old)/(grad_old.dot(grad_old));
      // hestenes-stiefel
      //beta = -grad.dot(grad - grad_old)/Si_old.dot(grad - grad_old);
      // dai-yuan
      //beta = -grad.dot(grad)/Si_old.dot(grad - grad_old);
      beta = MAX(beta, 0); // optional, direction reset
      Si = -grad + beta*Si_old;
    }

    const Scalar rate = this->linesearch(x, Si);

    x = x + rate * Si;

    iter++;
    grad_old = grad;
    Si_old = Si;

    if (0) { 
      stopwatch.stop();
      std::cout << "iteration " << iter << " " << _functor.f(x) 
                << " dt=" << stopwatch.elapsed()/1e3 << std::endl;
      stopwatch.start();
    }
    settings.numIters = iter;
    _functor.gradient(x, grad);
  }
}
}

/* namespace pwie */
