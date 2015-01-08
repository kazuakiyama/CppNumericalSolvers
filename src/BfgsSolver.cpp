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

#include "BfgsSolver.hpp"
#include <iostream>
namespace pwie
{

template <typename Func>
BfgsSolver<Func>::BfgsSolver() : ISolver<Func>()
{

}


template <typename Func>
void
BfgsSolver<Func>::BfgsSolver::internalSolve(InputType & x)
{
    const size_t DIM = x.rows();
    size_t iter = 0;
    HessianType H = HessianType::Identity(DIM, DIM);
    JacobianType grad(DIM);

    InputType x_old = x;

    do
    {
        this->gradient(x, grad);
        JacobianType p = -1 * H * grad;
        const double rate = this->linesearch(x, p);
        x = x + rate * p;
        JacobianType grad_old = grad;
        this->gradient(x, grad);

        InputType s = x - x_old;
        JacobianType y = grad - grad_old;

        const double rho = 1.0 / y.dot(s);
        if(iter == 0)
        {
            H = ((y.dot(s)) / (y.dot(y)) * HessianType::Identity(DIM, DIM));
        }
        H = H - rho * (s * (y.transpose() * H) + (H * y) * s.transpose())
          + rho * rho * (y.dot(H * y) + 1.0 / rho) * (s * s.transpose());
        x_old = x;
        iter++;

    }
    while((grad.template lpNorm<Eigen::Infinity>() > settings.gradTol) &&
          (iter < settings.maxIter));
}
}

/* namespace pwie */
