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

#include "Meta.hpp"
#include <cmath>
#include <iostream>
#include <complex>

namespace pwie
{

template <typename Func>
void computeGradient(const Func & func, const typename Func::InputType & x,
                     typename Func::JacobianType & grad, const double eps)
{
    const size_t DIM = x.rows();
    typename Func::JacobianType finite(DIM);
    for(size_t i = 0; i < DIM; i++)
    {
        typename Func::InputType xx = x;
        xx[i] += eps;
        typename Func::InputType xy = x;
        xy[i] -= eps;
        finite[i] = (func.f(xx) - func.f(xy)) / (2.0 * eps);
    }
    grad = finite;
}

template <typename Func>
void computeGradientC(const Func & func, const typename Func::InputType & x,
                      typename Func::JacobianType & grad, const double eps)
{
    const size_t DIM = x.rows();
    typename Func::JacobianType finite(DIM);
    typename Func::DualInputType xx = x.template cast<typename Func::DualScalar > ();
    for(size_t i = 0; i < DIM; i++)
    {
        xx[i] += Func::DualScalar(0,1) * eps;
        finite[i] = imag(func.f(xx)) / eps;
        xx[i] = x[i];
    }
    grad = finite;
}

template <typename Func>
double checkGradient(const Func & func, const typename Func::InputType & x,
                      const typename Func::JacobianType & grad, const double eps)
{
    const size_t DIM = x.rows();
    typename Func::JacobianType finite(DIM);
    for(size_t i = 0; i < DIM; i++)
    {
        typename Func::InputType xx = x;
        xx[i] += eps;
        typename Func::InputType xy = x;
        xy[i] -= eps;
        finite[i] = (func.f(xx) - func.f(xy)) / (2.0 * eps);
    }
    double error = static_cast<typename Func::JacobianType>((finite - grad)).norm() 
                 / static_cast<typename Func::JacobianType>((finite + grad)).norm();
    if (error > eps)
      std::cerr << "(checkGradient error)=" << error << " eps=" << eps << "\n";
    return error;
}

template <typename Func>
double checkGradientC(const typename Func::InputType & x,
                      const typename Func::JacobianType & grad, const double eps)
{
    const size_t DIM = x.rows();
    typename Func::JacobianType finite(DIM);
    computeGradientC<Func>(x, finite, eps);
    //double error = static_cast<Vector>((finite - grad)).norm() / static_cast<Vector>((finite + grad)).norm();
    size_t row, col;
    typename Func::JacobianType diff = (finite - grad);
    //double error = diff.norm();
    double error = diff.cwiseAbs().maxCoeff(&row, &col);
    if (error > 1e-13)
      std::cerr << "(checkGradientC error)=" << error
                << "(" << row << "," << col << ") g=" << grad(row,col) << " gest=" << finite(row,col)
                << " eps=" << eps << "\n";
    return error;
}

template <typename Func>
void computeHessian(const Func & func, const typename Func::InputType & x, 
                    typename Func::HessianType & hessian, const double eps)
{
    Assert(x.rows() == hessian.rows(), "hessian has wrong dimension (number of rows)");
    Assert(x.rows() == hessian.cols(), "hessian has wrong dimension (number of cols)");
    const size_t DIM = x.rows();
    typename Func::InputType xx = x;
    for(size_t i = 0; i < DIM; i++)
    {
        for(size_t j = 0; j < DIM; j++)
        {
            double f4 = func.f(xx);
            xx[i] += eps;
            xx[j] += eps;

            double f1 = func.f(xx);

            xx[j] -= eps;
            double f2 = func.f(xx);
            xx[j] += eps;
            xx[i] -= eps;
            double f3 = func.f(xx);

            hessian(i, j) = (f1 - f2 - f3 + f4) / (eps * eps);
            xx[i] = x[i];
            xx[j] = x[j];
        }
    }
}

}

