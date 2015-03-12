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

#ifndef LBFGSBSOLVER_H_
#define LBFGSBSOLVER_H_
#include "ISolver.hpp"
#include <vector>
namespace pwie
{

/* coded from scratch !!!
 * based on the paper
 * A LIMITED MEMORY ALGORITHM FOR BOUND CONSTRAINED OPTIMIZATION
 * (Byrd, Lu, Nocedal, Zhu)
 */
template <typename Func>
class LbfgsbSolver : public ISolver<Func>
{
    typedef typename Func::Scalar Scalar;
    typedef typename Func::InputType InputType;
    typedef typename Func::JacobianType JacobianType;
    typedef typename ISolver<Func>::HessianType HessianType;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixX;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorX;

private:
    MatrixX W, M;
    Scalar theta;
    InputType _ub;
    InputType _lb;
    int _DIM;
public:
    using ISolver<Func>::settings;
    using ISolver<Func>::_functor;
    using ISolver<Func>::LineSearch;
private:
    /// <summary>
    /// find cauchy point in x
    /// </summary>
    /// <parameter name="x">start in x</parameter>
    void GetGeneralizedCauchyPoint(const InputType & x, JacobianType & g, InputType & x_cauchy, VectorX & c);
    /// <summary>
    /// find valid alpha for (8.5)
    /// </summary>
    /// <parameter name="x_cp">cauchy point</parameter>
    /// <parameter name="du">unconstrained solution of subspace minimization</parameter>
    /// <parameter name="FreeVariables">flag (1 if is free variable and 0 if is not free variable)</parameter>
    Scalar FindAlpha(const InputType & x_cp, const VectorX & du, const std::vector<int> & FreeVariables);
    /// <summary>
    /// direct primal approach
    /// </summary>
    /// <parameter name="x">start in x</parameter>
    void SubspaceMinimization(InputType & x_cauchy, InputType & x, const VectorX & c,
                              JacobianType & g, InputType & SubspaceMin);

public:
    LbfgsbSolver(const Func & func);
    virtual ~LbfgsbSolver() {};
    void internalSolve(InputType & x0);
};

} /* namespace pwie */

#include "CppNumericalSolvers/src/LbfgsbSolver.cpp"

#endif /* LBFGSBSOLVER_H_ */
