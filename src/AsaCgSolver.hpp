/**
 * Copyright (c) 2014 Patrick Wieschollek
 * Copyright (c) 2014 Michael Tesch
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

#ifndef ASACGSOLVER_H_
#define ASACGSOLVER_H_
#include "ISolver.hpp"
namespace pwie
{

template <typename Func>
class AsaCgSolver : public ISolver<Func>
{
    typedef typename Func::Scalar Scalar;
    typedef typename Func::InputType InputType;
    typedef typename Func::JacobianType JacobianType;
    typedef typename ISolver<Func>::HessianType HessianType;

    InputType _upper;
    InputType _lower;
public:
    AsaCgSolver();
    void internalSolve(InputType & x0);
};

} /* namespace pwie */

#include "CppNumericalSolvers/src/AsaCgSolver.cpp"

#endif /* ASACGSOLVER_H_ */