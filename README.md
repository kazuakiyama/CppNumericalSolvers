CppNumericalSolvers (C++11 implementation with MATLAB bindings)   
=================================================================

[![Build Status](https://api.travis-ci.org/PatWie/CppNumericalSolvers.svg?branch=master)](http://travis-ci.org/PatWie/CppNumericalSolvers)
[![Build Status](https://img.shields.io/github/release/PatWie/CppNumericalSolvers.svg)](https://github.com/PatWie/CppNumericalSolvers/releases)
[![Build Status](https://img.shields.io/github/issues/PatWie/CppNumericalSolvers.svg)](https://github.com/PatWie/CppNumericalSolvers/issues)

**Faster than Matlab-built-in** optimization algorithms, when gradient information is available.


Quick Introduction
-----------
- run `make install` to download and build dependencies
- run `make test` to verify the results by unit tests
- run `make main` to build cpp examples
- run `make.m` within MATLAB to build the MATLAB wrapper

Extensive Introduction
-----------

This repository contains several solvers implemented in C++11 using the [Eigen3][eigen3] library for all matrix computations. All implementations were written from scratch. 
You can use this library in **C++ and [Matlab][matlab]** in an comfortable way.

The library currently contains the following solvers:

- gradient descent solver
- conjugate gradient descent solver
- Newton descent solver
- BFGS solver
- L-BFGS solver
- L-BFGS-B solver

Additional helpful included functions are

- gradient check
- gradient approximation by finite differences
- Hessian matrix approximation by finite differences

See the unit-test for the minimization of the Rosenbrock function when starting in [-1.2,1]. By convention these solvers minimize a given objective function.

### Usage in C++
First, make sure that you have the [Eigen3][eigen3] in your include paths and your compiler is C++11-ready (e.g. use g++-4.8). To minimize a function it is sufficient to define the objective function as a C++11 functional:

	auto function_value = [] (const Vector &x) -> double {};
	auto gradient_value = [] (const Vector x, Vector &grad) -> void {}; // this definition is optional, but gainful

If you only implement the objective function then the library will compute a numerical approximation of the gradient (this is slow!). If you want to use `Newton Descent` without explicitly computing the gradient or Hessian matrix you should use `computeGradient`, `computeHessian` from `Meta.h` (check the examples).
To optimize (minimize) the defined objective function you need to select a solver and provide a initial guess:


	Vector x0;                  // initial guess
    x0 << ..foo ... ;
	BfgsSolver alg;          // choose a solver

	// solve without explicitly declaring a gradient
	alg.Solve(x0,function_value);

	// or solve with given gradient
	alg.Solve(x0,function_value, gradient_value);

	// or use the Newton method
	auto hessian_value = [&](const Vector x, Matrix & hes) -> void
	{
	    hes = Matrix::Zero(x.rows(), x.rows());
	    computeHessian(function_value, x, hes);
	};
	NewtonDescentSolver newton;      // or use newton descent
	newton.Solve(x0,function_value,gradient_value,hessian_value);


I encourage you to check you gradient once by 

	checkGradient(YOUR-OBJECTIVE-FUNCTION, X, A-SAMPLE-GRADIENT-VECTOR_FROM-YOUR-FUNCTIONAL);

to make sure that there are not mistakes in your objective or gradient function.

	// full sample
	int main(void) {

		// create function
		auto rosenbrock = [] (const Vector &x) -> double {
			const double t1 = (1-x[0]);
			const double t2 = (x[1]-x[0]*x[0]);
			return   t1*t1 + 100*t2*t2;
		};

		// create derivative of function
		auto Drosenbrock = [] (const Vector x, Vector &grad) -> void {
			grad[0]  = -2*(1-x[0])+200*(x[1]-x[0]*x[0])*(-2*x[0]);
			grad[1]  = 200*(x[1]-x[0]*x[0]);
		};

		// initial guess
		Vector x0(2);x0 << -1.2,1;
		// get derivative
		Vector dx0(2);
		Drosenbrock(x0,dx0);
		// verify gradient by finite differences
		checkGradient(rosenbrock,x0,dx0);

		// use solver (GradientDescentSolver,BfgsSolver,LbfgsSolver,LbfgsbSolver)
		BfgsSolver alg;
		alg.Solve(x0,rosenbrock,Drosenbrock);

		std::cout << std::endl<<std::endl<< x0.transpose();

		return 0;
	}

## tiny benchmark in C++
All algorithms start in [-1.2,1] and [15, 8] when minimizing the Rosenbrock objective function.

    [==========] Running 11 tests from 6 test cases.
    [----------] 2 tests from GradientDescentTest (85 ms total)
    [----------] 2 tests from ConjugateGradientTest (0 ms total)
    [----------] 2 tests from NewtonDescentTest (1 ms total)
    [----------] 2 tests from BfgsTest (0 ms total)
    [----------] 2 tests from LbfgsTest (0 ms total)
    [----------] 1 test from LbfgsbTest (11 ms total)
    [==========] 11 tests from 6 test cases ran. (97 ms total)



### Usage in MATLAB (experimental)
If something goes wrong MATLAB will tell it to you by simply crashing. Hence, it is an experimental version. 

- Download [Eigen][eigen3] and copy the `Eigen` folder into `matlab-bindings`. The compilation script will check if Eigen exists.
- Mex files are no fun! But if your configs are correct ("make.m" will use "-std=c++11") you can run `make` inside the `matlab-bindings` folder. This should create a file `cppsolver.mexw64` or `cppsolver.mexw32`.

To solve your problem just define the objective as a function like

	function y = rosenbrock( x )
	  t1 = (1 - x(1));
	  t2 = (x(2) - x(1) * x(1));
	  y  = t1 * t1 + 100 * t2 * t2; 
	end

an save it under `rosenbrock.m`. If everything is correct you should be able to minimize your function:

	[solution, value] = cppsolver([-1;2], @rosenbrock)
	% or if you specify a gradient:
	[solution, value] = cppsolver([-1;2], @rosenbrock, 'gradient', @rosenbrock_grad)


If you pass a gradient function the call of the mex-file, it will check the gradient function once in the initial guess to make sure you have no typos in your gradient function. You can skip this check by calling 

	[solution, value] = cppsolver([-1;2],@rosenbrock,'gradient',@rosenbrock_grad,'skip_gradient_check','true')


The default solver is `BFGS`. To specify another solver you can call

	[solution, value] = cppsolver([-1;2],@rosenbrock,'gradient',@rosenbrock_grad,'solver','THE-SOLVER-TYPE')

There is an example in the folder `matlab-bindings`. 

## tiny benchmark in Matlab

All algorithms start in [-1.2,1] when minimizing the Rosenbrock objective function.

| algorithm     | elapsed time      | available information         | library         |
| ------------- | -------------     | -------------                 | -------------   |
| L-BFGS        | 0.005792 seconds. | (function, gradient)          | this library    |
| L-BFGS        | 0.011297 seconds. | (function)                    | this library    |
| NEWTON        | 0.004525 seconds. | (function, gradient, hessian) | this library    |
| BFGS          | 0.004113 seconds. | (function, gradient)          | this library    |
| conj. grad    | 0.080964 seconds. | (function, gradient)          | this library    |
| fminsearch    | 0.007286 seconds. | (function)                    | built-in Matlab |
| minunc        | 0.028901 seconds. | (function)                    | built-in Matlab |
| minimize      | 0.057358 seconds. | (function, gradient)          | (C.Rasmussen)   |



Notice, that fminsearch probably use caches if you run this file multiple times. If you do not specify at least one bound for the L-BFGS-B algorithm, the mex-file will use call the L-BFGS-algorithm instead.

# License

	Copyright (c) 2014-2015 Patrick Wieschollek
	Url: https://github.com/PatWie/CppNumericalSolvers

	Permission is hereby granted, free of charge, to any person obtaining a copy
	of this software and associated documentation files (the "Software"), to deal
	in the Software without restriction, including without limitation the rights
	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
	copies of the Software, and to permit persons to whom the Software is
	furnished to do so, subject to the following conditions:

	The above copyright notice and this permission notice shall be included in all
	copies or substantial portions of the Software.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
	SOFTWARE.



[eigen3]: http://eigen.tuxfamily.org/
[matlab]: http://www.mathworks.de/products/matlab/
[matlab-binding]: https://github.com/PatWie/CppNumericalSolvers/archive/matlab.zip