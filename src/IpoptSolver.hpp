/**
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

#ifndef IPOPTSOLVER_H_
#define IPOPTSOLVER_H_

#include <assert.h>
#include "IpIpoptApplication.hpp"
#include "CppNumericalSolvers/src/ISolver.hpp"

using namespace Ipopt;

namespace pwie {

template <class Func>
class IpoptSolver : public TNLP, public ISolver<Func> {
  typedef typename Func::Scalar Scalar;
  typedef typename Func::InputType InputType;
  typedef typename Func::JacobianType JacobianType;

protected:
  
  using ISolver<Func>::settings;
  InputType _x;
  Referencer _dumbPtr; // work around for Ipopt's stupid SmartPtr

public:
  IpoptSolver() { AddRef(&_dumbPtr); }
  virtual ~IpoptSolver() { ReleaseRef(&_dumbPtr); }

  virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, 
                            Index& nnz_h_lag, IndexStyleEnum& index_style)
  {
    n = Func::InputsAtCompileTime;
    m = 0; // TODO- detect if Func has constraints defined
    nnz_jac_g = 0;
    nnz_h_lag = n * n;
    index_style = TNLP::C_STYLE;
    return true;
  }

  virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                               Index m, Number* g_l, Number* g_u)
  {
    // the variables bounds
    InputType lower = getLowerBound(n);
    InputType upper = getUpperBound(n);

    for (Index i = 0; i < n; i++) {
      x_l[i] = lower[i];
      x_u[i] = upper[i];
    }

#if 0
    // TODO: detect if Func has constraints
    for (Index jj = 0; jj < m; jj++) {
      g_l[jj] = -2e19;
      g_u[jj] = 2e19;
    }
#endif
    return true;
  }

  virtual bool get_starting_point(Index n,
                                  bool init_x, Number* x,
                                  bool init_z, Number* z_L, Number* z_U,
                                  Index m, bool init_lambda,
                                  Number* lambda)
  {
    assert(init_x == true);
    assert(init_z == false);
    assert(init_lambda == false);

    // initialize to the given starting point
    for (Index i = 0; i < n; i++)
      x[i] = _x[i];

    return true;
  }

  // objective function
  virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
  {
    assert(n == (Index)Func::InputsAtCompileTime);
    Eigen::Map<const InputType> xmap(x, n);
    // TODO: change objective function argument to accept MatrixBase, get rid of xx, pass xmap
    InputType xx = xmap;
    obj_value = this->f(xx);
    return true;
  }

  // gradient of objective function
  virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
  {
    assert(n == (Index)Func::InputsAtCompileTime);
    Eigen::Map<const InputType> xmap(x, n);
    Eigen::Map<JacobianType> gmap(grad_f, n);
    // TODO: change objective function argument to accept MatrixBase:
    //  get rid of xx,gg, pass xmap,gmap
    InputType xx = xmap;
    JacobianType gg(n);
    this->gradient(xx, gg);
    gmap = gg;
    return true;
  }

  // constraints function
  virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
  {
    assert(n == (Index)Func::InputsAtCompileTime);
    // TODO
    return false;
  }

  // jacobian
  virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
                          Index m, Index nele_jac, Index* iRow, Index *jCol,
                          Number* values)
  {
    assert(n == (Index)Func::InputsAtCompileTime);
    assert(m == (Index)0);
    // TODO
    return false;
  }

  // hessian
  virtual bool eval_h(Index n, const Number* x, bool new_x,
                      Number obj_factor, Index m, const Number* lambda,
                      bool new_lambda, Index nele_hess, Index* iRow,
                      Index* jCol, Number* values)
  {
    // TODO
    return false;
  }

  void finalize_solution(SolverReturn status,
                         Index n, const Number* x, const Number* z_L,
                         const Number* z_U, Index m, const Number* g,
                         const Number* lambda, Number obj_value,
                         const IpoptData* ip_data, IpoptCalculatedQuantities* ip_cq)
  {
    assert(n == (Index)Func::InputsAtCompileTime);
    assert(m == (Index)0);
    for (Index i=0; i<n; i++)
      _x[i] = x[i];
  }

  virtual void internalSolve(InputType & x0) {
    _x = x0;

    // because (this) already has a reference _dumbPtr, it wont be deleted when
    // mynlp goes out of scope.
    SmartPtr<TNLP> mynlp = this;

    // Use the Ipopt application factory
    SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

    // Adjust options
    // Note: The following choices are only examples, they might not be
    //       suitable for your optimization problem.
    app->Options()->SetNumericValue("tol", 1e-10);
    app->Options()->SetNumericValue("acceptable_tol", 1e-8);
    app->Options()->SetIntegerValue("max_iter", settings.maxIter);
    app->Options()->SetIntegerValue("print_level", 0);
    //app->Options()->SetNumericValue("print_frequency_time", 2);
    //app->Options()->SetNumericValue("print_timing_statistics", yes);
    //app->Options()->SetStringValue("derivative_test", "first-order");
    app->Options()->SetStringValue("linear_solver", "ma27"); // 57,77,86,97,pardiso,wsmp,mumps,
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    app->Options()->SetStringValue("limited_memory_update_type", "bfgs");
    app->Options()->SetIntegerValue("limited_memory_max_history", 10);

    // Intialize the IpoptApplication and process the options
    ApplicationReturnStatus status;
    status = app->Initialize();
    if (status != Solve_Succeeded) {
      std::cerr << "\n*** IpoptSolver: Error during initialization!\n";
      return;
    }

    // Ipopt solves the problem
    status = app->OptimizeTNLP(this);

    if (status == Solve_Succeeded) {
      //std::cerr << "*** IpoptSolver: Problem solved!\n";
      //std::cerr << x0;
    }
    else {
      std::cerr << "*** IpoptSolver: Problem FAILED!\n";
    }

    x0 = _x;
  }
};

}

#endif
