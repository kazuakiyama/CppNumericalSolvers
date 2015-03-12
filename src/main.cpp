/*
 * main.cpp
 *      Author: Patrick Wieschollek
 */
#include <iostream>
#include <memory>
#include "CppNumericalSolvers.hpp"

int example000(void);
int example001(void);
int example002(void);
int example003(void);
int example004(void);

using namespace pwie;

int main(void)
{
    example000();
    example001();
    example002();
    example003();
    //example004();

    return 0;
}

static const char * solver_names[] = {
  "gradient descent",
  "newton descent",
  "conjugate gradient",
  "asa conjugate gradient",
  "BFGS",
  "L-BFGS",
  "L-BFGS-B",
  "Ipopt",
};

class example000Func {
public:
    typedef double Scalar;
    enum {
      InputsAtCompileTime = 2,
      ValuesAtCompileTime = 1
    };
    typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
  //typedef Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType; todo - use this
    typedef Eigen::Matrix<Scalar,InputsAtCompileTime,ValuesAtCompileTime> JacobianType;
    typedef Eigen::Matrix<Scalar,InputsAtCompileTime,InputsAtCompileTime> HessianType;

    // create function
    template <typename eT>
    eT f(const Eigen::Matrix<eT,InputsAtCompileTime,1> & x) const {
        const eT t1 = (1.0 - x[0]);
        const eT t2 = (x[1] - x[0] * x[0]);
        return t1 * t1 + 100.0 * t2 * t2;
    }

    // create derivative of function
    void gradient(const InputType & x, JacobianType & grad) const {
      grad[0]  = -2.0 * (1.0 - x[0]) + 200.0 * (x[1] - x[0] * x[0]) * (-2.0 * x[0]);
      //grad[0]  = -400.0 * x[0] * (x[1] - x[0] * x[0]);
      grad[1]  = 200.0 * (x[1] - x[0] * x[0]);
    }
};

struct example000FuncBADG : public example000Func {
    // create derivative of function
    void gradient(const InputType & x, JacobianType & grad) const {
        // wrong derivative of function
        grad[0]  = -2.0 * (1.0 - x[0]) + 200 * (x[1] - x[0] * x[0]) * (2.0 * x[0]); // <<-- last term wrong sign
        grad[1]  = 200.0 * (x[1] - x[0] * x[0]);
    }
};

int example000(void)
{
    std::cout << "------------------------------------------" << std::endl;
    std::cout << "example000 checkgradient:    " << std::endl
              << "------------------------------------------" << std::endl;

    // check your gradient!
    example000Func f_;
    Functor<example000Func> f(f_);

    // random x
    example000Func::InputType x(2);
    x << 15, 8;

    // get "derivatives"
    example000Func::JacobianType dx(2);
    f.gradient(x, dx);

    // verify gradient by finite differences
    std::cout << "check correct gradient: ";

    double gerr = f.checkGradient(x);
    if(gerr < 1e-13)
    {
        std::cout << "gradient probably correct" << std::endl;
    }

    std::cout << "check wrong gradient:   ";
    example000FuncBADG fb_;
    Functor<example000Func> fb(fb_);
    gerr = fb.checkGradient(x);
    if (gerr > 1e-13)
    {
        std::cout << "gradient probably wrong" << std::endl;
    }

    return 0;
}

class example001Func {
public:
    typedef double Scalar;
    enum {
      InputsAtCompileTime = 2,
      ValuesAtCompileTime = 1
    };
    typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
    typedef Eigen::Matrix<Scalar,InputsAtCompileTime,ValuesAtCompileTime> JacobianType;
    typedef Eigen::Matrix<Scalar,InputsAtCompileTime,InputsAtCompileTime> HessianType;

    // create function
    template <typename eT>
    eT f(const Eigen::Matrix<eT,InputsAtCompileTime,1> & x) const {
        const eT t1 = (1.0 - x[0]);
        const eT t2 = (x[1] - x[0] * x[0]);
        return   t1 * t1 + 100.0 * t2 * t2;
    }
};

int example001(void)
{
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "example001 optimization without gradient:    " << std::endl
            << "------------------------------------------" << std::endl;

  // initial guess
  example001Func::InputType x0(2);
  example001Func f;
  for (int ii = 0; ii < SOLVER_COUNT; ii++) {
    x0 << 3, 2;
    // use solver (GradientDescentSolver,BfgsSolver,LbfgsSolver,LbfgsbSolver)
    std::unique_ptr<ISolver<example001Func> > g(getSolver<example001Func>(f, (solver_id)ii));
    try {
      g->solve(x0);
    } catch (std::exception & ex) {
      std::cout << "excepted\n";
    }

    std::cout << "result:    " << x0.transpose() << " : " << solver_names[ii];
    std::cout << std::endl;
    std::cout << "should be: " << "1 1" << std::endl;
  }

  return 0;
}

int example002(void)
{
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "example002 optimization with  gradient:    " << std::endl
            << "------------------------------------------" << std::endl;

  // initial guess
  example000Func::InputType x0(2);
  example000Func f;
  for (int ii = 0; ii < SOLVER_COUNT; ii++) {
    // use solver (GradientDescentSolver,BfgsSolver,LbfgsSolver,LbfgsbSolver)
    std::unique_ptr<ISolver<example000Func> > g(getSolver<example000Func>(f, (solver_id)ii));
    x0 << 15, 8;
    x0 << -1.2,1;

    try {
      g->solve(x0);
    } catch (std::exception & ex) {
      std::cout << "excepted\n";
    }

    std::cout << "result:    " << x0.transpose() << " : " << solver_names[ii];
    std::cout << std::endl;
    std::cout << "should be: " << "1 1" << std::endl;
  }
  return 0;
}

class example003Func {
public:
    typedef double Scalar;
    enum {
      InputsAtCompileTime = 2,
      ValuesAtCompileTime = 1
    };
    typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
    typedef Eigen::Matrix<Scalar,InputsAtCompileTime,ValuesAtCompileTime> JacobianType;
    typedef Eigen::Matrix<Scalar,InputsAtCompileTime,InputsAtCompileTime> HessianType;

    // create derivative of function
    void gradient(const InputType & x, JacobianType & grad) const {
        grad = JacobianType(x.rows());
        grad = 2 * A.transpose() * (A * x) - 2 * A.transpose() * y;
    }
    example003Func() {
        A << 0.2, 0.25, 0.4, 0.5, 0.4, 0.25;
        y << 0.9, 1.7, 1.2;
    }

    template <typename eT>
    eT f(const Eigen::Matrix<eT,InputsAtCompileTime,1> & x) const {
        Eigen::Matrix<eT,3,1> tmp = A * x - y;
        return tmp.dot(tmp);
    }
private:
    Eigen::Matrix<double,3,2> A;
    Eigen::Matrix<double,3,1> y;
};

int example003(void)
{

  std::cout << "------------------------------------------" << std::endl;
  std::cout << "example003 optimization least squares:    " << std::endl
            << "------------------------------------------" << std::endl;

  // initial guess
  example003Func::InputType x0(2);
  example003Func f;
  for (int ii = 0; ii < SOLVER_COUNT; ii++) {
    // use solver (GradientDescentSolver,BfgsSolver,LbfgsSolver,LbfgsbSolver)
    std::unique_ptr<ISolver<example003Func> > g(getSolver<example003Func>(f, (solver_id)ii));
    x0 << 15, 8;

    try {
      g->solve(x0);
    } catch (std::exception & ex) {
      std::cout << "excepted\n";
    }
    
    std::cout << "result:    " << x0.transpose() << " : " << solver_names[ii];
    std::cout << std::endl;
    std::cout << "should be:  " << "1.7 2.08" << std::endl;
  }

  return 0;
}

// super-simple case: y = x^2
class example004Func {
public:
  typedef double Scalar;
  enum {
    InputsAtCompileTime = 1,
    ValuesAtCompileTime = 1
  };
  typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
  typedef Eigen::Matrix<Scalar,InputsAtCompileTime,ValuesAtCompileTime> JacobianType;
  typedef Eigen::Matrix<Scalar,InputsAtCompileTime,InputsAtCompileTime> HessianType;

  // create derivative of function
  void gradient(const InputType & x, JacobianType & grad) const {
    grad = JacobianType(x.rows());
    grad[0] = -sin(x[0]);
  }

  template <typename eT>
  eT f(const Eigen::Matrix<eT,InputsAtCompileTime,1> & x) const {
    return cos(x[0]);
    //return x * x;
  }
private:
};

int example004(void)
{
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "example004 function checks:    " << std::endl
            << "------------------------------------------" << std::endl;

  // initial guess
  example004Func::InputType x0(1);
  double y0;
  example004Func::JacobianType gradFD(1), gradD(1), grad(1);
  example004Func func;
  Functor<example004Func> functor(func);
  pwie::LbfgsbSolver<example004Func> diff(func);

  diff.settings.verbosity = 1;

  for (int ii = -100; ii < 100; ii++) {
    x0[0] = ii/10.0;
    y0 = functor.f(x0);
    functor.gradient(x0, grad);
    functor.gradientFiniteDiff(x0, gradFD);
    functor.gradientDual(x0, gradD);
    std::cout << "x=" << x0 << " y=" << y0 << " f.d. grad=" << gradFD
              << " dual=" << gradD << " grad=" << grad << "\n";
  }

  diff.solve(x0);
  std::cout << "xmin @ " << x0 << "\n";

  return 0;
}
