#include "gtest/gtest.h"
#include <iostream>
#include <functional>
#include <list>
#include "Meta.hpp"
#include "BfgsSolver.hpp"
#include "LbfgsSolver.hpp"
#include "LbfgsbSolver.hpp"
#include "GradientDescentSolver.hpp"
#include "ConjugateGradientSolver.hpp"
#include "NewtonDescentSolver.hpp"

#define PRECISION 1e-7

class rosenbrock {
public:
    typedef double Scalar;
    enum {
      InputsAtCompileTime = 2,
      ValuesAtCompileTime = 1
    };
    typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
    typedef Eigen::Matrix<Scalar,InputsAtCompileTime,ValuesAtCompileTime> JacobianType;
    typedef Eigen::Matrix<Scalar,InputsAtCompileTime,InputsAtCompileTime> HessianType;

    // function
    template <typename eT>
    eT f(const Eigen::Matrix<eT,InputsAtCompileTime,1> & x) const {
        const eT t1 = (1.0 - x[0]);
        const eT t2 = (x[1] - x[0] * x[0]);
        return   t1 * t1 + 100.0 * t2 * t2;
    }
};
class rosenbrockD : public rosenbrock {
public:
    // derivative of function
    void gradient(const InputType & x, JacobianType & grad) const {
        grad[0]  = -2.0 * (1.0 - x[0]) + 200 * (x[1] - x[0] * x[0]) * (-2.0 * x[0]);
        grad[1]  = 200.0 * (x[1] - x[0] * x[0]);
    }

};
class rosenbrockDD : public rosenbrockD {
public:
    // hessian of function
    void hessian(const InputType & x, HessianType & hes) const {
        hes = HessianType::Zero(x.rows(), x.rows());
        hes(0, 0) = 1200.0 * x[0] * x[0] - 400.0 * x[1] + 1.0;
        hes(1, 0) = -400.0 * x[0];
        hes(0, 1) = -400.0 * x[0];
        hes(1, 1) = 200.0;
    }
};


#define SOLVE_0stOrder(sol,func,a,b,fx)         \
  func::InputType x(2);x(0) = a;x(1) = b;       \
  func f; sol<func> solver(f); solver.solve(x); \
  EXPECT_NEAR(fx, f.f(x), PRECISION);
#define SOLVE_1stOrder(sol,func,a,b,fx)         \
  func::InputType x(2);x(0) = a;x(1) = b;       \
  func##D f;                                    \
  sol<func##D> solver(f);                       \
  solver.solve(x);                              \
  EXPECT_NEAR(fx, f.f(x), PRECISION);
#define SOLVE_2ndOrder(sol,func,a,b,fx)          \
  func::InputType x(2);x(0) = a;x(1) = b;        \
  func##DD f;                                    \
  sol<func##DD> solver(f);                       \
  solver.solve(x);                               \
  EXPECT_NEAR(fx, f.f(x), PRECISION);

TEST(GradientDescentTest, RosenbrockFar)   {
  SOLVE_1stOrder(pwie::GradientDescentSolver, rosenbrock, 15.0, 8.0, 0.0)
    }
TEST(GradientDescentTest, DISABLED_RosenbrockNear)  {
  SOLVE_1stOrder(pwie::GradientDescentSolver, rosenbrock, -1.2, 1.0, 0.0)
    }
TEST(GradientDescentTest, RosenbrockNearish)  {
  SOLVE_1stOrder(pwie::GradientDescentSolver, rosenbrock, 1.0, 3.0, 0.0)
    }
TEST(GradientDescentTest, RosenbrockNearFar)  {
  SOLVE_1stOrder(pwie::GradientDescentSolver, rosenbrock, -1.2, 100.0, 0.0)
    }

TEST(ConjugateGradientTest, RosenbrockFar) {
  SOLVE_1stOrder(pwie::ConjugateGradientSolver, rosenbrock, 15.0, 8.0, 0.0) }
TEST(ConjugateGradientTest, RosenbrockNear){
  SOLVE_1stOrder(pwie::ConjugateGradientSolver, rosenbrock, -1.2, 1.0, 0.0) }
TEST(ConjugateGradientTest, DISABLED_RosenbrockNearFar)  {
  SOLVE_1stOrder(pwie::ConjugateGradientSolver, rosenbrock, -1.2, 100.0, 0.0)
    }

TEST(NewtonDescentTest, RosenbrockFar)     {
  SOLVE_2ndOrder(pwie::NewtonDescentSolver, rosenbrock, 15.0, 8.0, 0.0) }
TEST(NewtonDescentTest, RosenbrockNear)    {
  SOLVE_2ndOrder(pwie::NewtonDescentSolver, rosenbrock, 1.0, 3.0, 0.0) }

TEST(BfgsTest, RosenbrockFar)              {
  SOLVE_1stOrder(pwie::BfgsSolver, rosenbrock, 15.0, 8.0, 0.0) }
TEST(BfgsTest, RosenbrockNear)             {
  SOLVE_1stOrder(pwie::BfgsSolver, rosenbrock, 1.0, 3.0, 0.0) }

TEST(LbfgsTest, RosenbrockFar)             {
  SOLVE_1stOrder(pwie::LbfgsSolver, rosenbrock, 15.0, 8.0, 0.0) }
TEST(LbfgsTest, RosenbrockNear)            {
  SOLVE_1stOrder(pwie::LbfgsSolver, rosenbrock, 1.0, 3.0, 0.0) }

TEST(LbfgsbTest, RosenbrockFar)            { 
  SOLVE_1stOrder(pwie::LbfgsbSolver, rosenbrock, 15.0, 8.0, 0.0) }
TEST(LbfgsbTest, RosenbrockNear)           { 
  SOLVE_1stOrder(pwie::LbfgsbSolver, rosenbrock, -1.2, 1.0, 0.0) }
TEST(LbfgsbTest, RosenbrockNearFar)  {
  SOLVE_1stOrder(pwie::LbfgsbSolver, rosenbrock, -1.2, 100.0, 0.0)
    }

int main (int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
