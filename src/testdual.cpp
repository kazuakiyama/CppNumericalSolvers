#include "gtest/gtest.h"
#include <iostream>
#include <complex>
#include <functional>
#include <list>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/StdVector>
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/DualNum>

#define PRECISION (1e-15)

typedef std::complex<double> complexd;
typedef std::complex<float> complexf;

using namespace Eigen;

TEST (DualTest, Construct)
{
}

template <typename DUALTYPE, typename Scalar>
void construct()
{
  Matrix<DUALTYPE,1,1> M;
  M << (Scalar)1.1;
  EXPECT_EQ(M(0), (Scalar)1.1);
}

template <typename DUALTYPE, typename Scalar>
void equality()
{
  Matrix<DUALTYPE,2,2> a = Matrix<DUALTYPE,2,2>::Random();
  Matrix<DUALTYPE,2,2> b = a;
  EXPECT_EQ(a,b);
}

template <typename DUALTYPE, typename Scalar>
void compare()
{
  Matrix<DUALTYPE,1,1> a;
  
}

template <typename DUALTYPE, typename Scalar>
void transcendental()
{
  Matrix<DUALTYPE,1,1> a;
  
}

template <typename DUALTYPE, typename Scalar>
void arithmetic()
{
  Matrix<DUALTYPE,1,1> a;
  
}

#define TESTFUNC(func) \
  TEST (Duald, func) { func<Duald, double>(); } \
  TEST (Dualf, func) { func<Dualf, float>(); } \
  TEST (Dualcd, func) { func<Dualcd, complexd>(); } \
  TEST (Dualdf, func) { func<Dualcf, complexf>(); }

TESTFUNC(construct)
TESTFUNC(equality)
TESTFUNC(compare)
TESTFUNC(arithmetic)
TESTFUNC(transcendental)

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  std::cout.precision(20);
  std::cerr.precision(20);
  return RUN_ALL_TESTS();
}
