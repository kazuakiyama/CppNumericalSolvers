#include "gtest/gtest.h"
#include <iostream>
#include <functional>
#include <list>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/StdVector>
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/DualNum>

#define PRECISION (1e-9)

using namespace Eigen;

TEST (DualTest, Construct)
{
  Matrix<double,1,1> M;
  M << 1.1;
  EXPECT_NEAR(M(0), 1.1, PRECISION);
  EXPECT_EQ(M(0), 1.1);

  DualNum<double> x = 1.1;
  EXPECT_EQ(x.real(), 1.1);
  EXPECT_EQ(real(x), 1.1);
  EXPECT_EQ(x.dual(), 0.0);
  EXPECT_EQ(dual(x), 0.0);

  DualNum<double> z(1.1);
  EXPECT_EQ(z.real(), 1.1);
  EXPECT_EQ(real(z), 1.1);
  EXPECT_EQ(z.dual(), 0.0);
  EXPECT_EQ(dual(z), 0.0);

  DualNum<double> y(1.1, 2.2);
  EXPECT_EQ(y.real(), 1.1);
  EXPECT_EQ(real(y), 1.1);
  EXPECT_EQ(y.dual(), 2.2);
  EXPECT_EQ(dual(y), 2.2);

  DualNum<double> w = {1.1, 2.2};
  EXPECT_EQ(w.real(), 1.1);
  EXPECT_EQ(w.dual(), 2.2);

  DualNum<double> a{1.1, 2.2};
  EXPECT_EQ(a.real(), 1.1);
  EXPECT_EQ(a.dual(), 2.2);
}

template <typename DUALTYPE>
void compare()
{
  typedef typename DUALTYPE::Scalar Scalar;
  DUALTYPE d{1.1, 2.2};
  DUALTYPE e{1.1, 2.2};
  DUALTYPE f{1.1, 2.21};
  DUALTYPE g{2.2, 2.21};
  DUALTYPE h{-2.2, -2.21};
  DUALTYPE j{3, 0};
  DUALTYPE k{9, 6};
  Matrix<DUALTYPE,1,1> a;
  a << d;
  EXPECT_EQ(a[0], d);
  EXPECT_EQ(e, d);
  EXPECT_EQ(e, f);
  EXPECT_NE(f, g);
  EXPECT_NE(f, h);
  EXPECT_NE(f, 1.2);
  EXPECT_NE(1.2, f);
  EXPECT_EQ(a[0], (Scalar)1.1);
  EXPECT_EQ((Scalar)1.1, a[0]);
  EXPECT_GT(1.2, a[0]);
  EXPECT_GT(a[0],1.0);
  EXPECT_GT(g,f);
  EXPECT_LT(1.0, a[0]);
  EXPECT_LT(a[0],(Scalar)1.2);
  EXPECT_LT(f,g);
  EXPECT_LE(a[0],(Scalar)1.1);
  EXPECT_LE((Scalar)1.1,a[0]);
  EXPECT_LE(f,g);
  EXPECT_LE(f,e);
  EXPECT_GE(a[0],(Scalar)1.1);
  EXPECT_GE((Scalar)1.1,a[0]);
  EXPECT_GE(g,f);
  EXPECT_GE(e,f);
  EXPECT_TRUE(max(f,2.2) == (Scalar)2.2);
  EXPECT_TRUE(max(2.2,f) == (Scalar)2.2);
  EXPECT_TRUE(min(f,(Scalar)2.2) == f);
  EXPECT_TRUE(min((Scalar)2.2,f) == f);
  EXPECT_TRUE(max(f,g) == g);
  EXPECT_TRUE(max(g,f) == g);
  EXPECT_TRUE(min(f,g) == f);
  EXPECT_TRUE(min(g,f) == f);
  EXPECT_EQ(pow(j,2.0), k);
  EXPECT_EQ(sqrt(k), j);
  EXPECT_EQ(abs(h), g);
  EXPECT_EQ(abs(-h), g);
  EXPECT_EQ(-abs(h), h);
  EXPECT_EQ(abs(h), (Scalar)2.2);
  EXPECT_EQ(+abs(h), (Scalar)2.2);
  EXPECT_EQ(-abs(h), (Scalar)-2.2);
}

TEST (Duald, Comparison)
{
  compare<Duald>();
}

TEST (Dualf, Comparison)
{
  compare<Dualf>();
}

TEST (Dualcd, Comparison)
{
  compare<Dualcd>();
}

TEST (Dualcf, Comparison)
{
  compare<Dualcf>();
}

TEST (Duald, Transcendental)
{
  
}

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
