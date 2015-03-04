#include "gtest/gtest.h"
#include <iostream>
#include <complex>
#include <type_traits>
#include <functional>
#include <list>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/StdVector>
#include <unsupported/Eigen/MatrixFunctions>
#include <unsupported/Eigen/DualNum>

#define PRECISION (1e-15)

typedef std::complex<float> complexf;
typedef std::complex<double> complexd;
typedef std::complex<long double> complexld;

using namespace Eigen;
using namespace cxxduals;

#if 0
template <typename _Tp> _Tp random();
template <> dual<std::complex<float> > random<dual<std::complex<float> > >() {
  return dual<std::complex<float> >({(float)drand48(), (float)drand48()});
}
template <> dual<std::complex<double> > random<dual<std::complex<double> > >() {
  std::complex<double> cx(drand48(), drand48());
  return dual<std::complex<double> >(cx);
}
#endif

TEST (DualTest, Construct)
{
  dual<double> a;
}

template <typename DUALTYPE, typename Scalar>
void construct()
{
  Matrix<DUALTYPE,1,1> M;
  M << (Scalar)1.1;
  EXPECT_EQ(M(0), (Scalar)1.1);

  bool IsInteger = NumTraits<DUALTYPE>::IsInteger;
  bool IsSigned = NumTraits<DUALTYPE>::IsSigned;
  bool IsComplex = NumTraits<DUALTYPE>::IsComplex;
  bool RequireInitialization = NumTraits<DUALTYPE>::RequireInitialization;
  int ReadCost = NumTraits<DUALTYPE>::ReadCost;
  int AddCost = NumTraits<DUALTYPE>::AddCost;
  int MulCost = NumTraits<DUALTYPE>::MulCost;
  std::stringstream s;
  s << "ICSR" << IsInteger << IsComplex << IsSigned << RequireInitialization
    << "R" << ReadCost
    << "A" << AddCost
    << "M" << MulCost;
  //std::cout << s.str() << "\n";
}

template <typename DUALTYPE, typename Scalar>
void equality()
{
  Matrix<DUALTYPE,2,2> a(Matrix<DUALTYPE,2,2>::Ones());
  Matrix<DUALTYPE,2,2> b = a;
  Matrix<DUALTYPE,2,2> c = Matrix<DUALTYPE,2,2>::Identity();
  a = b * 1.0;
  EXPECT_EQ(a,b);
  a = b * c;
  EXPECT_EQ(a,b);
  a = 1.0 * b;
  EXPECT_EQ(a,b);
  b = 1.0 * b;
  EXPECT_EQ(a,b);
  b = Matrix<DUALTYPE,2,2>::Zero();
  EXPECT_NE(a,b);
  a.setZero();
  EXPECT_EQ(a,b);
}

template <typename DUALTYPE, typename Scalar>
void compare()
{
  using std::abs;
  Matrix<DUALTYPE,2,2> a = Matrix<DUALTYPE,2,2>::Random();
  Matrix<DUALTYPE,2,2> b = a;
  Matrix<DUALTYPE,2,2> c = a;
  a = b * 1.0;
  a = 2.0 * b;
  b = a / 2.0;
  EXPECT_LT(abs((a-b*2).norm()), 
            abs(std::numeric_limits<typename DUALTYPE::basic_value_type>::epsilon()));
  a = b + c;
  a = b * 2.0;
  a = 2.0 * b;
  a = b - c;
  a = b * c;
  a = b / 2;
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

template <typename T1, typename T2>
void casting()
{
  typedef Matrix<DualNum<T1>, 2, 2> dt1;
  typedef Matrix<DualNum<T2>, 2, 2> dt2;
  dt1 m1 = dt1::Zero();
  dt1 m1p = dt1::Ones();
  dt2 m2 = dt2::Ones();
  m1 = m2.template cast<DualNum<T1> >();
  EXPECT_TRUE(m1 == m1p);
  m2 = dt2::Identity();
  m1p = dt1::Identity();
  m1 = m2.template cast<DualNum<T1> >();
  EXPECT_TRUE(m1 == m1p);
}

#define TESTFUNC_REAL(func) \
  TEST (Dualf, func) { func<Dualf, float>(); } \
  TEST (Duald, func) { func<Duald, double>(); } \
  TEST (Dualld, func) { func<Dualld, long double>(); }

#define TESTFUNC(func) \
  TEST (Dualf, func) { func<Dualf, float>(); } \
  TEST (Duald, func) { func<Duald, double>(); } \
  TEST (Dualld, func) { func<Dualld, long double>(); } \
  TEST (Dualcf, func) { func<Dualcf, complexf>(); } \
  TEST (Dualcd, func) { func<Dualcd, complexd>(); } \
  TEST (Dualcld, func) { func<Dualcld, complexld>(); }

#define TEST_INTER_TYPE(func) \
  TEST (float_double, func) { func<float,double>(); } \
  TEST (double_float, func) { func<double,float>(); } \
  TEST (Duald_Dualf, func) { func<Duald,Dualf>(); }

TESTFUNC(construct)
TESTFUNC(equality)
TESTFUNC_REAL(compare)
TESTFUNC(arithmetic)
TESTFUNC(transcendental)

TEST_INTER_TYPE(casting)

int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  std::cout.precision(20);
  std::cerr.precision(20);
  return RUN_ALL_TESTS();
}
