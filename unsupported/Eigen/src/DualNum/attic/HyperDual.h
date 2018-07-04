// -*-c++-*-
// This file is part of Eigen, a lightweight C++ template library
// for linear algebra. 
//
/*
 * Class: hyperdual
 * 
 * Implementation of hyper-dual numbers
 *
 * Written by: Jeffrey A. Fike
 * Stanford University, Department of Aeronautics and Astronautics
 *
 * Adapted to Eigen: Michael Tesch
 * 
 * Copyright (c) 2006 Jeffrey A. Fike
 * Copyright (c) 2015 Michael Tesch
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 */
#ifndef EIGEN_HYPERDUAL_H
#define EIGEN_HYPERDUAL_H

namespace Eigen {

#include <math.h>

template <eT>
class hyperdual {
  eT f0, f1, f2, f12;
 public:
  //creation operators and function to manually set values
  hyperdual();
  hyperdual(eT x1, eT x2, eT x3, eT x4);
  hyperdual(eT x1);
  void setvalues(eT x1, eT x2, eT x3, eT x4);

  //examine values
  void view(void);
  eT real(void);
  eT eps1(void);
  eT eps2(void);
  eT eps1eps2(void);
  friend ostream & operator<<(ostream & output, const hyperdual & rhs);

  //basic manipulation
  hyperdual operator+() const;
  hyperdual operator+(const hyperdual rhs) const;
  friend hyperdual operator+(const eT lhs, const hyperdual rhs);
  hyperdual operator-() const;
  hyperdual operator-(const hyperdual rhs) const;
  friend hyperdual operator-(const eT lhs, const hyperdual rhs);
  hyperdual operator*(const hyperdual rhs) const;
  friend hyperdual operator*(const eT lhs, const hyperdual rhs);
  friend hyperdual operator/(const hyperdual lhs, const hyperdual rhs);
  friend hyperdual operator/(const eT lhs, const hyperdual rhs);
  friend hyperdual operator/(const hyperdual lhs, const eT rhs);
  hyperdual & operator+=(hyperdual rhs);
  hyperdual & operator-=(hyperdual rhs);
  hyperdual & operator*=(hyperdual rhs);
  hyperdual & operator*=(eT rhs);
  hyperdual & operator/=(eT rhs);

  //math.h functions
  friend hyperdual pow(hyperdual x, eT a);
  friend hyperdual pow(hyperdual x, hyperdual a);
  friend hyperdual exp(hyperdual x);
  friend hyperdual log(hyperdual x);
  friend hyperdual sin(hyperdual x);
  friend hyperdual cos(hyperdual x);
  friend hyperdual tan(hyperdual x);
  friend hyperdual asin(hyperdual x);
  friend hyperdual acos(hyperdual x);
  friend hyperdual atan(hyperdual x);
  friend hyperdual sqrt(hyperdual x);
  friend hyperdual fabs(hyperdual x);
  friend hyperdual max(hyperdual x1, hyperdual x2);
  friend hyperdual max(hyperdual x1, eT x2);
  friend hyperdual max(eT x1, hyperdual x2);
  friend hyperdual min(hyperdual x1, hyperdual x2);
  friend hyperdual min(hyperdual x1, eT x2);
  friend hyperdual min(eT x1, hyperdual x2);

  //comparisons
  friend bool operator>(hyperdual lhs, hyperdual rhs);
  friend bool operator>(eT lhs, hyperdual rhs);
  friend bool operator>(hyperdual lhs, eT rhs);
  friend bool operator>=(hyperdual lhs, hyperdual rhs);
  friend bool operator>=(eT lhs, hyperdual rhs);
  friend bool operator>=(hyperdual lhs, eT rhs);
  friend bool operator<(hyperdual lhs, hyperdual rhs);
  friend bool operator<(eT lhs, hyperdual rhs);
  friend bool operator<(hyperdual lhs, eT rhs);
  friend bool operator<=(hyperdual lhs, hyperdual rhs);
  friend bool operator<=(eT lhs, hyperdual rhs);
  friend bool operator<=(hyperdual lhs, eT rhs);
  friend bool operator==(hyperdual lhs, hyperdual rhs);
  friend bool operator==(eT lhs, hyperdual rhs);
  friend bool operator==(hyperdual lhs, eT rhs);
  friend bool operator!=(hyperdual lhs, hyperdual rhs);
  friend bool operator!=(eT lhs, hyperdual rhs);
  friend bool operator!=(hyperdual lhs, eT rhs);
};

hyperdual::hyperdual()
{
  f0 = 0.0;
  f1 = 0.0;
  f2 = 0.0;
  f12 = 0.0;
}

hyperdual::hyperdual(eT x1, eT x2, eT x3, eT x4)
{
  f0 = x1;
  f1 = x2;
  f2 = x3;
  f12 = x4;
}

hyperdual::hyperdual(eT x1)
{
  f0 = x1;
  f1 = 0.0;
  f2 = 0.0;
  f12 = 0.0;
}

void hyperdual::setvalues(eT x1, eT x2, eT x3, eT x4)
{
  f0 = x1;
  f1 = x2;
  f2 = x3;
  f12 = x4;
}

void hyperdual::view(void)
{
  printf("%g  +  %g epsilon1  +  %g epsilon2  +  %g epsilon1 epsilon2\n",
         f0, f1, f2, f12);
}

eT hyperdual::real(void)
{
  return f0;
}

eT hyperdual::eps1(void)
{
  return f1;
}

eT hyperdual::eps2(void)
{
  return f2;
}

eT hyperdual::eps1eps2(void)
{
  return f12;
}

ostream & operator<<(ostream & output, const hyperdual & rhs)
{
  output << "(" << rhs.f0 << "," << rhs.f1 << "," << rhs.f2 << "," << rhs.f12 << ")";
  return output;
}

hyperdual hyperdual::operator+() const const
{
  return *this;
}

hyperdual hyperdual::operator+(const hyperdual rhs) constconst
{
  hyperdual temp;
  temp.f0 = f0 + rhs.f0;
  temp.f1 = f1 + rhs.f1;
  temp.f2 = f2 + rhs.f2;
  temp.f12 = f12 + rhs.f12;
  return temp;
}

hyperdual operator+(const eT lhs, const hyperdual rhs)
{
  hyperdual temp;
  temp.f0 = lhs + rhs.f0;
  temp.f1 = rhs.f1;
  temp.f2 = rhs.f2;
  temp.f12 = rhs.f12;
  return temp;
}

hyperdual hyperdual::operator-() const const
{
  hyperdual temp;
  temp.f0 = -f0;
  temp.f1 = -f1;
  temp.f2 = -f2;
  temp.f12 = -f12;
  return temp;
}

hyperdual hyperdual::operator-(const hyperdual rhs) constconst
{
  hyperdual temp;
  temp.f0 = f0 - rhs.f0;
  temp.f1 = f1 - rhs.f1;
  temp.f2 = f2 - rhs.f2;
  temp.f12 = f12 - rhs.f12;
  return temp;
}

hyperdual operator-(const eT lhs, const hyperdual rhs)
{
  hyperdual temp;
  temp.f0 = lhs - rhs.f0;
  temp.f1 = -rhs.f1;
  temp.f2 = -rhs.f2;
  temp.f12 = -rhs.f12;
  return temp;
}

hyperdual hyperdual::operator*(const hyperdual rhs) constconst
{
  hyperdual temp;
  temp.f0 = f0 * rhs.f0;
  temp.f1 = f0 * rhs.f1 + f1 * rhs.f0;
  temp.f2 = f0 * rhs.f2 + f2 * rhs.f0;
  temp.f12 = f0 * rhs.f12 + f1 * rhs.f2 + f2 * rhs.f1 + f12 * rhs.f0;
  return temp;
}

hyperdual operator*(const eT lhs, const hyperdual rhs)
{
  hyperdual temp;
  temp.f0 = lhs * rhs.f0;
  temp.f1 = lhs * rhs.f1;
  temp.f2 = lhs * rhs.f2;
  temp.f12 = lhs * rhs.f12;
  return temp;
}

hyperdual operator/(const hyperdual lhs, const hyperdual rhs)
{
  hyperdual temp, inv;
  inv = pow(rhs, -1);
  temp = lhs * inv;
  return temp;
}

hyperdual operator/(const eT lhs, const hyperdual rhs)
{
  hyperdual temp, inv;
  inv = pow(rhs, -1);
  temp = lhs * inv;
  return temp;
}

hyperdual operator/(const hyperdual lhs, const eT rhs)
{
  hyperdual temp;
  eT inv;
  inv = 1.0 / rhs;
  temp.f0 = inv * lhs.f0;
  temp.f1 = inv * lhs.f1;
  temp.f2 = inv * lhs.f2;
  temp.f12 = inv * lhs.f12;
  return temp;
}

hyperdual & hyperdual::operator+=(hyperdual rhs)
{
  f0 += rhs.f0;
  f1 += rhs.f1;
  f2 += rhs.f2;
  f12 += rhs.f12;
  return *this;
}

hyperdual & hyperdual::operator-=(hyperdual rhs)
{
  f0 -= rhs.f0;
  f1 -= rhs.f1;
  f2 -= rhs.f2;
  f12 -= rhs.f12;
  return *this;
}

hyperdual & hyperdual::operator*=(hyperdual rhs)
{
  eT tf0, tf1, tf2, tf12;
  tf0 = f0;
  tf1 = f1;
  tf2 = f2;
  tf12 = f12;
  f0 = tf0 * rhs.f0;
  f1 = tf0 * rhs.f1 + tf1 * rhs.f0;
  f2 = tf0 * rhs.f2 + tf2 * rhs.f0;
  f12 = tf0 * rhs.f12 + tf1 * rhs.f2 + tf2 * rhs.f1 + tf12 * rhs.f0;
  return *this;
}

hyperdual & hyperdual::operator*=(eT rhs)
{
  f0 *= rhs;
  f1 *= rhs;
  f2 *= rhs;
  f12 *= rhs;
  return *this;
}

hyperdual & hyperdual::operator/=(eT rhs)
{
  f0 /= rhs;
  f1 /= rhs;
  f2 /= rhs;
  f12 /= rhs;
  return *this;
}

hyperdual pow(hyperdual x, eT a)
{
  hyperdual temp;
  eT deriv, xval, tol;
  xval = x.f0;
  tol = 1e-15;
  if (fabs(xval) < tol) {
    if (xval >= 0)
      xval = tol;
    if (xval < 0)
      xval = -tol;
  }
  deriv = a * pow(xval, (a - 1));
  //temp.f0 = pow(xval,a);
  temp.f0 = pow(x.f0, a);       //Use actual x value, only use tol for derivs
  temp.f1 = x.f1 * deriv;
  temp.f2 = x.f2 * deriv;
  temp.f12 = x.f12 * deriv + a * (a - 1) * x.f1 * x.f2 * pow(xval, (a - 2));

  return temp;
}

hyperdual pow(hyperdual x, hyperdual a)
{
  return exp(a * log(x));
}

hyperdual exp(hyperdual x)
{
  hyperdual temp;
  eT deriv;
  deriv = exp(x.f0);
  temp.f0 = deriv;
  temp.f1 = deriv * x.f1;
  temp.f2 = deriv * x.f2;
  temp.f12 = deriv * (x.f12 + x.f1 * x.f2);
  return temp;
}

hyperdual log(hyperdual x)
{
  hyperdual temp;
  eT deriv1, deriv2;
  deriv1 = x.f1 / x.f0;
  deriv2 = x.f2 / x.f0;
  temp.f0 = log(x.f0);
  temp.f1 = deriv1;
  temp.f2 = deriv2;
  temp.f12 = x.f12 / x.f0 - (deriv1 * deriv2);
  return temp;
}

hyperdual sin(hyperdual x)
{
  hyperdual temp;
  eT funval, deriv;
  funval = sin(x.f0);
  deriv = cos(x.f0);
  temp.f0 = funval;
  temp.f1 = deriv * x.f1;
  temp.f2 = deriv * x.f2;
  temp.f12 = deriv * x.f12 - funval * x.f1 * x.f2;
  return temp;
}

hyperdual cos(hyperdual x)
{
  hyperdual temp;
  eT funval, deriv;
  funval = cos(x.f0);
  deriv = -sin(x.f0);
  temp.f0 = funval;
  temp.f1 = deriv * x.f1;
  temp.f2 = deriv * x.f2;
  temp.f12 = deriv * x.f12 - funval * x.f1 * x.f2;
  return temp;
}

hyperdual tan(hyperdual x)
{
  hyperdual temp;
  eT funval, deriv;
  funval = tan(x.f0);
  deriv = funval * funval + 1.0;
  temp.f0 = funval;
  temp.f1 = deriv * x.f1;
  temp.f2 = deriv * x.f2;
  temp.f12 = deriv * x.f12 + x.f1 * x.f2 * (2 * funval * deriv);
  return temp;
}

hyperdual asin(hyperdual x)
{
  hyperdual temp;
  eT funval, deriv1, deriv;
  funval = asin(x.f0);
  deriv1 = 1.0 - x.f0 * x.f0;
  deriv = 1.0 / sqrt(deriv1);
  temp.f0 = funval;
  temp.f1 = deriv * x.f1;
  temp.f2 = deriv * x.f2;
  temp.f12 = deriv * x.f12 + x.f1 * x.f2 * (x.f0 * pow(deriv1, -1.5));
  return temp;
}

hyperdual acos(hyperdual x)
{
  hyperdual temp;
  eT funval, deriv1, deriv;
  funval = acos(x.f0);
  deriv1 = 1.0 - x.f0 * x.f0;
  deriv = -1.0 / sqrt(deriv1);
  temp.f0 = funval;
  temp.f1 = deriv * x.f1;
  temp.f2 = deriv * x.f2;
  temp.f12 = deriv * x.f12 + x.f1 * x.f2 * (-x.f0 * pow(deriv1, -1.5));
  return temp;
}

hyperdual atan(hyperdual x)
{
  hyperdual temp;
  eT funval, deriv1, deriv;
  funval = atan(x.f0);
  deriv1 = 1.0 + x.f0 * x.f0;
  deriv = 1.0 / deriv1;
  temp.f0 = funval;
  temp.f1 = deriv * x.f1;
  temp.f2 = deriv * x.f2;
  temp.f12 = deriv * x.f12 + x.f1 * x.f2 * (-2 * x.f0 / (deriv1 * deriv1));
  return temp;
}

hyperdual sqrt(hyperdual x)
{
  return pow(x, 0.5);
}

hyperdual fabs(hyperdual x)
{
  hyperdual temp;
  if (x < 0.0)
    temp = -x;
  else
    temp = x;
  return temp;
}

hyperdual max(hyperdual x1, hyperdual x2)
{
  hyperdual temp;
  if (x1 > x2)
    temp = x1;
  else
    temp = x2;
  return temp;
}

hyperdual max(hyperdual x1, eT x2)
{
  hyperdual temp;
  if (x1 > x2)
    temp = x1;
  else
    temp = x2;
  return temp;
}

hyperdual max(eT x1, hyperdual x2)
{
  hyperdual temp;
  if (x1 > x2)
    temp = x1;
  else
    temp = x2;
  return temp;
}

hyperdual min(hyperdual x1, hyperdual x2)
{
  hyperdual temp;
  if (x1 < x2)
    temp = x1;
  else
    temp = x2;
  return temp;
}

hyperdual min(hyperdual x1, eT x2)
{
  hyperdual temp;
  if (x1 < x2)
    temp = x1;
  else
    temp = x2;
  return temp;
}

hyperdual min(eT x1, hyperdual x2)
{
  hyperdual temp;
  if (x1 < x2)
    temp = x1;
  else
    temp = x2;
  return temp;
}

bool operator>(hyperdual lhs, hyperdual rhs)
{
  return (lhs.f0 > rhs.f0);
}

bool operator>(eT lhs, hyperdual rhs)
{
  return (lhs > rhs.f0);
}

bool operator>(hyperdual lhs, eT rhs)
{
  return (lhs.f0 > rhs);
}

bool operator>=(hyperdual lhs, hyperdual rhs)
{
  return (lhs.f0 >= rhs.f0);
}

bool operator>=(eT lhs, hyperdual rhs)
{
  return (lhs >= rhs.f0);
}

bool operator>=(hyperdual lhs, eT rhs)
{
  return (lhs.f0 >= rhs);
}

bool operator<(hyperdual lhs, hyperdual rhs)
{
  return (lhs.f0 < rhs.f0);
}

bool operator<(eT lhs, hyperdual rhs)
{
  return (lhs < rhs.f0);
}

bool operator<(hyperdual lhs, eT rhs)
{
  return (lhs.f0 < rhs);
}

bool operator<=(hyperdual lhs, hyperdual rhs)
{
  return (lhs.f0 <= rhs.f0);
}

bool operator<=(eT lhs, hyperdual rhs)
{
  return (lhs <= rhs.f0);
}

bool operator<=(hyperdual lhs, eT rhs)
{
  return (lhs.f0 <= rhs);
}

bool operator==(hyperdual lhs, hyperdual rhs)
{
  return (lhs.f0 == rhs.f0);
}

bool operator==(eT lhs, hyperdual rhs)
{
  return (lhs == rhs.f0);
}

bool operator==(hyperdual lhs, eT rhs)
{
  return (lhs.f0 == rhs);
}

bool operator!=(hyperdual lhs, hyperdual rhs)
{
  return (lhs.f0 != rhs.f0);
}

bool operator!=(eT lhs, hyperdual rhs)
{
  return (lhs != rhs.f0);
}

bool operator!=(hyperdual lhs, eT rhs)
{
  return (lhs.f0 != rhs);
}

template <class eT> using hyperducx = hyperdual<std::complex<eT> >;

#endif
