//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Moose.h"
#include "libmesh/libmesh.h"
#include "libmesh/utility.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/compare_types.h"

namespace MathUtils
{

inline Real
round(Real x)
{
  return ::round(x); // use round from math.h
}

inline Real
sign(Real x)
{
  return x >= 0.0 ? 1.0 : -1.0;
}

Real poly1Log(Real x, Real tol, int deriv);
Real poly2Log(Real x, Real tol, int deriv);
Real poly3Log(Real x, Real tol, int order);
Real poly4Log(Real x, Real tol, int order);
Real taylorLog(Real x);

Real pow(Real x, int e);
Real poly(std::vector<Real> c, const Real x, const bool deriv);

inline Real
heavyside(Real x)
{
  return x < 0.0 ? 0.0 : 1.0;
}
inline Real
positivePart(Real x)
{
  return x > 0.0 ? x : 0.0;
}
inline Real
negativePart(Real x)
{
  return x < 0.0 ? x : 0.0;
}

template <typename T,
          typename T2,
          typename T3,
          typename std::enable_if<ScalarTraits<T>::value && ScalarTraits<T2>::value &&
                                      ScalarTraits<T3>::value,
                                  int>::type = 0>
void
addScaled(const T & a, const T2 & b, T3 & result)
{
  result += a * b;
}

template <typename T,
          typename T2,
          typename T3,
          typename std::enable_if<ScalarTraits<T>::value, int>::type = 0>
void
addScaled(const T & scalar, const NumericVector<T2> & numeric_vector, NumericVector<T3> & result)
{
  result.add(scalar, numeric_vector);
}

template <
    typename T,
    typename T2,
    template <typename> class W,
    template <typename> class W2,
    typename std::enable_if<std::is_same<typename W<T>::index_type, unsigned int>::value &&
                                std::is_same<typename W2<T2>::index_type, unsigned int>::value,
                            int>::type = 0>
typename CompareTypes<T, T2>::supertype
dotProduct(const W<T> & a, const W2<T2> & b)
{
  return a * b;
}

template <typename T,
          typename T2,
          template <typename> class W,
          template <typename> class W2,
          typename std::enable_if<std::is_same<typename W<T>::index_type,
                                               std::tuple<unsigned int, unsigned int>>::value &&
                                      std::is_same<typename W2<T2>::index_type,
                                                   std::tuple<unsigned int, unsigned int>>::value,
                                  int>::type = 0>
typename CompareTypes<T, T2>::supertype
dotProduct(const W<T> & a, const W2<T2> & b)
{
  return a.contract(b);
}

template <typename T>
T
poly(std::vector<Real> c, const T x, const bool derivative)
{
  const unsigned int size = c.size();
  if (size == 0)
    return 0.0;

  T value = c[0];
  if (derivative)
  {
    value *= size - 1;
    for (unsigned int i = 1; i < size - 1; i++)
      value = value * x + c[i] * (size - i - 1);
  }
  else
  {
    for (unsigned int i = 1; i < size; i++)
      value = value * x + c[i];
  }

  return value;
}

template <typename T>
T
clamp(const T & x, Real lowerlimit, Real upperlimit)
{
  if (x < lowerlimit)
    return lowerlimit;
  if (x > upperlimit)
    return upperlimit;
  return x;
}

template <typename T>
T
smootherStep(T x, T start, T end, bool derivative = false)
{
  if (end == start)
    return 0.0;
  x = clamp((x - start) / (end - start), 0.0, 1.0);
  if (x == 0.0)
    return 0.0;
  if (derivative)
  {
    if (x == 1.0)
      return 0.0;
    return 30.0 * Utility::pow<2>(x) * (x * (x - 2.0) + 1.0) / (end - start);
  }
  if (x == 1.0)
    return 1.0;
  return Utility::pow<3>(x) * (x * (x * 6.0 - 15.0) + 10.0);
}

} // namespace MathUtils
