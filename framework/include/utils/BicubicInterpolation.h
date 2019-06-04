//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseTypes.h"

/**
 * This class interpolates tabulated data with a bicubic function. In order to
 * minimize the computational expense of each sample, the coefficients at each
 * point in the tabulated data are computed once in advance, and then accessed
 * during the interpolation.
 *
 * As a result, this bicubic interpolation can be much faster than the bicubic
 * spline interpolation method in BicubicSplineInterpolation.
 *
 * Adapted from Numerical Recipes in C (section 3.6). The terminology used is
 * consistent with that used in Numerical Recipes, where moving over a column
 * corresponds to moving over the x1 coord. Likewise, moving over a row means
 * moving over the x2 coord.
 */
class BicubicInterpolation
{
public:
  BicubicInterpolation(const std::vector<Real> & x1,
                       const std::vector<Real> & x2,
                       const std::vector<std::vector<Real>> & y);

  virtual ~BicubicInterpolation() = default;

  /**
   * Sanity checks on input data
   */
  void errorCheck();

  /**
   * Samples value at point (x1, x2)
   */
  Real sample(Real x1, Real x2);

  /**
   * Samples value and first derivatives at point (x1, x2)
   * Use this function for speed when computing both value and derivatives,
   * as it minimizes the amount of time spent locating the point in the
   * tabulated data
   */
  void sampleValueAndDerivatives(Real x1, Real x2, Real & y, Real & dy1, Real & dy2);

  /**
   * Samples first derivative at point (x1, x2)
   */
  Real sampleDerivative(Real x1, Real x2, unsigned int deriv_var);

  /**
   * Samples second derivative at point (x1, x2)
   */
  Real sample2ndDerivative(Real x1, Real x2, unsigned int deriv_var);

  /**
   * Precompute all of the coefficients for the bicubic interpolation to avoid
   * calculating them repeatedly
   */
  void precomputeCoefficients();

protected:
  /**
   * Find the indices of the dependent values axis which bracket the point xi
   */
  void findInterval(
      const std::vector<Real> & x, Real xi, unsigned int & klo, unsigned int & khi, Real & xs);

  /**
   * Provides the values of the first derivatives in each direction at all
   * points in the table of data, as well as the mixed second derivative.
   * This is implemented using finite differencing, but could be supplied through
   * other means (such as by sampling with cubic splines)
   */
  void tableDerivatives(std::vector<std::vector<Real>> & dy_dx1,
                        std::vector<std::vector<Real>> & dy_dx2,
                        std::vector<std::vector<Real>> & d2y_dx1x2);

  /// Independent values in the x1 direction
  std::vector<Real> _x1;
  /// Independent values in the x2 direction
  std::vector<Real> _x2;
  /// The dependent values at (x1, x2) points
  std::vector<std::vector<Real>> _y;

  /// Matrix of precomputed coefficients
  /// There are four coefficients in each direction at each dependent value
  std::vector<std::vector<std::vector<std::vector<Real>>>> _bicubic_coeffs;

  /// Matrix used to calculate bicubic interpolation coefficients
  /// (from Numerical Recipes)
  const std::vector<std::vector<int>> _wt{
      {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
      {-3, 0, 0, 3, 0, 0, 0, 0, -2, 0, 0, -1, 0, 0, 0, 0},
      {2, 0, 0, -2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0},
      {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
      {0, 0, 0, 0, -3, 0, 0, 3, 0, 0, 0, 0, -2, 0, 0, -1},
      {0, 0, 0, 0, 2, 0, 0, -2, 0, 0, 0, 0, 1, 0, 0, 1},
      {-3, 3, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, -2, -1, 0, 0},
      {9, -9, 9, -9, 6, 3, -3, -6, 6, -6, -3, 3, 4, 2, 1, 2},
      {-6, 6, -6, 6, -4, -2, 2, 4, -3, 3, 3, -3, -2, -1, -1, -2},
      {2, -2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 1, 1, 0, 0},
      {-6, 6, -6, 6, -3, -3, 3, 3, -4, 4, 2, -2, -2, -2, -1, -1},
      {4, -4, 4, -4, 2, 2, -2, -2, 2, -2, -2, 2, 1, 1, 1, 1}};
};

