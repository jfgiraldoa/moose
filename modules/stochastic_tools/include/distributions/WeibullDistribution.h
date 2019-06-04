//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Distribution.h"

class WeibullDistribution;

template <>
InputParameters validParams<WeibullDistribution>();

/**
 * A class used to generate a three-parameter Weibull distribution
 */
class WeibullDistribution : public Distribution
{
public:
  WeibullDistribution(const InputParameters & parameters);

  virtual Real pdf(const Real & x) const override;
  virtual Real cdf(const Real & x) const override;
  virtual Real quantile(const Real & p) const override;

  Real pdf(const Real & x, const Real & location, const Real & scale, const Real & shape) const;
  Real cdf(const Real & x, const Real & location, const Real & scale, const Real & shape) const;
  Real
  quantile(const Real & p, const Real & location, const Real & scale, const Real & shape) const;

protected:
  /// The location parameter (a or low)
  const Real & _a;

  /// The scale parameter (b or lambda)
  const Real & _b;

  /// The shape parameter (c or k)
  const Real & _c;
};

