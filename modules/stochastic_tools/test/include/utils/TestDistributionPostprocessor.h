//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralPostprocessor.h"

class TestDistributionPostprocessor;
class Distribution;

template <>
InputParameters validParams<TestDistributionPostprocessor>();

/**
 * Test object for testing distribution capabilities.
 *
 * WARNING! This object is only for testing and should not be used in general.
 */
class TestDistributionPostprocessor : public GeneralPostprocessor
{
public:
  TestDistributionPostprocessor(const InputParameters & parameters);
  virtual void initialize() override {}
  virtual void execute() override {}
  virtual PostprocessorValue getValue() override;

protected:
  /// Object of statistics distribution
  const Distribution & _distribution;

  /// The value to supply to method
  const Real & _value;

  /// The distribution method to call
  const MooseEnum & _distribution_method;
};

