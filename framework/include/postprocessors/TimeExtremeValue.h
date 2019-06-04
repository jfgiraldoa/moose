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

// Forward Declarations
class TimeExtremeValue;

// Input parameters
template <>
InputParameters validParams<TimeExtremeValue>();

/// A postprocessor for reporting the max/min value of another postprocessor over time
class TimeExtremeValue : public GeneralPostprocessor
{
public:
  /// What type of extreme value we are going to compute
  enum ExtremeType
  {
    MAX,
    MIN,
    ABS_MAX,
    ABS_MIN
  };

  /**
   * Class constructor
   * @param parameters The input parameters
   */
  TimeExtremeValue(const InputParameters & parameters);
  virtual void initialize() override {}
  virtual void execute() override;
  virtual Real getValue() override;

protected:
  const PostprocessorValue & _postprocessor;

  /// The extreme value type ("max", "min", etc.)
  ExtremeType _type;

  /// The extreme value
  Real & _value;
};

