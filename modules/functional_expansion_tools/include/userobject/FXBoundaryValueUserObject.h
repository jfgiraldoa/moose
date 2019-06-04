//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "FXBoundaryBaseUserObject.h"

class FXBoundaryValueUserObject;

template <>
InputParameters validParams<FXBoundaryValueUserObject>();

/**
 * This boundary FX evaluator calculates the values
 */
class FXBoundaryValueUserObject final : public FXBoundaryBaseUserObject
{
public:
  FXBoundaryValueUserObject(const InputParameters & parameters);
};

