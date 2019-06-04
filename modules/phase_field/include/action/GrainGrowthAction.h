//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "Action.h"

#include "libmesh/fe_type.h"

// Forward declaration
class GrainGrowthAction;

template <>
InputParameters validParams<GrainGrowthAction>();

class GrainGrowthAction : public Action
{
public:
  GrainGrowthAction(const InputParameters & params);

  virtual void act();

protected:
  /// number of variables and variable name base for variable creation
  const unsigned int _op_num;
  const std::string _var_name_base;

  /// FEType for the variable being created
  const FEType _fe_type;
};

