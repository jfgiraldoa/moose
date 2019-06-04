//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Action.h"

// Forward Declarations
class CHPFCRFFSplitKernelAction;

template <>
InputParameters validParams<CHPFCRFFSplitKernelAction>();

/**
 * \todo Needs documentation.
 */
class CHPFCRFFSplitKernelAction : public Action
{
public:
  CHPFCRFFSplitKernelAction(const InputParameters & params);

  virtual void act();

private:
  const unsigned int _num_L;
  const std::string _L_name_base;
  const NonlinearVariableName _n_name;
};

