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
class NumDOFs;

// libMesh forward declarations
namespace libMesh
{
class System;
class EquationSystems;
}

template <>
InputParameters validParams<NumDOFs>();

class NumDOFs : public GeneralPostprocessor
{
public:
  NumDOFs(const InputParameters & parameters);

  virtual void initialize() override {}
  virtual void execute() override {}
  virtual Real getValue() override;

protected:
  enum SystemEnum
  {
    NL,
    AUX,
    ALL
  };

  const SystemEnum _system_enum;

  const System * _system_pointer;
  const EquationSystems * _es_pointer;
};

