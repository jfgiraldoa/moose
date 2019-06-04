//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "DivergenceBC.h"

class MatDivergenceBC;

template <>
InputParameters validParams<MatDivergenceBC>();

/**
 * Extends DivergenceBC by multiplication of material property
 */
class MatDivergenceBC : public DivergenceBC
{
public:
  MatDivergenceBC(const InputParameters & parameters);
  virtual ~MatDivergenceBC();

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

  const MaterialProperty<Real> & _mat;
};

