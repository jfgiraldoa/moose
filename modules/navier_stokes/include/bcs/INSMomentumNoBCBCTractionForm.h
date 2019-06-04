//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "INSMomentumNoBCBCBase.h"

// Forward Declarations
class INSMomentumNoBCBCTractionForm;

template <>
InputParameters validParams<INSMomentumNoBCBCTractionForm>();

/**
 * This class implements the "No BC" boundary condition based on the
 * "traction" form of the viscous stress tensor.
 */
class INSMomentumNoBCBCTractionForm : public INSMomentumNoBCBCBase
{
public:
  INSMomentumNoBCBCTractionForm(const InputParameters & parameters);

  virtual ~INSMomentumNoBCBCTractionForm() {}

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned jvar);
};

