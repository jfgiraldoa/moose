//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Kernel.h"

// Forward Declarations
class StressDivergenceTensorsTruss;

template <>
InputParameters validParams<StressDivergenceTensorsTruss>();

class StressDivergenceTensorsTruss : public Kernel
{
public:
  StressDivergenceTensorsTruss(const InputParameters & parameters);

protected:
  virtual void initialSetup() override;
  virtual void computeResidual() override;
  virtual Real computeQpResidual() override { return 0.0; }
  virtual Real computeStiffness(unsigned int i, unsigned int j);
  virtual void computeJacobian() override;
  virtual void computeOffDiagJacobian(MooseVariableFEBase & jvar) override;
  using Kernel::computeOffDiagJacobian;

  std::string _base_name;

  const MaterialProperty<Real> & _axial_stress;
  const MaterialProperty<Real> & _e_over_l;

private:
  const unsigned int _component;
  const unsigned int _ndisp;
  const bool _temp_coupled;

  const unsigned int _temp_var;
  std::vector<unsigned int> _disp_var;

  const VariableValue & _area;
  const std::vector<RealGradient> * _orientation;
};

