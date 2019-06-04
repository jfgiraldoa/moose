//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADKernel.h"

#define usingTemplKernelValueMembers(T) usingTemplKernelMembers(T)
#define usingKernelValueMembers usingTemplKernelValueMembers(Real)
#define usingVectorKernelValueMembers usingTemplKernelMembers(RealVectorValue)

template <typename, ComputeStage>
class ADKernelValueTempl;

template <ComputeStage compute_stage>
using ADKernelValue = ADKernelValueTempl<Real, compute_stage>;
template <ComputeStage compute_stage>
using ADVectorKernelValue = ADKernelValueTempl<RealVectorValue, compute_stage>;

declareADValidParams(ADKernelValue);
declareADValidParams(ADVectorKernelValue);

template <typename T, ComputeStage compute_stage>
class ADKernelValueTempl : public ADKernelTempl<T, compute_stage>
{
public:
  ADKernelValueTempl(const InputParameters & parameters);

  // See KernelBase base for documentation of these overridden methods
  virtual void computeResidual() override;
  virtual void computeJacobian() override;
  virtual void computeADOffDiagJacobian() override;

protected:
  /**
   * Called before forming the residual for an element
   */
  virtual typename OutputTools<typename Moose::ValueType<T, compute_stage>::type>::OutputValue
  precomputeQpResidual() = 0;

  virtual ADResidual computeQpResidual() override final;

  usingTemplKernelMembers(T);
};

