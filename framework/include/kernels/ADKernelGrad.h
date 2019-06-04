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

#define usingTemplKernelGradMembers(T) usingTemplKernelMembers(T)
#define usingKernelGradMembers usingTemplKernelGradMembers(Real)
#define usingVectorKernelGradMembers usingTemplKernelMembers(RealVectorValue)

template <typename, ComputeStage>
class ADKernelGradTempl;

template <ComputeStage compute_stage>
using ADKernelGrad = ADKernelGradTempl<Real, compute_stage>;
template <ComputeStage compute_stage>
using ADVectorKernelGrad = ADKernelGradTempl<RealVectorValue, compute_stage>;

declareADValidParams(ADKernelGrad);
declareADValidParams(ADVectorKernelGrad);

template <typename T, ComputeStage compute_stage>
class ADKernelGradTempl : public ADKernelTempl<T, compute_stage>
{
public:
  ADKernelGradTempl(const InputParameters & parameters);

  virtual void computeResidual() override;
  virtual void computeJacobian() override;
  virtual void computeADOffDiagJacobian() override;

protected:
  /**
   * Called before forming the residual for an element
   */
  virtual typename OutputTools<typename Moose::ValueType<T, compute_stage>::type>::OutputGradient
  precomputeQpResidual() = 0;

  virtual ADResidual computeQpResidual() override final;

  usingTemplKernelMembers(T);
};

