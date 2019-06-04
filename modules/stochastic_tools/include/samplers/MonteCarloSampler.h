//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Sampler.h"

class MonteCarloSampler;

template <>
InputParameters validParams<MonteCarloSampler>();
/**
 * A class used to perform Monte Carlo Sampling
 */
class MonteCarloSampler : public Sampler
{
public:
  MonteCarloSampler(const InputParameters & parameters);

protected:
  virtual std::vector<DenseMatrix<Real>> sample() override;

  /// Number of matrices
  const dof_id_type _num_matrices;

  /// Number of monte carlo samples to create for each distribution
  const dof_id_type _num_samples;
};
