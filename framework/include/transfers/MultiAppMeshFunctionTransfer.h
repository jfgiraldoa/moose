//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MultiAppFieldTransfer.h"

// Forward declarations
class MultiAppMeshFunctionTransfer;

template <>
InputParameters validParams<MultiAppMeshFunctionTransfer>();

/**
 * Transfers a vector of variables. For each individual one,
 * samples the variable's value in the Master domain at the point where
 * the MultiApp is. Copies that value into a postprocessor in the MultiApp.
 * The source and destination vectors (of variables) should be ordered consistently.
 */
class MultiAppMeshFunctionTransfer : public MultiAppFieldTransfer
{
public:
  MultiAppMeshFunctionTransfer(const InputParameters & parameters);

  virtual void execute() override;

protected:
  /// The number of variables to transfer
  unsigned int _var_size;
  bool _error_on_miss;

private:
  /**
   * Performs the transfer for the variable of index i
   */
  void transferVariable(unsigned int i);

  /// To send points to other processors
  std::vector<std::vector<Parallel::Request>> _send_points;
  /// To send values to other processors
  std::vector<std::vector<Parallel::Request>> _send_evals;
  /// To send app ids to other processors
  std::vector<std::vector<Parallel::Request>> _send_ids;
};

