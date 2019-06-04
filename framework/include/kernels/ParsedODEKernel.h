//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ODEKernel.h"
#include "FunctionParserUtils.h"

// Forward Declarations
class ParsedODEKernel;

template <>
InputParameters validParams<ParsedODEKernel>();

/**
 *
 */
class ParsedODEKernel : public ODEKernel, public FunctionParserUtils
{
public:
  ParsedODEKernel(const InputParameters & parameters);

protected:
  void updateParams();

  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// function expression
  std::string _function;

  /// coupled variables
  unsigned int _nargs;
  std::vector<VariableValue *> _args;
  std::vector<std::string> _arg_names;

  /// function parser object for the residual and on-diagonal Jacobian
  ADFunctionPtr _func_F;
  ADFunctionPtr _func_dFdu;

  /// function parser objects for the Jacobian
  std::vector<ADFunctionPtr> _func_dFdarg;

  /// number of non-linear variables in the problem
  const unsigned int _number_of_nl_variables;

private:
  /// Vector to look up the internal coupled variable index into _arg_*  through the libMesh variable number
  std::vector<unsigned int> _arg_index;
};

