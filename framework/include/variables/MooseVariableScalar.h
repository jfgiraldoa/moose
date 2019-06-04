//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseVariableBase.h"
#include "SystemBase.h"

// libMesh forward declarations
namespace libMesh
{
template <typename T>
class NumericVector;
}

class Assembly;

/**
 * Class for scalar variables (they are different).
 */
class MooseVariableScalar : public MooseVariableBase
{
public:
  MooseVariableScalar(unsigned int var_num,
                      const FEType & fe_type,
                      SystemBase & sys,
                      Assembly & assembly,
                      Moose::VarKindType var_kind,
                      THREAD_ID tid);
  virtual ~MooseVariableScalar();

  void reinit();

  //
  VariableValue & sln() { return _u; }
  VariableValue & slnOld() { return _u_old; }
  VariableValue & slnOlder() { return _u_older; }
  VariableValue & vectorTagSln(TagID tag)
  {
    _need_vector_tag_u[tag] = true;
    return _vector_tag_u[tag];
  }
  VariableValue & matrixTagSln(TagID tag)
  {
    _need_matrix_tag_u[tag] = true;
    return _matrix_tag_u[tag];
  }

  VariableValue & uDot()
  {
    if (_sys.solutionUDot())
    {
      _need_u_dot = true;
      return _u_dot;
    }
    else
      mooseError("MooseVariableScalar: Time derivative of solution (`u_dot`) is not stored. Please "
                 "set uDotRequested() to true in FEProblemBase before requesting `u_dot`.");
  }

  VariableValue & uDotDot()
  {
    if (_sys.solutionUDotDot())
    {
      _need_u_dotdot = true;
      return _u_dotdot;
    }
    else
      mooseError("MooseVariableScalar: Second time derivative of solution (`u_dotdot`) is not "
                 "stored. Please set uDotDotRequested() to true in FEProblemBase before requesting "
                 "`u_dotdot`.");
  }

  VariableValue & uDotOld()
  {
    if (_sys.solutionUDotOld())
    {
      _need_u_dot_old = true;
      return _u_dot_old;
    }
    else
      mooseError("MooseVariableScalar: Old time derivative of solution (`u_dot_old`) is not "
                 "stored. Please set uDotOldRequested() to true in FEProblemBase before requesting "
                 "`u_dot_old`.");
  }

  VariableValue & uDotDotOld()
  {
    if (_sys.solutionUDotDotOld())
    {
      _need_u_dotdot_old = true;
      return _u_dotdot_old;
    }
    else
      mooseError("MooseVariableScalar: Old second time derivative of solution (`u_dotdot_old`) is "
                 "not stored. Please set uDotDotOldRequested() to true in FEProblemBase before "
                 "requesting `u_dotdot_old`.");
  }

  VariableValue & duDotDu()
  {
    _need_du_dot_du = true;
    return _du_dot_du;
  }

  VariableValue & duDotDotDu()
  {
    _need_du_dotdot_du = true;
    return _du_dotdot_du;
  }

  /**
   * Set the nodal value for this variable (to keep everything up to date
   */
  void setValue(unsigned int i, Number value);

  /**
   * Set all of the values of this scalar variable to the same value
   */
  void setValues(Number value);

  void insert(NumericVector<Number> & soln);

protected:
  /// The assembly
  Assembly & _assembly;

  /// The value of scalar variable
  VariableValue _u;
  /// The old value of scalar variable
  VariableValue _u_old;
  /// The older value of scalar variable
  VariableValue _u_older;
  /// Tagged vectors
  std::vector<VariableValue> _vector_tag_u;
  /// Only cache data when need it
  std::vector<bool> _need_vector_tag_u;
  /// Tagged matrices
  std::vector<VariableValue> _matrix_tag_u;
  /// Only cache data when need it
  std::vector<bool> _need_matrix_tag_u;

  VariableValue _u_dot;
  VariableValue _u_dotdot;
  VariableValue _u_dot_old;
  VariableValue _u_dotdot_old;
  VariableValue _du_dot_du;
  VariableValue _du_dotdot_du;

  bool _need_u_dot;
  bool _need_u_dotdot;
  bool _need_u_dot_old;
  bool _need_u_dotdot_old;
  bool _need_du_dot_du;
  bool _need_du_dotdot_du;
};

