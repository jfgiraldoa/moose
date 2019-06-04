//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "DerivativeMaterialPropertyNameInterface.h"
#include "Material.h"

#include <string>
#include <vector>

class MooseObject;

/**
 * Material properties get fully described using this structure, including their dependent
 * variables and derivation state.
 */
class FunctionMaterialPropertyDescriptor : public DerivativeMaterialPropertyNameInterface
{
public:
  /*
   * The descriptor is constructed with an expression that describes the
   * material property.
   * Examples:
   *   'F'               A material property called 'F' with no declared variable
   *                     dependencies (i.e. vanishing derivatives)
   *   'F(c,phi)'        A material property called 'F' with declared dependence
   *                     on 'c' and 'phi' (uses DerivativeFunctionMaterial rules to
   *                     look up the derivatives)
   *   'a:=D[x(t),t,t]'  The second time derivative of the t-dependent material property 'x'
   *                     which will be referred to as 'a' in the function expression.
   */
  FunctionMaterialPropertyDescriptor(const std::string &, MooseObject *);

  /// default constructor
  FunctionMaterialPropertyDescriptor();

  /// copy constructor
  FunctionMaterialPropertyDescriptor(const FunctionMaterialPropertyDescriptor &);

  /// copy constructor assigning new parent
  FunctionMaterialPropertyDescriptor(const FunctionMaterialPropertyDescriptor &, MooseObject *);

  /// construct a vector of FunctionMaterialPropertyDescriptors from a vector of strings
  static std::vector<FunctionMaterialPropertyDescriptor>
  parseVector(const std::vector<std::string> &, MooseObject *);

  /// get the fparser symbol name
  const std::string & getSymbolName() const { return _fparser_name; };

  /// set the fparser symbol name
  void setSymbolName(const std::string & n) { _fparser_name = n; };

  /// get the property name
  const std::string getPropertyName() const
  {
    return derivativePropertyName(_base_name, _derivative_vars);
  };

  /// get the property reference
  const MaterialProperty<Real> & value() const;

  /// take another derivative
  void addDerivative(const VariableName & var);

  /**
   * Check if a material property depends on a given variable.
   * A dependency is indicated by either directly specifying it, or by requesting a
   * derivative w.r.t. that variable using the D[x,a] syntax
   */
  bool dependsOn(const std::string & var) const;

  /// builds a list of dependent variables (exactly all variabled for which depends on returns true)
  std::vector<VariableName> getDependentVariables();

  // output the internal state of this descriptor for debugging purposes
  void printDebug();

private:
  void parseDerivative(const std::string &);
  void parseDependentVariables(const std::string &);

  /// name used in function expression
  std::string _fparser_name;

  /// function material property base name
  std::string _base_name;

  std::vector<VariableName> _dependent_vars;
  std::vector<VariableName> _derivative_vars;

  /// material property value (this is lazily updated and cached when read through value())
  mutable const MaterialProperty<Real> * _value;

  /// material object that owns this descriptor
  MooseObject * _parent;
};

