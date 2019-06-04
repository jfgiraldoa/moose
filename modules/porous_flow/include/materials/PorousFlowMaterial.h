//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"
#include "MaterialData.h"
#include "PorousFlowDictator.h"

class PorousFlowMaterial;

template <>
InputParameters validParams<PorousFlowMaterial>();

/**
 * PorousFlowMaterial is the base class for all PorousFlow Materials
 * It allows users to specify that the Material should be a "nodal"
 * Material, in which Material Properties will be evaluated at
 * nodes (using the Variable's nodal values rather than their quadpoint
 * values).  In a derived class's computeQpProperties,
 * _qp must be recognized as a label for a quadpoint (for
 * ordinary Materials) or a node (for nodal Materials).
 *
 * For the nodal Material case, the Material Properties are sized
 * to max(number of nodes, number of quadpoints).  Only "number of nodes"
 * of these will ever be computed and used: the remaining ones
 * (if any) exist just to make sure that the vectors are correctly
 * sized in MOOSE's copying operations (etc).
 *
 * If number of quadpoints < number of nodes (eg for boundary elements)
 * care should be taken to store the required nodal information in
 * the first number_of_quadpoint elements in the std::vector!
 */
class PorousFlowMaterial : public Material
{
public:
  PorousFlowMaterial(const InputParameters & parameters);

protected:
  /// Correctly sizes nodal materials, then initialises using Material::initStatefulProperties
  virtual void initStatefulProperties(unsigned int n_points) override;

  /// Correctly sizes nodal materials, then computes using Material::computeProperties
  virtual void computeProperties() override;

  /**
   * Makes property with name prop_name to be size equal to
   * max(number of nodes, number of quadpoints) in the current element
   */
  void sizeNodalProperty(const std::string & prop_name);

  /**
   * Makes all supplied properties for this material to be size
   * equal to max(number of nodes, number of quadpoints) in the current element
   */
  void sizeAllSuppliedProperties();

  /**
   * Find the nearest quadpoint to the node labelled by nodenum
   * in the current element
   * @param nodenum the node number in the current element
   * @return the nearest quadpoint
   */
  unsigned nearestQP(unsigned nodenum) const;

  /// Whether the derived class holds nodal values
  const bool _nodal_material;

  /// The variable names UserObject for the PorousFlow variables
  const PorousFlowDictator & _dictator;

  /// Names of variables used to declare/get derivatives in the
  /// DerivativeMaterialInterface to ensure consistency
  const VariableName _pressure_variable_name;
  const VariableName _saturation_variable_name;
  const VariableName _temperature_variable_name;
  const VariableName _mass_fraction_variable_name;
};

