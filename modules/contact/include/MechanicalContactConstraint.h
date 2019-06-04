//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "NodeFaceConstraint.h"
#include "ContactMaster.h"

// Forward Declarations
class MechanicalContactConstraint;
class ContactLineSearchBase;

template <>
InputParameters validParams<MechanicalContactConstraint>();

/**
 * A MechanicalContactConstraint forces the value of a variable to be the same on both sides of an
 * interface.
 */
class MechanicalContactConstraint : public NodeFaceConstraint
{
public:
  MechanicalContactConstraint(const InputParameters & parameters);

  virtual void timestepSetup() override;
  virtual void jacobianSetup() override;
  virtual void residualEnd() override;

  virtual bool AugmentedLagrangianContactConverged();

  virtual void updateAugmentedLagrangianMultiplier(bool beginning_of_step = false);

  virtual void updateContactStatefulData(bool beginning_of_step = false);

  virtual Real computeQpSlaveValue() override;

  virtual Real computeQpResidual(Moose::ConstraintType type) override;

  /**
   * Computes the jacobian for the current element.
   */
  virtual void computeJacobian() override;

  /**
   * Compute off-diagonal Jacobian entries
   * @param jvar The index of the coupled variable
   */
  virtual void computeOffDiagJacobian(unsigned int jvar) override;

  virtual Real computeQpJacobian(Moose::ConstraintJacobianType type) override;

  /**
   * Compute off-diagonal Jacobian entries
   * @param type The type of coupling
   * @param jvar The index of the coupled variable
   */
  virtual Real computeQpOffDiagJacobian(Moose::ConstraintJacobianType type,
                                        unsigned int jvar) override;

  /**
   * Get the dof indices of the nodes connected to the slave node for a specific variable
   * @param var_num The number of the variable for which dof indices are gathered
   * @return bool indicating whether the coupled variable is one of the displacement variables
   */
  virtual void getConnectedDofIndices(unsigned int var_num) override;

  /**
   * Determine whether the coupled variable is one of the displacement variables,
   * and find its component
   * @param var_num The number of the variable to be checked
   * @param component The component index computed in this routine
   * @return bool indicating whether the coupled variable is one of the displacement variables
   */
  bool getCoupledVarComponent(unsigned int var_num, unsigned int & component);

  virtual bool addCouplingEntriesToJacobian() override { return _master_slave_jacobian; }

  bool shouldApply() override;
  void computeContactForce(PenetrationInfo * pinfo, bool update_contact_set);

protected:
  MooseSharedPointer<DisplacedProblem> _displaced_problem;
  FEProblem & _fe_problem;
  Real nodalArea(PenetrationInfo & pinfo);
  Real getPenalty(PenetrationInfo & pinfo);
  Real getTangentialPenalty(PenetrationInfo & pinfo);

  const unsigned int _component;
  ContactModel _model;
  const ContactFormulation _formulation;
  const bool _normalize_penalty;

  const Real _penalty;
  Real _penalty_tangential;
  const Real _friction_coefficient;
  const Real _tension_release;
  const Real _capture_tolerance;
  const unsigned int _stick_lock_iterations;
  const Real _stick_unlock_factor;
  bool _update_stateful_data;

  NumericVector<Number> & _residual_copy;
  //  std::map<Point, PenetrationInfo *> _point_to_info;

  const unsigned int _mesh_dimension;

  std::vector<unsigned int> _vars;
  std::vector<MooseVariable *> _var_objects;

  MooseVariable * _nodal_area_var;
  SystemBase & _aux_system;
  const NumericVector<Number> * _aux_solution;

  /// Whether to include coupling between the master and slave nodes in the Jacobian
  const bool _master_slave_jacobian;
  /// Whether to include coupling terms with the nodes connected to the slave nodes in the Jacobian
  const bool _connected_slave_nodes_jacobian;
  /// Whether to include coupling terms with non-displacement variables in the Jacobian
  const bool _non_displacement_vars_jacobian;

  /// The tolerance of the penetration for augmented Lagrangian method
  Real _al_penetration_tolerance;
  /// The tolerance of the incremental slip for augmented Lagrangian method
  Real _al_incremental_slip_tolerance;
  /// The tolerance of the frictional force for augmented Lagrangian method
  Real _al_frictional_force_tolerance;

  std::shared_ptr<ContactLineSearchBase> _contact_linesearch;
  std::set<dof_id_type> _current_contact_state;
  std::set<dof_id_type> _old_contact_state;

  const bool _print_contact_nodes;
  static Threads::spin_mutex _contact_set_mutex;
};

