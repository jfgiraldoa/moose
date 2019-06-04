//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "TimeNodalKernel.h"

// Forward Declarations
class NodalTranslationalInertia;

template <>
InputParameters validParams<NodalTranslationalInertia>();

/**
 * Calculates the inertial force and mass proportional damping for a nodal mass
 */
class NodalTranslationalInertia : public TimeNodalKernel
{
public:
  NodalTranslationalInertia(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;

  /// Booleans for validity of params
  const bool _has_mass;
  const bool _has_beta;
  const bool _has_gamma;
  const bool _has_velocity;
  const bool _has_acceleration;
  const bool _has_nodal_mass_file;

  /// Mass associated with the node
  const Real _mass;

  /// Old value of displacement
  const VariableValue * _u_old;

  /// Newmark time integration parameter
  const Real _beta;

  /// Newmark time integration parameter
  const Real _gamma;

  /// Mass proportional Rayliegh damping
  const Real & _eta;

  /// HHT time integration parameter
  const Real & _alpha;

  /// Auxiliary system object
  AuxiliarySystem * _aux_sys;

  /// Variable number corresponding to the velocity aux variable
  unsigned int _vel_num;

  /// Variable number corresponding to the acceleration aux variable
  unsigned int _accel_num;

  /// Map between boundary nodes and nodal mass
  std::map<dof_id_type, Real> _node_id_to_mass;

  /// Velocity variable value
  const MooseArray<Number> * _vel;

  /// Old velocity variable value
  const MooseArray<Number> * _vel_old;

  /// Acceleration variable value
  const MooseArray<Number> * _accel;

  /// du_dot_du variable value
  const MooseArray<Number> * _du_dot_du;

  /// du_dotdot_du variable value
  const MooseArray<Number> * _du_dotdot_du;
};

