//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef EXTERNALSIDEINDICATOR_H
#define EXTERNALSIDEINDICATOR_H

// local includes
#include "Indicator.h"
#include "Coupleable.h"
#include "ScalarCoupleable.h"
#include "MooseVariableInterface.h"

// Forward Declarations
class ExternalSideIndicator;
template <typename>
class MooseVariableFE;
typedef MooseVariableFE<Real> MooseVariable;
typedef MooseVariableFE<VectorValue<Real>> VectorMooseVariable;

template <>
InputParameters validParams<ExternalSideIndicator>();

/**
 * The ExternalSideIndicator class is responsible for calculating the residuals for various
 * physics on internal sides (edges/faces).
 *
 */
class ExternalSideIndicator : public Indicator,
                              public Coupleable,
                              public ScalarCoupleable,
                              public MooseVariableInterface<Real>
{
public:
  /**
   * Factory constructor initializes all internal references needed for indicator computation.
   *
   */
  ExternalSideIndicator(const InputParameters & parameters);

  /**
   * Computes the indicator for the current side.
   */
  virtual void computeIndicator() override;

  virtual void finalize() override;

protected:
  MooseVariable & _field_var;

  const Elem *& _current_elem;

  /// Current side
  unsigned int & _current_side;
  /// Current side element
  const Elem *& _current_side_elem;

  /// Coordinate system
  const Moose::CoordinateSystemType & _coord_sys;
  unsigned int _qp;
  const MooseArray<Point> & _q_point;
  QBase *& _qrule;
  const MooseArray<Real> & _JxW;
  const MooseArray<Real> & _coord;

  BoundaryID _boundary_id;

  MooseVariable & _var;

  /// Holds the current solution at the current quadrature point on the face.
  const VariableValue & _u;

  /// Holds the current solution gradient at the current quadrature point on the face.
  const VariableGradient & _grad_u;

  /// Normal vectors at the quadrature points
  const MooseArray<Point> & _normals;

  /// Holds the boundary IDs on which this indicator applies
  std::set<BoundaryID> _sideset_ids;

  /**
   * The virtual function you will want to override to compute error contributions.
   * This is called once per quadrature point on each interior side of every element.
   *
   * You should return the error^2
   */
  virtual Real computeQpIntegral() = 0;

public:
  // boundary id used for internal edges (all DG kernels lives on this boundary id -- a made-up
  // number)
  static const BoundaryID InternalBndId;
};

#endif // EXTERNALSIDEINDICATOR_H
