//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "GeneralUserObject.h"
#include "ElementPairLocator.h"
#include "MooseVariableFE.h"

#include "libmesh/vector_value.h"

// Forward Declarations
class PointValueAtXFEMInterface;
class XFEM;
class LineSegmentCutSetUserObject;

template <>
InputParameters validParams<PointValueAtXFEMInterface>();

class PointValueAtXFEMInterface : public GeneralUserObject
{
public:
  PointValueAtXFEMInterface(const InputParameters & parameters);

  virtual ~PointValueAtXFEMInterface() {}

  virtual void initialize() override;
  virtual void execute() override;
  virtual void finalize() override;

  /**
   * get the map that stores the point index and its values at the positive level set side
   */
  std::map<unsigned int, Real> getValueAtPositiveLevelSet() const
  {
    return _values_positive_level_set_side;
  };

  /**
   * get the map that stores the point index and its values at the negative level set side
   */
  std::map<unsigned int, Real> getValueAtNegativeLevelSet() const
  {
    return _values_negative_level_set_side;
  };

  /**
   * get the map that stores the point index and its gradient at the positive level set side
   */
  std::map<unsigned int, RealVectorValue> getGradientAtPositiveLevelSet() const
  {
    return _grad_values_positive_level_set_side;
  };

  /**
   * get the map that stores the point index and its graident at the negative level set side
   */
  std::map<unsigned int, RealVectorValue> getGradientAtNegativeLevelSet() const
  {
    return _grad_values_negative_level_set_side;
  };

  unsigned int numberPoints() const { return _points.size(); };

protected:
  /**
   * Find the element in the element pairs that contains the point in its physical domain.
   * @param p The point in physical space
   * @param positive_level_set True if the physical domain is in positive level set region
   * @return The Elem containing the point or NULL if this processor doesn't contain an element that
   * contains this point.
   */
  const Elem * getElemContainingPoint(const Point & p, bool positive_level_set);

  /// The Mesh we're using
  MooseMesh & _mesh;

  /// The points to evaluate at
  std::vector<Point> _points;

  /// Pointer to PointLocatorBase object
  std::unique_ptr<PointLocatorBase> _pl;

  /// Pointer to the XFEM controller object
  std::shared_ptr<XFEM> _xfem;

  /// Pointer to ElementPairList object
  const ElementPairLocator::ElementPairList * _elem_pairs;

  /// Pointer to LineSegmentCutSetUserObject object
  const LineSegmentCutSetUserObject * _geo_cut;

  /// Pointer to MooseVariableFEBase object
  MooseVariableFEBase * _var;

  /// The variable number of the level set variable we are operating on
  const unsigned int _level_set_var_number;

  /// System reference
  const System & _system;

  /// The subproblem solution vector
  const NumericVector<Number> * _solution;

  /// Mapping from point index and its values at the positive level set side
  std::map<unsigned int, Real> _values_positive_level_set_side;

  /// Mapping from point index and its values at the negative level set side
  std::map<unsigned int, Real> _values_negative_level_set_side;

  /// Mapping from point index and its gradient at the positive level set side
  std::map<unsigned int, RealVectorValue> _grad_values_positive_level_set_side;

  /// Mapping from point index and its gradient at the negative level set side
  std::map<unsigned int, RealVectorValue> _grad_values_negative_level_set_side;
};

