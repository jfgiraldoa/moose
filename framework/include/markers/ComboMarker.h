//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Marker.h"

class ComboMarker;

/**
 * Combines multiple marker fields.  The most conservative wins.
 */
template <>
InputParameters validParams<ComboMarker>();

class ComboMarker : public Marker
{
public:
  ComboMarker(const InputParameters & parameters);

protected:
  virtual MarkerValue computeElementMarker() override;

  std::vector<MarkerName> _names;

  std::vector<const VariableValue *> _markers;
};

