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

class RandomHitUserObject;

class RandomHitMarker;
template <>
InputParameters validParams<RandomHitMarker>();

class RandomHitMarker : public Marker
{
public:
  RandomHitMarker(const InputParameters & parameters);
  virtual ~RandomHitMarker(){};

protected:
  virtual MarkerValue computeElementMarker();

  const RandomHitUserObject & _random_hits;
};

