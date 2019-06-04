//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "HEVPStrengthUOBase.h"

class HEVPRambergOsgoodHardening;

template <>
InputParameters validParams<HEVPRambergOsgoodHardening>();

/**
 * This user object classs
 * Computes power law  hardening
 */
class HEVPRambergOsgoodHardening : public HEVPStrengthUOBase
{
public:
  HEVPRambergOsgoodHardening(const InputParameters & parameters);

  virtual bool computeValue(unsigned int, Real &) const;
  virtual bool computeDerivative(unsigned int, const std::string &, Real &) const;

protected:
  Real _sig0;
  Real _peeq0;
  Real _exponent;
};

