//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseTypes.h"

/**
 * Rogers-Stallybrass-Clements version of effective saturation as a function of capillary pressure.
 * valid for residual saturations = 0, and viscosityOil = 2*viscosityWater.  (the "2" is important
 * here!).
 * C Rogers, MP Stallybrass and DL Clements "On two phase filtration under gravity and with boundary
 * infiltration: application of a Backlund transformation" Nonlinear Analysis Theory Methods and
 * Applications 7 (1983) 785--799.
 */
namespace PorousFlowRogersStallybrassClements
{
/**
 * Effective saturation as a function of capillary pressure
 * @param pc capillary pressure
 * @param shift RSC's shift parameter
 * @param scale RSC's scale parameter
 * @return effective saturation
 */
Real effectiveSaturation(Real pc, Real shift, Real scale);

/**
 * Derivative of effective saturation wrt capillary pressure
 * @param pc capillary pressure
 * @param shift RSC's shift parameter
 * @param scale RSC's scale parameter
 * @return derivative of effective saturation wrt capillary pressure
 */
Real dEffectiveSaturation(Real pc, Real shift, Real scale);

/**
 * Second derivative of effective saturation wrt capillary pressure
 * @param pc capillary pressure
 * @param shift RSC's shift parameter
 * @param scale RSC's scale parameter
 * @return second derivative of effective saturation wrt capillary pressure
 */
Real d2EffectiveSaturation(Real pc, Real shift, Real scale);
}

