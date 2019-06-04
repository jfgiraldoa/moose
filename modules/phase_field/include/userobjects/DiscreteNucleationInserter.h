//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "DiscreteNucleationInserterBase.h"

class DiscreteNucleationInserter;

template <>
InputParameters validParams<DiscreteNucleationInserter>();

/**
 * This UserObject manages the insertion and expiration of nuclei in the simulation
 * domain it manages a list of nuclei with their insertion times and their center
 * positions. A DiscreteNucleationMap is needed to enable the DiscreteNucleation
 * material to look up if a nucleus is present at a given element/qp.
 */
class DiscreteNucleationInserter : public DiscreteNucleationInserterBase
{
public:
  DiscreteNucleationInserter(const InputParameters & parameters);

  virtual void initialize();
  virtual void execute();
  virtual void threadJoin(const UserObject & y);
  virtual void finalize();

  const Real & getRate() const { return _nucleation_rate; }

protected:
  /// Nucleation rate density (should be a material property implementing nucleation theory)
  const MaterialProperty<Real> & _probability;

  /// Duration of time each nucleus is kept active after insertion
  Real _hold_time;

  /// the local nucleus list of nuclei centered in the domain of the current processor
  NucleusList & _local_nucleus_list;

  /// total nucleation rate
  Real _nucleation_rate;
};

