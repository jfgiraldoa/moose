//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElementUserObject.h"

class DiscreteNucleationInserterBase;

template <>
InputParameters validParams<DiscreteNucleationInserterBase>();

/**
 * This UserObject manages the insertion and expiration of nuclei in the simulation
 * domain it manages a list of nuclei with their insertion times and their center
 * positions. A DiscreteNucleationMap is needed to enable the DiscreteNucleation
 * material to look up if a nucleus is present at a given element/qp.
 */
class DiscreteNucleationInserterBase : public ElementUserObject
{
public:
  DiscreteNucleationInserterBase(const InputParameters & parameters);

  /// A nucleus has an expiration time and a location
  using NucleusLocation = std::pair<Real, Point>;

  /// Every MPI task should keep a full list of nuclei (in case they cross domains with their finite radii)
  using NucleusList = std::vector<NucleusLocation>;

  // counter pair to track insertions and deletion in the current timestep
  using NucleusChanges = std::pair<unsigned int, unsigned int>;

  virtual bool isMapUpdateRequired() const { return _update_required; }
  virtual const NucleusList & getNucleusList() const { return _global_nucleus_list; }
  virtual const NucleusChanges & getInsertionsAndDeletions() const { return _changes_made; }

  virtual const Real & getRate() const = 0;

protected:
  /// the global list of all nuclei over all processors
  NucleusList & _global_nucleus_list;

  /// count the number of nucleus insertions and deletions
  NucleusChanges _changes_made;

  /// is a map update required
  bool _update_required;
};

