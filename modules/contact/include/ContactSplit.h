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
#include "Split.h"

class ContactSplit;

template <>
InputParameters validParams<ContactSplit>();

/**
 * Split-based preconditioner for contact problems.
 */
class ContactSplit : public Split
{
public:
  ContactSplit(const InputParameters & params);
  virtual void setup(const std::string & prefix = "-") override;

#if defined(LIBMESH_HAVE_PETSC) && !PETSC_VERSION_LESS_THAN(3, 3, 0)
protected:
  std::vector<std::string> _contact_master;
  std::vector<std::string> _contact_slave;
  std::vector<int> _contact_displaced;
  std::vector<std::string> _uncontact_master;
  std::vector<std::string> _uncontact_slave;
  std::vector<int> _uncontact_displaced;
  bool _include_all_contact_nodes;
#endif // defined(LIBMESH_HAVE_PETSC) && !PETSC_VERSION_LESS_THAN(3,3,0)
};

