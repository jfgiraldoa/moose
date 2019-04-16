//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ExternalSideIndicator.h"

// MOOSE includes
#include "Assembly.h"
#include "MooseTypes.h"
#include "MooseVariableFE.h"
#include "Problem.h"
#include "SubProblem.h"
#include "SystemBase.h"

#include "libmesh/dof_map.h"
#include "libmesh/dense_vector.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/libmesh_common.h"
#include "libmesh/quadrature.h"

template <>
InputParameters
validParams<ExternalSideIndicator>()
{
  InputParameters params = validParams<Indicator>();
  params.addRequiredParam<VariableName>(
      "variable", "The name of the variable that this side indicator applies to");
  params.addRequiredParam<std::vector<BoundaryName>>("sidesets",
    "Names of the sidesets to integrate over");
  return params;
}

ExternalSideIndicator::ExternalSideIndicator(const InputParameters & parameters)
  : Indicator(parameters),
    Coupleable(this, false),
    ScalarCoupleable(this),
    MooseVariableInterface(
        this, false, "variable", Moose::VarKindType::VAR_ANY, Moose::VarFieldType::VAR_FIELD_STANDARD),
    _field_var(_subproblem.getStandardVariable(_tid, name())),

    _current_elem(_assembly.elem()),

    _current_side(_assembly.side()),
    _current_side_elem(_assembly.sideElem()),

    _coord_sys(_assembly.coordSystem()),
    _q_point(_assembly.qPointsFace()),
    _qrule(_assembly.qRuleFace()),
    _JxW(_assembly.JxWFace()),
    _coord(_assembly.coordTransformation()),

    _var(_subproblem.getStandardVariable(_tid, parameters.get<VariableName>("variable"))),

    _u(_var.sln()),
    _grad_u(_var.gradSln()),

    _normals(_field_var.normals())
{
  const std::vector<MooseVariableFEBase *> & coupled_vars = getCoupledMooseVars();
  for (const auto & var : coupled_vars)
    addMooseVariableDependency(var);

  addMooseVariableDependency(mooseVariable());

  if (isParamValid("sidesets"))
  {
    auto & sideset_names = getParam<std::vector<BoundaryName>>("sidesets");
    for (auto & sideset_name : sideset_names)
    {
      _sideset_ids.insert(_mesh.getBoundaryID(sideset_name));
      //std::cout  << "_sideset_ids, name="<<sideset_name<< ", ID="<<_mesh.getBoundaryID(sideset_name)<<std::endl;
    }
  }
  else
  {
    mooseError("'sidesets' param input is not valid.");
  }
}

void
ExternalSideIndicator::computeIndicator()
{
  Real sum = 0;

  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
    sum += _JxW[_qp] * _coord[_qp] * computeQpIntegral();

  {
    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);

    _solution.add(_field_var.nodalDofIndex(), sum * _current_elem->hmax());
  }
}

void
ExternalSideIndicator::finalize()
{
  // The 0 is because CONSTANT MONOMIALS only have one coefficient per element...
  Real value = _field_var.dofValues()[0];

  {
    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    _solution.set(_field_var.nodalDofIndex(), std::sqrt(value));
  }
}
