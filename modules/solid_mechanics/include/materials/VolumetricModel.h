//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"

class SymmTensor;

class VolumetricModel;

template <>
InputParameters validParams<VolumetricModel>();

class VolumetricModel : public Material
{
public:
  VolumetricModel(const InputParameters & parameters);
  virtual ~VolumetricModel();

  virtual void modifyStrain(const unsigned int qp,
                            const Real scale_factor,
                            SymmTensor & strain_increment,
                            SymmTensor & dstrain_increment_dT) = 0;

  virtual std::vector<std::string> getDependentMaterialProperties() const
  {
    return std::vector<std::string>(1, "");
  }

private:
  using Material::_qp;
};

