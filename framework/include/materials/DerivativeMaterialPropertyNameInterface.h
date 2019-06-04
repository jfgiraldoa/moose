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
#include "MooseError.h"

class DerivativeMaterialPropertyNameInterface
{
public:
  /**
   * Helper functions to generate the material property names for the
   * arbitrary derivatives.
   */
  const MaterialPropertyName derivativePropertyName(const MaterialPropertyName & base,
                                                    const std::vector<VariableName> & c) const;

  /**
   * Helper functions to generate the material property names for the
   * first derivatives.
   */
  const MaterialPropertyName derivativePropertyNameFirst(const MaterialPropertyName & base,
                                                         const VariableName & c1) const;

  /**
   * Helper functions to generate the material property names for the
   * second derivatives.
   */
  const MaterialPropertyName derivativePropertyNameSecond(const MaterialPropertyName & base,
                                                          const VariableName & c1,
                                                          const VariableName & c2) const;

  /**
   * Helper functions to generate the material property names for the
   * third derivatives.
   */
  const MaterialPropertyName derivativePropertyNameThird(const MaterialPropertyName & base,
                                                         const VariableName & c1,
                                                         const VariableName & c2,
                                                         const VariableName & c3) const;

  ///@{ aliases for the deprecated old function names
  const MaterialPropertyName propertyName(const MaterialPropertyName & base,
                                          const std::vector<VariableName> & c) const
  {
    mooseDeprecated("This function was renamed to 'derivativePropertyName'");
    return derivativePropertyName(base, c);
  }
  const MaterialPropertyName propertyNameFirst(const MaterialPropertyName & base,
                                               const VariableName & c1) const
  {
    mooseDeprecated("This function was renamed to 'derivativePropertyNameFirst'");
    return derivativePropertyNameFirst(base, c1);
  }
  const MaterialPropertyName propertyNameSecond(const MaterialPropertyName & base,
                                                const VariableName & c1,
                                                const VariableName & c2) const
  {
    mooseDeprecated("This function was renamed to 'derivativePropertyNameSecond'");
    return derivativePropertyNameSecond(base, c1, c2);
  }
  const MaterialPropertyName propertyNameThird(const MaterialPropertyName & base,
                                               const VariableName & c1,
                                               const VariableName & c2,
                                               const VariableName & c3) const
  {
    mooseDeprecated("This function was renamed to 'derivativePropertyNameThird'");
    return derivativePropertyNameThird(base, c1, c2, c3);
  }
  ///@}
};

