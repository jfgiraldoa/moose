# This is a basic test of the system for continuum damage mechanics
# materials. It uses ScalarMaterialDamage for the damage model,
# which simply gets its damage index from another material. In this
# case, we prescribe the evolution of the damage index using a
# function. A single element has a fixed prescribed displacement
# on one side that puts the element in tension, and then the
# damage index evolves from 0 to 1 over time, and this verifies
# that the stress correspondingly drops to 0.

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 1
  ny = 1
  nz = 1
  elem_type = HEX8
[]

[AuxVariables]
  [damage_index]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Modules/TensorMechanics/Master]
  [all]
    strain = SMALL
    incremental = true
    add_variables = true
    generate_output = 'stress_xx strain_xx'
  []
[]

[AuxKernels]
  [damage_index]
    type = MaterialRealAux
    variable = damage_index
    property = damage_index_prop
    execute_on = timestep_end
  []
[]

[BCs]
  [symmy]
    type = PresetBC
    variable = disp_y
    boundary = bottom
    value = 0
  []
  [symmx]
    type = PresetBC
    variable = disp_x
    boundary = left
    value = 0
  []
  [symmz]
    type = PresetBC
    variable = disp_z
    boundary = back
    value = 0
  []
  [axial_load]
    type = PresetBC
    variable = disp_x
    boundary = right
    value = 0.01
  []
[]

[Functions]
  [damage_evolution]
    type = PiecewiseLinear
    xy_data = '0.0   0.0
               0.1   0.0
               2.1   2.0'
  []
[]

[Materials]
  [damage_index]
    type = GenericFunctionMaterial
    prop_names = damage_index_prop
    prop_values = damage_evolution
  []
  [damage]
    type = ScalarMaterialDamage
    damage_index = damage_index_prop
  []
  [stress]
    type = ComputeDamageStress
    damage_model = damage
  []
  [elasticity]
    type = ComputeIsotropicElasticityTensor
    poissons_ratio = 0.2
    youngs_modulus = 10e9
  []
[]

[Postprocessors]
  [stress_xx]
    type = ElementAverageValue
    variable = stress_xx
  []
  [strain_xx]
    type = ElementAverageValue
    variable = strain_xx
  []
  [damage_index]
    type = ElementAverageValue
    variable = damage_index
  []
  [time_step_limit]
    type = MaterialTimeStepPostprocessor
  []
[]

[Executioner]
  type = Transient

  l_max_its  = 50
  l_tol      = 1e-8
  nl_max_its = 20
  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-8

  dt = 0.1
  dtmin = 0.001
  end_time = 1.1
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.1
    growth_factor = 2.0
    cutback_factor = 0.5
    timestep_limiting_postprocessor = time_step_limit
  []
[]

[Outputs]
  csv=true
[]
