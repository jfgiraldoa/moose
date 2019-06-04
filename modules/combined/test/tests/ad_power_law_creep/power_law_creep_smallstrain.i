# 1x1x1 unit cube with uniform pressure on top face for the case of small strain.
#  This test does not have a solid mechanics analog because there is not an equvialent
#  small strain with rotations strain calculator material in solid mechanics

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 1
  ny = 1
  nz = 1
[]

[Variables]
  [./temp]
    order = FIRST
    family = LAGRANGE
    initial_condition = 1000.0
  [../]
[]

[Modules/TensorMechanics/Master]
  [./all]
    strain = SMALL
    incremental = true
    add_variables = true
    generate_output = 'stress_yy creep_strain_xx creep_strain_yy creep_strain_zz elastic_strain_yy'
    use_automatic_differentiation = true
  [../]
[]

[Functions]
  [./top_pull]
    type = PiecewiseLinear
    x = '0 1'
    y = '1 1'
  [../]
[]

[Kernels]
  [./heat]
    type = ADHeatConduction
    variable = temp
  [../]
  [./heat_ie]
    type = ADHeatConductionTimeDerivative
    variable = temp
  [../]
[]

[BCs]
  [./u_top_pull]
    type = ADPressure
    variable = disp_y
    component = 1
    boundary = top
    constant = -10.0e6
    function = top_pull
  [../]
  [./u_bottom_fix]
    type = ADPresetBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]
  [./u_yz_fix]
    type = ADPresetBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]
  [./u_xy_fix]
    type = ADPresetBC
    variable = disp_z
    boundary = back
    value = 0.0
  [../]
  [./temp_fix]
    type = DirichletBC
    variable = temp
    boundary = 'bottom top'
    value = 1000.0
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 2e11
    poissons_ratio = 0.3
  [../]
  [./radial_return_stress]
    type = ADComputeMultipleInelasticStress
    inelastic_models = 'power_law_creep'
  [../]
  [./power_law_creep]
    type = ADPowerLawCreepStressUpdate
    coefficient = 1.0e-15
    n_exponent = 4
    activation_energy = 3.0e5
    temperature = temp
  [../]

  [./thermal]
    type = HeatConductionMaterial
    specific_heat = 1.0
    thermal_conductivity = 100.
  [../]
  [./density]
    type = ADDensity
    density = 1.0
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient

  solve_type = 'PJFNK'

  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-ksp_gmres_restart'
  petsc_options_value = '101'

  line_search = 'none'

  l_max_its = 20
  nl_max_its = 20
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-6
  l_tol = 1e-5
  start_time = 0.0
  end_time = 1.0
  num_steps = 10
  dt = 0.1
[]

[Outputs]
  exodus = true
[]
