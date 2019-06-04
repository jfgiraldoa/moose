#
# Tests material model IsotropicPlasticity with material based time stepper
# Boundary conditions from NAFEMS test NL1
#
[GlobalParams]
  disp_x = disp_x
  disp_y = disp_y
  order = FIRST
  family = LAGRANGE
  volumetric_locking_correction = true
[]

[Mesh]#Comment
  file = one_elem2.e
  displacements = 'disp_x disp_y'
[] # Mesh

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
[] # Variables

[AuxVariables]
  [./stress_xx]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./stress_xy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vonmises]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./elastic_strain_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./plastic_strain_eff]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./tot_strain_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
[] # AuxVariables

[SolidMechanics]
  [./solid]
  [../]
[]

[AuxKernels]
  [./stress_xx]
    type = MaterialTensorAux
    tensor = stress
    variable = stress_xx
    index = 0
  [../]
  [./stress_yy]
    type = MaterialTensorAux
    tensor = stress
    variable = stress_yy
    index = 1
  [../]
  [./stress_zz]
    type = MaterialTensorAux
    tensor = stress
    variable = stress_zz
    index = 2
  [../]
  [./stress_xy]
    type = MaterialTensorAux
    tensor = stress
    variable = stress_xy
    index = 3
  [../]
  [./vonmises]
    type = MaterialTensorAux
    tensor = stress
    variable = vonmises
    quantity = vonmises
    execute_on = timestep_end
  [../]
  [./elastic_strain_yy]
    type = MaterialTensorAux
    tensor = elastic_strain
    variable = elastic_strain_yy
    index = 1
  [../]
  [./plastic_strain_eff]
    type = MaterialRealAux
    property = effective_plastic_strain
    variable = plastic_strain_eff
  [../]
  [./tot_strain_yy]
    type = MaterialTensorAux
    tensor = total_strain
    variable = tot_strain_yy
    index = 1
  [../]
[] # AuxKernels

[Functions]
  [./appl_dispx]
    type = PiecewiseLinear
    x = '0   1.0   2.0   3.0  4.0   5.0  6.0  7.0  8.0'
    y = '0.0 0.25e-4 0.50e-4 0.50e-4 0.50e-4 0.25e-4 0.0 0.0 0.0'
  [../]
  [./appl_dispy]
    type = PiecewiseLinear
    x = '0   1.0   2.0   3.0  4.0   5.0  6.0  7.0  8.0'
    y = '0.0 0.0  0.0 0.25e-4 0.50e-4 0.50e-4 0.50e-4  0.25e-4 0.0 '
  [../]
[]

[BCs]
  [./side_x]
    type = DirichletBC
    variable = disp_x
    boundary = 101
    value = 0.0
  [../]
  [./origin_x]
    type = DirichletBC
    variable = disp_x
    boundary = 103
    value = 0.0
  [../]
  [./bot_y]
    type = DirichletBC
    variable = disp_y
    boundary = 102
    value = 0.0
  [../]
  [./origin_y]
    type = DirichletBC
    variable = disp_y
    boundary = 103
    value = 0.0
  [../]
  [./top_y]
    type = FunctionPresetBC
    variable = disp_y
    boundary = 1
    function = appl_dispy
  [../]
  [./right_x]
    type = FunctionPresetBC
    variable = disp_x
    boundary = 2
    function = appl_dispx
  [../]
[] # BCs

[Materials]
  [./stiffStuff1]
    type = SolidModel
    block = 1
    youngs_modulus = 250e9
    poissons_ratio = 0.25
    constitutive_model = combined
    formulation = NonlinearPlaneStrain
    large_strain = true
  [../]
  [./combined]
    type = CombinedCreepPlasticity
    block = 1
    submodels = 'isoplas'
    absolute_tolerance = 1e-8
  [../]
  [./isoplas]
    type = IsotropicPlasticity
    block = 1
    yield_stress = 5e6
    hardening_constant = 0.0
    relative_tolerance = 1e-20
    absolute_tolerance = 1e-8
    max_inelastic_increment = 0.000001
  [../]
[] # Materials

[Executioner]
  type = Transient
  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-12
  l_tol = 1e-4
  l_max_its = 100
  nl_max_its = 20

  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.1
    time_t = '1.0  2.0  3.0  4.0  5.0  6.0  7.0  8.0'
    time_dt = '0.1 0.1  0.1  0.1  0.1  0.1  0.1  0.1'
    optimal_iterations = 30
    iteration_window = 9
    growth_factor = 2.0
    cutback_factor = 0.5
    timestep_limiting_postprocessor = matl_ts_min
  [../]

  start_time = 0.0
  num_steps = 1000
  end_time = 8.0
[] # Executioner

[Postprocessors]
  [./matl_ts_min]
    type = MaterialTimeStepPostprocessor
  [../]
  [./stress_xx]
    type = ElementAverageValue
    variable = stress_xx
  [../]
  [./stress_yy]
    type = ElementAverageValue
    variable = stress_yy
  [../]
  [./stress_zz]
    type = ElementAverageValue
    variable = stress_zz
  [../]
  [./vonmises]
    type = ElementAverageValue
    variable = vonmises
  [../]
  [./el_strain_yy]
    type = ElementAverageValue
    variable = elastic_strain_yy
  [../]
  [./plas_strain_eff]
    type = ElementAverageValue
    variable = plastic_strain_eff
  [../]
  [./tot_strain_yy]
    type = ElementAverageValue
    variable = tot_strain_yy
  [../]
  [./disp_x1]
    type = NodalVariableValue
    nodeid = 0
    variable = disp_x
  [../]
  [./disp_x4]
    type = NodalVariableValue
    nodeid = 3
    variable = disp_x
  [../]
  [./disp_y1]
    type = NodalVariableValue
    nodeid = 0
    variable = disp_y
  [../]
  [./disp_y4]
    type = NodalVariableValue
    nodeid = 3
    variable = disp_y
  [../]
  [./_dt]
    type = TimestepSize
  [../]
[]

[Outputs]
  exodus = true
  csv = true
  [./console]
    type = Console
    output_linear = true
  [../]
[] # Outputs
