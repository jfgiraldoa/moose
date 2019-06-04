#
# Test the conserved action with split solve and 1 variable
#
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 50
  ny = 50
  xmax = 50
  ymax = 50
  elem_type = QUAD
[]

[Modules]
  [./PhaseField]
    [./Conserved]
      [./cv]
        solve_type = REVERSE_SPLIT
        free_energy = F
        kappa = 2.0
        mobility = 1.0
      [../]
    [../]
  [../]
[]

[ICs]
  [./InitialCondition]
    type = CrossIC
    x1 = 5.0
    y1 = 5.0
    x2 = 45.0
    y2 = 45.0
    variable = cv
  [../]
[]

[Materials]
  [./free_energy]
    type = DerivativeParsedMaterial
    f_name = F
    args = 'cv'
    function = '(1-cv)^2 * (1+cv)^2'
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  scheme = 'bdf2'

  solve_type = 'NEWTON'

  l_max_its = 30
  l_tol = 1.0e-5
  nl_max_its = 10
  nl_rel_tol = 1.0e-12

  start_time = 0.0
  num_steps = 5
  dt = 0.7
[]

[Outputs]
  exodus = true
[]
