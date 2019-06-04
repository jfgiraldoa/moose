[Mesh]
  file = 2blk-gap.e
[]

[MeshModifiers]
  [slave]
    type = LowerDBlockFromSideset
    sidesets = '101'
    new_block_id = '10001'
    new_block_name = 'slave_lower'
  []
  [master]
    type = LowerDBlockFromSideset
    sidesets = '100'
    new_block_id = '10000'
    new_block_name = 'master_lower'
  []
[]

[Problem]
  kernel_coverage_check = false
  material_coverage_check = false
[]

[Variables]
  [./temp]
    order = FIRST
    family = LAGRANGE
    block = '1 2'
  [../]

  [./lm]
    order = FIRST
    family = LAGRANGE
    block = 'slave_lower'
  [../]
[]

[Materials]
  [./left]
    type = HeatConductionMaterial
    block = 1
    thermal_conductivity = 1000
    specific_heat = 1
  [../]

  [./right]
    type = HeatConductionMaterial
    block = 2
    thermal_conductivity = 500
    specific_heat = 1
  [../]
[]

[Kernels]
  [./hc]
    type = HeatConduction
    variable = temp
    use_displaced_mesh = false
    block = '1 2'
  [../]
[]

[Constraints]
  [./ced]
    type = GapConductanceConstraint
    variable = lm
    slave_variable = temp
    k = 100
    master_boundary = 100
    master_subdomain = 10000
    slave_boundary = 101
    slave_subdomain = 10001
  [../]
[]

[BCs]
  [./left]
    type = DirichletBC
    variable = temp
    boundary = 'left'
    value = 1
  [../]

  [./right]
    type = DirichletBC
    variable = temp
    boundary = 'right'
    value = 0
  [../]
[]

[Preconditioning]
  [./fmp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Steady
  solve_type = 'PJFNK'
  nl_rel_tol = 1e-11
[]

[Outputs]
  exodus = true
[]
