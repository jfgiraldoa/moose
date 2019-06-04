[Mesh]
  type = GeneratedMesh
  dim = 2
  xmin = 0.01 # to make sure the meshes don't align
  xmax = 0.49 # to make sure the meshes don't align
  ymax = 1
  nx = 10
  ny = 10
[]

[MeshModifiers]
  [block1]
    type = SubdomainBoundingBox
    block_id = 1
    bottom_left = '0.2 0.2 0'
    top_right = '0.3 0.8 0'
  []
[]

[Variables]
  [sink]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[Functions]
  [sink_func]
    type = ParsedFunction
    value = '5e2*x*(0.5-x)+5e1'
  []
[]

[Kernels]
  [reaction]
    type = Reaction
    variable = sink
  []

  [coupledforce]
    type = BodyForce
    variable = sink
    function = sink_func
  []
[]

[AuxVariables]
  [from_master]
    block = 1
  []
[]

[Executioner]
  type = Steady
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Postprocessors]
  [sink]
    type = ElementIntegralVariablePostprocessor
    block = 1
    variable = sink
    execute_on = 'transfer nonlinear TIMESTEP_END'
  []
  [from_master_pp]
    type = ElementIntegralVariablePostprocessor
    block = 1
    variable = from_master
    execute_on = 'transfer nonlinear TIMESTEP_END'
  []
[]

[Outputs]
  exodus = true
[]
