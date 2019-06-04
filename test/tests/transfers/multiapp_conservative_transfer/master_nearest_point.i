[Mesh]
  type = GeneratedMesh
  dim = 2
  xmax = 1
  ymax = 1
  nx = 10
  ny = 10
[]

[MeshModifiers]
  [block1]
    type = SubdomainBoundingBox
    block_id = 1
    bottom_left = '0.5 0 0'
    top_right = '1 1 0'
  []
[]

[Variables]
  [power_density]
  []
[]

[Functions]
  [pwr_func]
    type = ParsedFunction
    value = '1e3*x*(1-x)+5e2'
  []
[]

[Kernels]
  [diff]
    type = Diffusion
    variable = power_density
  []

  [coupledforce]
    type = BodyForce
    variable = power_density
    function = pwr_func
  []
[]

[BCs]
  [left]
    type = DirichletBC
    variable = power_density
    boundary = left
    value = 0
  []
  [right]
    type = DirichletBC
    variable = power_density
    boundary = right
    value = 1e3
  []
[]

[AuxVariables]
  [from_sub]
  []
[]

[VectorPostprocessors]
  [from_nearest_point]
    type = NearestPointIntegralVariablePostprocessor
    variable = power_density
    points = '0 0.5 0 1 0.5 0'
    execute_on = 'transfer nonlinear TIMESTEP_END'
  []

  [to_nearest_point]
    type = NearestPointIntegralVariablePostprocessor
    variable = from_sub
    points = '0 0.5 0 1 0.5 0'
    execute_on = 'transfer nonlinear TIMESTEP_END'
  []
[]

[Executioner]
  type = Steady
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[MultiApps]
  [sub]
    type = FullSolveMultiApp
    input_files = sub_nearest_point.i
    positions = '0 0 0 0.5 0 0'
    execute_on = timestep_end
  []
[]

[Transfers]
  [to_sub]
    type = MultiAppMeshFunctionTransfer
    direction = to_multiapp
    source_variable = power_density
    variable = from_master
    multi_app = sub
    execute_on = timestep_end

    # The following inputs specify what postprocessors should be conserved
    # 1 NearestPointIntegralVariablePostprocessor is specified on the master
    # side with N points, where N is the number of subapps
    # 1 pp is specified on the subapp side
    from_postprocessors_to_be_preserved = 'from_nearest_point'
    to_postprocessors_to_be_preserved = 'from_master_pp'
  []

  [from_sub]
    type = MultiAppMeshFunctionTransfer
    direction = from_multiapp
    source_variable = sink
    variable = from_sub
    multi_app = sub
    execute_on = timestep_end

    # The following inputs specify what postprocessors should be conserved
    # 1 NearestPointIntegralVariablePostprocessor is specified on the master
    # with N points, where N is the number of subapps
    # 1 pp is specified on the subapp side
    to_postprocessors_to_be_preserved = 'to_nearest_point'
    from_postprocessors_to_be_preserved = 'sink'
  []
[]

[Outputs]
  csv = true
  exodus = true
[]
