[Mesh]
  type = MeshGeneratorMesh
  displacements = 'disp_x disp_y'
[]

[MeshGenerators]
  [./left]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 2
    ny = 3
    xmin = -3
    xmax = 0
    ymin = -5
    ymax = 5
  [../]
  [./right]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 2
    ny = 3
    xmin = 3
    xmax = 6
    ymin = -5
    ymax = 5
  [../]

  [./left_and_right]
    type = MeshCollectionGenerator
    inputs = 'left right'
  [../]
[]

[MeshModifiers]
  [./leftleft]
    type = SideSetsAroundSubdomain
    block = 0
    new_boundary = leftleft
    normal = '-1 0 0'
  [../]
  [./leftright]
    type = SideSetsAroundSubdomain
    block = 0
    new_boundary = leftright
    normal = '1 0 0'
  [../]

  [./right]
    type = SubdomainBoundingBox
    top_right = '6 5 0'
    bottom_left = '3 -5 0'
    block_id = 1
  [../]

  [./rightleft]
    type = SideSetsAroundSubdomain
    depends_on = right
    block = 1
    new_boundary = rightleft
    normal = '-1 0 0'
  [../]
  [./rightright]
    type = SideSetsAroundSubdomain
    depends_on = right
    block = 1
    new_boundary = rightright
    normal = '1 0 0'
  [../]
[]

[Variables]
  [./temp]
  [../]
[]

[AuxVariables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./gap_conductance]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[Functions]
  [./disp_x]
    type = ParsedFunction
    value = -3+t
  [../]
  [./left_temp]
    type = ParsedFunction
    value = 1000+t
  [../]
[]

[Kernels]
  [./hc]
    type = HeatConduction
    variable = temp
  [../]
[]

[AuxKernels]
  [./disp_x]
    type = FunctionAux
    block = 1
    variable = disp_x
    function = disp_x
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
  [./gap_conductivity]
    type = MaterialRealAux
    boundary = leftright
    property = gap_conductance
    variable = gap_conductance
    execute_on = 'INITIAL TIMESTEP_END'
  [../]
[]

[BCs]
  [./left]
    type = FunctionDirichletBC
    variable = temp
    boundary = leftleft
    function = left_temp
  [../]
  [./right]
    type = DirichletBC
    variable = temp
    boundary = rightright
    value = 400
  [../]
[]

[ThermalContact]
  [./left_to_right]
    slave = leftright
    quadrature = true
    master = rightleft
    variable = temp
    min_gap = 1
    min_gap_order = 1
    type = GapHeatTransfer
  [../]
[]

[Materials]
  [./hcm]
    type = HeatConductionMaterial
    block = '0 1'
    specific_heat = 1
    thermal_conductivity = 1
    use_displaced_mesh = true
  [../]
[]

[Postprocessors]
  [./gap_conductance]
    type = PointValue
    point = '0 0 0'
    variable = gap_conductance
  [../]
[]

[Executioner]
  type = Transient
  dt = 0.25
  end_time = 3.0
  solve_type = 'PJFNK'
[]

[Outputs]
  csv = true
  execute_on = 'TIMESTEP_END'
[]
