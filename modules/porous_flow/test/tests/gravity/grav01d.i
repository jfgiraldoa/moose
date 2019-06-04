# Test illustrating that PorousFlow allows block-restricted relative permeabilities and capillarities
# and automatically adds appropriate Joiners.
# Physically, this test is checking that gravity head is established
# for 1phase, vanGenuchten, constant fluid-bulk, constant viscosity, constant permeability, Corey relative perm
# For better agreement with the analytical solution (ana_pp), just increase nx

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 100
  xmin = -1
  xmax = 0
[]

[MeshModifiers]
  [./define_block1]
    type = SubdomainBoundingBox
    block_id = 1
    bottom_left = '-1 -1 -1'
    top_right = '-0.5 1 1'
  [../]
[]

[GlobalParams]
  PorousFlowDictator = dictator
[]

[Variables]
  [./pp]
    [./InitialCondition]
      type = RandomIC
      min = -1
      max = 1
    [../]
  [../]
[]

[Kernels]
  [./dot]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = pp
  [../]
  [./flux0]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    variable = pp
    gravity = '-1 0 0'
  [../]
[]

[Functions]
  [./ana_pp]
    type = ParsedFunction
    vars = 'g B p0 rho0'
    vals = '1 2 -1 1'
    value = '-B*log(exp(-p0/B)+g*rho0*x/B)' # expected pp at base
  [../]
[]

[BCs]
  [./z]
    type = PresetBC
    variable = pp
    boundary = right
    value = -1
  [../]
[]

[UserObjects]
  [./dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pp'
    number_fluid_phases = 1
    number_fluid_components = 1
  [../]
  [./pc_0]
    type = PorousFlowCapillaryPressureVG
    m = 0.5
    alpha = 1
  [../]
  [./pc_1]
    type = PorousFlowCapillaryPressureVG
    m = 0.6
    alpha = 2
  [../]
[]

[Modules]
  [./FluidProperties]
    [./simple_fluid]
      type = SimpleFluidProperties
      bulk_modulus = 2
      density0 = 1
      viscosity = 1
      thermal_expansion = 0
    [../]
  [../]
[]

[Materials]
  [./temperature]
    type = PorousFlowTemperature
  [../]
  [./ppss_0]
    type = PorousFlow1PhaseP
    block = 0
    porepressure = pp
    capillary_pressure = pc_0
  [../]
  [./ppss_1]
    type = PorousFlow1PhaseP
    block = 1
    porepressure = pp
    capillary_pressure = pc_1
  [../]
  [./massfrac]
    type = PorousFlowMassFraction
  [../]
  [./simple_fluid]
    type = PorousFlowSingleComponentFluid
    fp = simple_fluid
    phase = 0
  [../]
  [./porosity]
    type = PorousFlowPorosityConst
    porosity = 0.1
  [../]
  [./permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1 0 0  0 2 0  0 0 3'
  [../]
  [./relperm_0]
    type = PorousFlowRelativePermeabilityCorey
    block = 0
    n = 1
    phase = 0
  [../]
  [./relperm_1]
    type = PorousFlowRelativePermeabilityCorey
    block = 1
    n = 2
    phase = 0
  [../]
[]

[Postprocessors]
  [./pp_base]
    type = PointValue
    variable = pp
    point = '-1 0 0'
  [../]
  [./pp_analytical]
    type = FunctionValuePostprocessor
    function = ana_pp
    point = '-1 0 0'
  [../]
[]

[Preconditioning]
  active = andy
  [./andy]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = Newton
  dt = 1E6
  end_time = 1E6
[]

[Outputs]
  execute_on = 'timestep_end'
  file_base = grav01d
  csv = true
[]
