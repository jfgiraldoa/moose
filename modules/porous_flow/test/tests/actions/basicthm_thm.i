# PorousFlowBasicTHM action with coupling_type = ThermoHydroMechanical

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 3
  xmax = 10
  ymax = 3
[]

[MeshModifiers]
  [./aquifer]
    type = SubdomainBoundingBox
    block_id = 1
    bottom_left = '0 1 0'
    top_right = '10 2 0'
  [../]
  [./injection_area]
    type = SideSetsAroundSubdomain
    block = 1
    new_boundary = 'injection_area'
    normal = '-1 0 0'
    depends_on = 'aquifer'
  [../]
  [./outflow_area]
    type = SideSetsAroundSubdomain
    block = 1
    new_boundary = 'outflow_area'
    normal = '1 0 0'
    depends_on = 'aquifer'
  [../]
  [./rename]
    type = RenameBlock
    old_block_id = '0 1'
    new_block_name = 'caprock aquifer'
    depends_on = 'injection_area'
  [../]
[]

[GlobalParams]
  PorousFlowDictator = dictator
  displacements = 'disp_x disp_y'
  biot_coefficient = 1.0
[]

[Variables]
  [./porepressure]
    initial_condition = 1e6
  [../]
  [./temperature]
    initial_condition = 293
    scaling = 1e-6
  [../]
  [./disp_x]
    scaling = 1e-6
  [../]
  [./disp_y]
    scaling = 1e-6
  [../]
[]

[PorousFlowBasicTHM]
  porepressure = porepressure
  temperature = temperature
  coupling_type = ThermoHydroMechanical
  gravity = '0 0 0'
  fp = simple_fluid
  thermal_eigenstrain_name = thermal_contribution
  use_displaced_mesh = false
  add_stress_aux = false
[]

[BCs]
  [./constant_injection_porepressure]
    type = PresetBC
    variable = porepressure
    value = 1.5e6
    boundary = injection_area
  [../]
  [./constant_injection_temperature]
    type = PresetBC
    variable = temperature
    value = 313
    boundary = injection_area
  [../]
  [./constant_outflow_porepressure]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure
    boundary = outflow_area
    pt_vals = '0 1e9'
    multipliers = '0 1e9'
    flux_function = 1e-6
    PT_shift = 1e6
  [../]
  [./constant_outflow_temperature]
    type = PresetBC
    variable = temperature
    value = 293
    boundary = outflow_area
  [../]
  [./top_bottom]
    type = PresetBC
    variable = disp_y
    value = 0
    boundary = 'top bottom'
  [../]
  [./right]
    type = PresetBC
    variable = disp_x
    value = 0
    boundary = right
  [../]
[]

[Modules]
  [./FluidProperties]
    [./simple_fluid]
      type = SimpleFluidProperties
    [../]
  [../]
[]

[Materials]
  [./porosity]
    type = PorousFlowPorosity
    porosity_zero = 0.1
  [../]
  [./biot_modulus]
    type = PorousFlowConstantBiotModulus
    biot_coefficient = 0.8
    solid_bulk_compliance = 2e-7
    fluid_bulk_modulus = 1e7
  [../]
  [./permeability_aquifer]
    type = PorousFlowPermeabilityConst
    block = aquifer
    permeability = '1e-13 0 0   0 1e-13 0   0 0 1e-13'
  [../]
  [./permeability_caprock]
    type = PorousFlowPermeabilityConst
    block = caprock
    permeability = '1e-15 0 0   0 1e-15 0   0 0 1e-15'
  [../]
  [./thermal_expansion]
    type = PorousFlowConstantThermalExpansionCoefficient
    drained_coefficient = 0.003
    fluid_coefficient = 0.0002
  [../]
  [./rock_internal_energy]
    type = PorousFlowMatrixInternalEnergy
    density = 2500.0
    specific_heat_capacity = 1200.0
  [../]
  [./thermal_conductivity]
    type = PorousFlowThermalConductivityIdeal
    dry_thermal_conductivity = '10 0 0  0 10 0  0 0 10'
    block = 'caprock aquifer'
  [../]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 5e9
    poissons_ratio = 0.0
  [../]
  [./strain]
    type = ComputeSmallStrain
    eigenstrain_names = thermal_contribution
  [../]
  [./thermal_contribution]
    type = ComputeThermalExpansionEigenstrain
    temperature = temperature
    thermal_expansion_coeff = 0.001
    eigenstrain_name = thermal_contribution
    stress_free_temperature = 293
  [../]
  [./stress]
    type = ComputeLinearElasticStress
  [../]
[]

[Preconditioning]
  [./basic]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = Newton
  end_time = 1e4
  dt = 1e3
  nl_abs_tol = 1e-12
  nl_rel_tol = 1E-10
[]

[Outputs]
  exodus = true
[]
