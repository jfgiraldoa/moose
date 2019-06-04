# Test for small strain euler beam vibration in y direction

# An impulse load is applied at the end of a cantilever beam of length 4m.
# The beam is massless with a lumped mass at the end of the beam. The lumped
# mass also has a moment of inertia associated with it.
# The properties of the cantilever beam are as follows:
# Young's modulus (E) = 1e4
# Shear modulus (G) = 4e7
# Shear coefficient (k) = 1.0
# Cross-section area (A) = 0.01
# Iy = 1e-4 = Iz
# Length (L)= 4 m
# mass (m) = 0.01899772
# Moment of inertia of lumped mass:
# Ixx = 0.2
# Iyy = 0.1
# Izz = 0.1
# mass proportional damping coefficient (eta) = 0.1

# For this beam, the dimensionless parameter alpha = kAGL^2/EI = 6.4e6
# Therefore, the beam behaves like a Euler-Bernoulli beam.

# The displacement time history from this analysis matches with that obtained from Abaqus.

# Values from the first few time steps are as follows:
# time   disp_y              vel_y               accel_y
# 0.0    0.0                 0.0                 0.0
# 0.1    0.001278249649738   0.025564992994761   0.51129985989521
# 0.2    0.0049813090917644  0.048496195845768  -0.052675802875074
# 0.3    0.0094704658873002  0.041286940064947  -0.091509312741339
# 0.4    0.013082280729802   0.03094935678508   -0.115242352856
# 0.5    0.015588313103503   0.019171290688959  -0.12031896906642

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 10
  xmin = 0.0
  xmax = 4.0
  displacements = 'disp_x disp_y disp_z'
[]

[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
  [../]
  [./rot_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./rot_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./rot_z]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[AuxVariables]
  [./vel_x]
  order = FIRST
  family = LAGRANGE
  [../]
  [./vel_y]
  order = FIRST
  family = LAGRANGE
  [../]
  [./vel_z]
  order = FIRST
  family = LAGRANGE
  [../]
  [./accel_x]
  order = FIRST
  family = LAGRANGE
  [../]
  [./accel_y]
  order = FIRST
  family = LAGRANGE
  [../]
  [./accel_z]
  order = FIRST
  family = LAGRANGE
  [../]
  [./rot_vel_x]
  order = FIRST
  family = LAGRANGE
  [../]
  [./rot_vel_y]
  order = FIRST
  family = LAGRANGE
  [../]
  [./rot_vel_z]
  order = FIRST
  family = LAGRANGE
  [../]
  [./rot_accel_x]
  order = FIRST
  family = LAGRANGE
  [../]
  [./rot_accel_y]
  order = FIRST
  family = LAGRANGE
  [../]
  [./rot_accel_z]
  order = FIRST
  family = LAGRANGE
  [../]
[]

[AuxKernels]
  [./accel_x] # These auxkernels are only to check output
    type = TestNewmarkTI
    displacement = disp_x
    variable = accel_x
    first = false
  [../]
  [./accel_y]
    type = TestNewmarkTI
    displacement = disp_y
    variable = accel_y
    first = false
  [../]
  [./accel_z]
    type = TestNewmarkTI
    displacement = disp_z
    variable = accel_z
    first = false
  [../]
  [./vel_x]
    type = TestNewmarkTI
    displacement = disp_x
    variable = vel_x
  [../]
  [./vel_y]
    type = TestNewmarkTI
    displacement = disp_y
    variable = vel_y
  [../]
  [./vel_z]
    type = TestNewmarkTI
    displacement = disp_z
    variable = vel_z
  [../]
  [./rot_accel_x]
    type = TestNewmarkTI
    displacement = rot_x
    variable = rot_accel_x
    first = false
  [../]
  [./rot_accel_y]
    type = TestNewmarkTI
    displacement = rot_y
    variable = rot_accel_y
    first = false
  [../]
  [./rot_accel_z]
    type = TestNewmarkTI
    displacement = rot_z
    variable = rot_accel_z
    first = false
  [../]
  [./rot_vel_x]
    type = TestNewmarkTI
    displacement = rot_x
    variable = rot_vel_x
  [../]
  [./rot_vel_y]
    type = TestNewmarkTI
    displacement = rot_y
    variable = rot_vel_y
  [../]
  [./rot_vel_z]
    type = TestNewmarkTI
    displacement = rot_z
    variable = rot_vel_z
  [../]
[]

[BCs]
  [./fixx1]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]
  [./fixy1]
    type = DirichletBC
    variable = disp_y
    boundary = left
    value = 0.0
  [../]
  [./fixz1]
    type = DirichletBC
    variable = disp_z
    boundary = left
    value = 0.0
  [../]
  [./fixr1]
    type = DirichletBC
    variable = rot_x
    boundary = left
    value = 0.0
  [../]
  [./fixr2]
    type = DirichletBC
    variable = rot_y
    boundary = left
    value = 0.0
  [../]
  [./fixr3]
    type = DirichletBC
    variable = rot_z
    boundary = left
    value = 0.0
  [../]
[]

[NodalKernels]
  [./force_y2]
    type = UserForcingFunctionNodalKernel
    variable = disp_y
    boundary = right
    function = force
  [../]
  [./x_inertial]
    type = NodalTranslationalInertia
    variable = disp_x
    boundary = right
    mass = 0.01899772
    eta = 0.1
  [../]
  [./y_inertial]
    type = NodalTranslationalInertia
    variable = disp_y
    boundary = right
    mass = 0.01899772
    eta = 0.1
  [../]
  [./z_inertial]
    type = NodalTranslationalInertia
    variable = disp_z
    boundary = right
    mass = 0.01899772
    eta = 0.1
  [../]
  [./rot_x_inertial]
    type = NodalRotationalInertia
    variable = rot_x
    rotations = 'rot_x rot_y rot_z'
    boundary = right
    Ixx = 2e-1
    Iyy = 1e-1
    Izz = 1e-1
    eta = 0.1
    component = 0
  [../]
  [./rot_y_inertial]
    type = NodalRotationalInertia
    variable = rot_y
    rotations = 'rot_x rot_y rot_z'
    boundary = right
    Ixx = 2e-1
    Iyy = 1e-1
    Izz = 1e-1
    eta = 0.1
    component = 1
  [../]
  [./rot_z_inertial]
    type = NodalRotationalInertia
    variable = rot_z
    rotations = 'rot_x rot_y rot_z'
    boundary = right
    Ixx = 2e-1
    Iyy = 1e-1
    Izz = 1e-1
    eta = 0.1
    component = 2
  [../]
[]

[Functions]
  [./force]
    type = PiecewiseLinear
    x = '0.0 0.1 0.2 10.0'
    y = '0.0 1e-2  0.0  0.0'
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
  solve_type = NEWTON

  petsc_options_iname = '-ksp_type -pc_type'
  petsc_options_value = 'preonly   lu'

  start_time = 0.0
  dt = 0.1
  end_time = 5.0
  timestep_tolerance = 1e-6

  # Time integrator scheme
  scheme = "newmark-beta"
[]

[Kernels]
  [./solid_disp_x]
    type = StressDivergenceBeam
    block = '0'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    component = 0
    variable = disp_x
  [../]
  [./solid_disp_y]
    type = StressDivergenceBeam
    block = '0'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    component = 1
    variable = disp_y
  [../]
  [./solid_disp_z]
    type = StressDivergenceBeam
    block = '0'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    component = 2
    variable = disp_z
  [../]
  [./solid_rot_x]
    type = StressDivergenceBeam
    block = '0'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    component = 3
    variable = rot_x
  [../]
  [./solid_rot_y]
    type = StressDivergenceBeam
    block = '0'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    component = 4
    variable = rot_y
  [../]
  [./solid_rot_z]
    type = StressDivergenceBeam
    block = '0'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    component = 5
    variable = rot_z
  [../]
[]

[Materials]
  [./elasticity]
    type = ComputeElasticityBeam
    youngs_modulus = 1.0e4
    poissons_ratio = -0.999875
    shear_coefficient = 1.0
    block = 0
  [../]
  [./strain]
    type = ComputeIncrementalBeamStrain
    block = '0'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    area = 0.01
    Ay = 0.0
    Az = 0.0
    Iy = 1.0e-4
    Iz = 1.0e-4
    y_orientation = '0.0 1.0 0.0'
  [../]
  [./stress]
    type = ComputeBeamResultants
    block = 0
  [../]
[]

[Postprocessors]
  [./disp_x]
    type = PointValue
    point = '4.0 0.0 0.0'
    variable = disp_x
  [../]
  [./disp_y]
    type = PointValue
    point = '4.0 0.0 0.0'
    variable = disp_y
  [../]
  [./vel_y]
    type = PointValue
    point = '4.0 0.0 0.0'
    variable = vel_y
  [../]
  [./accel_y]
    type = PointValue
    point = '4.0 0.0 0.0'
    variable = accel_y
  [../]
[]

[Outputs]
  file_base = "dyn_euler_small_added_mass_inertia_damping_out"
  exodus = true
  csv = true
  perf_graph = true
[]
