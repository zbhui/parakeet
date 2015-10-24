[Mesh]
  file = grids/N0012-coarse-quad.msh

  block_id = 0
  block_name = 'fluid'

  boundary_id = '1 4 2 3'
  boundary_name = 'far_top far_bottom wall_top wall_bottom'
[]

[Problem]
  type = EulerProblem
  mach = 0.85
  attack = 1.0
  reynolds = 40.0
  jacobian_delay = 1


  [./Kernels]
    type = CFDCellKernel
    shock_indicator = indicator
  [../]

  [./DGKernels]
    type = CFDFaceKernel
  [../]
[]

[Variables]
  order = FIRST
  family = MONOMIAL
  variables = 'density momx momy momz rhoe'
  type = CFDInitialCondition
[]



[BoundaryCondition]
  [./euler_far_field]
    type = FarFieldPressure
    boundary = '1 4'
  [../]

  [./euler_wall]
    type = SlipWall
    boundary = '2 3'
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Adaptivity]
  [./Indicators]
    [./indicator]
      type = VariableJumpIndicator
      variables = 'density momx momy momz rhoe'
      variable = density
      scale = 1
    [../]
  [../]
  [./Markers]
    [./marker]
      type = ErrorFractionMarker
      indicator = indicator
      coarsen = 0.7
      refine = 0.5
    [../]
  [../]
[]

[Executioner]
  no_fe_reinit = true
  type = Transient
  solve_type = NEWTON
  num_steps = 100
  l_tol = 1e-01
  l_max_its = 10

  nl_max_its = 5
  nl_rel_tol = 1E-01

  [./TimeStepper]
    type = RatioTimeStepper
    dt = 1E-01
    ratio = 2
    step = 2
    max_dt = 20000
  [../]
[]

[Postprocessors]
  [./residuals]
    type = Residual
  [../]
[]

[Outputs]
  csv = true
  [./console]
    type = Console
    perf_log = true
    execute_on = 'linear nonlinear'
  [../]

  [./exodus]
    type = Exodus
    execute_on = 'initial timestep_end'
    interval = 1
  [../]

[]
