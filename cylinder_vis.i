[Mesh]
  file = grids/cylinder2.msh

  block_id = 10
  block_name = 'fluid'

  boundary_id = '8 9'
  boundary_name = 'far_field wall'
  uniform_refine = 1
[]

[Problem]
  type = NavierStokesProblem
  mach = 0.2
  reynolds = 40.0
  [./Variables]
    order = SECOND
    family = MONOMIAL
    variables = 'density momx momy momz rhoe'
  [../]
  [./AuxVariables]
     order = FIRST
     family = MONOMIAL
    type = NSAuxVariable
    aux_variables = 'p m'
  [../]

  [./Kernels]
    type = CFDCellKernel
  [../]

  [./DGKernels]
    type = CFDFaceKernel
  [../]

  [./BCs]
    [./far_field]
    type = FarFieldPressure
    boundary = far_field
  [../]

  [./wall]
    type = IsoThermalWall
    boundary = wall
  [../]

[]

[ICs]
  type = CFDInitialCondition
[]


[Preconditioning]
  [./SMP]
    type = SMP
    full = true
  [../]
[]

[Adaptivity]
  [./Indicators]
    [./error]
      type = VariableJumpIndicator
      variables = 'density momx momy momz rhoe'
      variable = density
      scale = 2
    [../]
  [../]
  [./Markers]
    [./marker]
      type = ErrorFractionMarker
      indicator = error
      coarsen = 0.7
      refine = 0.5
    [../]
  [../]
[]

[Executioner]
  type = Transient
  solve_type = newton
  num_steps = 100
  l_tol = 1e-01
  l_max_its = 30

  nl_max_its = 10
  nl_rel_tol = 1e-01

  [./TimeStepper]
    type = RatioTimeStepper
    dt = 1.0
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
