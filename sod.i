[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 100
[]

[Problem]
  type = Riemann1DProblem
  sub_type = sod
  mach = 0.85
  reynolds = 40.0
  [./Variables]
    order = FOURTH
    family = MONOMIAL
    variables = 'density momx momy momz rhoe'
  [../]

  [./Kernels]
    type = CFDCellKernel
  [../]

  [./DGKernels]
    type = CFDFaceKernel
  [../]


  [./BCs]
    [./exact_bc]
      type = CFDBC
      boundary = '0 1'
    [../]
  [../]
[]

[ICs]
  type = CFDInitialCondition 
[]

[Adaptivity]
  [./Indicators]
    [./error]
      type = VariableJumpIndicator
      variable = density
      execute_on = initial
    [../]
  [../]
  [./Markers]
    [./marker]
      type = ErrorFractionMarker
      indicator = error
      coarsen = 0.7
      refine = 0.9
    [../]
  [../]
[]


[Postprocessors]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_type  -pc_type'
    petsc_options_value = 'gmres       bjacobi'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = newton
  dt = 0.001
  num_steps = 2000
  scheme = bdf2
  l_tol = 1e-01
  l_max_its = 10
 	
  nl_max_its = 10
  nl_rel_tol = 1e-03
  end_time = 0.2
[]

[Outputs]
  [./exodus]
    type = Exodus
    interval = 1 					
  [../]
	
  [./console]
    type = Console	
    perf_log = true
    execute_on = 'linear nonlinear'
  [../]
[]

