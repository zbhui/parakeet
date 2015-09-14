[Mesh]
  type = GeneratedMesh
  dim = 2
  
  nx = 20
  ny = 10  
  
  xmin = 0
  xmax = 4

  ymin = 0
  ymax = 2
[]

[Problem]
  type = CouetteFlowProblem
  jacobian_delay = 1
  [./Variables]
    order = FIRST
    family = MONOMIAL
    variables = 'rho momx momy momz rhoe'
  [../]

  [./BCs]
	  [./euler_far_field]
		type = CFDBC
		boundary = '0 1 2 3'
	  [../]
  [../]

  [./Kernels]
    type = CFDCellKernel
  [../]

  [./DGKernels]
    type = CFDFaceKernel
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

[Executioner]
  type = Transient
  solve_type = newton
  scheme = bdf2
  dt = 0.02
  num_steps = 100
  l_tol = 1e-04
  l_max_its = 100
 	
  nl_max_its = 100
  nl_rel_tol = 1e-04

  petsc_options_iname = '-ksp_type  -pc_type '
  petsc_options_value = 'gmres       bjacobi '

  [./TimeStepper]
    type = RatioTimeStepper
    dt = 0.010
    ratio = 2
    step = 2
    max_dt = 20000	
  [../]
[]

[Functions]
  [./exact_rho]
    type = IsoVortexExact
  [../]
[]

[Postprocessors]
  [./l2_err]
    type = ProblemElementalL2Error
  [../]
[]

[Outputs]
 	csv = true
	[./exodus]
		type = Exodus
		output_on = 'initial timestep_end'
	[../]
	
	[./console]
		type = Console	
		perf_log = true
		output_on = 'linear nonlinear'
	[../]
[]



