[Mesh]
  type = GeneratedMesh
  dim = 2
  xmax  = 2
  nx = 80
  ny = 40
[]

[Problem]
  type = ShockVortexProblem
  jacobian_delay = 3
  [./Variables]
	order = FIRST
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
	  boundary = '0 1 2 3'
	[../]
  [../]
[]

[ICs]
  type = CFDInitialCondition 
[]

[Adaptivity]
  [./Indicators]
	[./error]
	  type = FluxJumpIndicator
	  variables = 'density momx momy momz rhoe'
	  variable = density
	  scale = 0.01
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
    dt = 0.01
    scheme = bdf2
    num_steps = 2000
    l_tol = 1e-01
    l_max_its = 10
  	
    nl_max_its = 10
    nl_rel_tol = 1e-01
    end_time = 1
  
[]

[Outputs]
  [./exodus]
	type = Exodus
	interval = 1 	
	refinements = 0				
  [../]
	
  [./console]
	type = Console	
	perf_log = true
	execute_on = 'linear nonlinear'
  [../]
[]

