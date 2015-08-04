[GlobalParams]
 	order = SECOND
 	family = MONOMIAL
  	
 	mach = 0.1
 	reynolds = 10.0
  	
 	variables = 'rho momentum_x momentum_y momentum_z rhoe'
[]

[Mesh]
  type = GeneratedMesh
  dim = 2
  
  nx = 20
  ny = 20  
  
  xmin = -10
  xmax = 0

  ymin = -10
  ymax = 0
  
  uniform_refine = 0
  second_order = false
  
  block_id = '0'
  block_name = 'fluid'

[]

[Problem]
  type = IsoVortexProblem
[]

[Variables]
	[./rho]
		[./InitialCondition] 
			type = CFDInitialCondition
      component = 0
		[../]
  [../]

 	[./momentum_x]
		[./InitialCondition] 
			type = CFDInitialCondition
      component = 1
		[../]
  [../]
  
 	[./momentum_y]
		[./InitialCondition] 
			type = CFDInitialCondition
      component = 2
		[../]
  [../]
  	
  [./momentum_z]
		[./InitialCondition] 
			type = CFDInitialCondition
      component = 3
		[../]
  [../]
  	
  [./rhoe]
		[./InitialCondition] 
			type = CFDInitialCondition
      component = 4
		[../]
  [../]	
		
[]

[Kernels]
	[./mass_time]
		type =TimeDerivative
		variable = rho
	[../]	
	
	[./x-momentum_time]
		type = TimeDerivative
		variable = momentum_x
	[../]
		
	[./y-momentum_time]
		type = TimeDerivative
		variable = momentum_y
	[../]
	
	[./z-momentum_time]
		type = TimeDerivative
		variable = momentum_z
	[../]	
		
	[./total-energy_time]
		type = TimeDerivative
		variable = rhoe
	[../]
	
	[./multi_kernel]
		type = CFDCellKernel
		variable = rhoe
	[../]		
[]

[DGKernels]
	[./multi_dg_kernel]
		type = CFDFaceKernel
		variable = rhoe
	[../]
[]

[BCs]
	[./euler_far_field]
		type = CFDBC
		boundary = '0 1 2 3'
		variable = rhoe
	[../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
    #off_diag_row ='rhoe '
    #off_diag_column = 'rhoe '
  [../]
[]

[Executioner]
  type = Transient
  solve_type = newton
  scheme = bdf2
  dt = 0.02
  num_steps = 1
  l_tol = 1e-04
  l_max_its = 100
 	
  nl_max_its = 100
  nl_rel_tol = 1e-04

  petsc_options_iname = '-ksp_type  -pc_type -snes_lag_jacobian -snes_lag_preconditioner'
  petsc_options_value = 'gmres       bjacobi 20 20'
[]

[Functions]
  [./exact_rho]
    type = IsoVortexExact
  [../]
[]

[Postprocessors]
  [./l2_err]
    type = ElementL2Error
    variable = rho
    function = exact_rho
  [../]
[]

[Outputs]
 	csv = true
	gnuplot = true	
	[./exodus]
		type = Exodus
		output_initial = true
	[../]
	
	[./console]
		type = Console	
		perf_log = true
		linear_residuals = true
	  	nonlinear_residuals =  true	
		#verbose = true
    	#setup_log_early = true
    	#time_precision = 6
    	#fit_mode = 100
	[../]
[]



