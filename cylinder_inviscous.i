[GlobalParams]
 	order = THIRD
 	family = MONOMIAL
  	
 	mach = 0.38
 	reynolds = 10.0
 	
 	variables = 'rho momentum_x momentum_y momentum_z rhoe'
[]

[Mesh]
  type = FileMesh
  file = grids/cylinder2.msh
  dim = 2
  
  block_id = 10
  block_name = 'fluid'
  
  boundary_id = '8 9'
  boundary_name = 'far_field wall'
	
	uniform_refine = 0 
  velocity = 0
[]

[Problem]
  fe_cache = true
  type = EulerProblem
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
		type = FarFieldPressure
		bc_type = far-field
		boundary = far_field 
		variable = rhoe
	[../]

	[./euler_wall]
		type = SlipWall
		boundary = wall 
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
  no_fe_reinit = true
  type = Transient
  solve_type = newton
  num_steps = 10
  l_tol = 1e-02
  #l_abs_step_tol = -1e-04
  l_max_its = 30
 	
  nl_max_its = 4
  nl_rel_tol = 1e-02

    petsc_options_iname = '-ksp_type  -pc_type -snes_lag_jacobian -snes_lag_preconditioner'
    petsc_options_value = 'gmres       bjacobi 10 10'
  [./TimeStepper]
    type = RatioTimeStepper
    dt = 0.001
    ratio = 2
    step = 2
    max_dt = 20	
  [../]
[]

[Postprocessors]
  [./runtime]
	  type = RunTime
	  time_type = alive
  [../]
  [./residuals]
    type = Residual
  [../]
[]

[Outputs]
	csv = true
	[./console]
		type = Console	
		perf_log = true
		linear_residuals = true
	  nonlinear_residuals =  true	
	[../]

	[./exodus]
		type = Exodus
		output_initial = true
	[../]
	
[]



