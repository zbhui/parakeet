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

[Variables]
	[./rho]
		[./InitialCondition] 
			type = IsoVortexIC
		[../]
  [../]

 	[./momentum_x]
		[./InitialCondition] 
			type = IsoVortexIC
		[../]
  [../]
  
 	[./momentum_y]
		[./InitialCondition] 
			type = IsoVortexIC
		[../]
  [../]
  	
  [./momentum_z]
		[./InitialCondition] 
			type = IsoVortexIC
		[../]
  [../]
  	
  [./rhoe]
		[./InitialCondition] 
			type = IsoVortexIC
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
		type = EulerCellKernel
		variable = rhoe
	[../]		
[]

[DGKernels]
	[./multi_dg_kernel]
		type = EulerFaceKernel
		variable = rhoe
	[../]
[]

[BCs]
	[./euler_far_field]
		type = IsoVortexBC
		boundary = '0 1 2 3'
		variable = rhoe
	[../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    #full = true
    off_diag_row ='rhoe '
    off_diag_column = 'rhoe '
  [../]
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  scheme = bdf2
  dt = 0.02
  num_steps = 10
  l_tol = 1e-04
  l_max_its = 100
 	
  nl_max_its = 100
  nl_rel_tol = 1e-02
  #nl_abs_tol = 1e-05
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



