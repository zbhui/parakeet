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
		variable = rho
	[../]		
[]

#面积分
[DGKernels]
	[./multi_dg_kernel]
		type = EulerFaceKernel
		variable = rho
	[../]
	[./multi_dg_kernel1]
		type = EulerFaceKernel
		variable = momentum_x
	[../]
	[./multi_dg_kernel2]
		type = EulerFaceKernel
		variable = momentum_y
	[../]
	[./multi_dg_kernel3]
		type = EulerFaceKernel
		variable = momentum_z
	[../]
	[./multi_dg_kernel4]
		type = EulerFaceKernel
		variable = rhoe
	[../]
[]

# 边界条件
[BCs]
	[./euler_far_field]
		type = IsoVortexBC
		boundary = '0 1 2 3'
		variable = rho
	[../]
	[./euler_far_field1]
		type = IsoVortexBC
		boundary = '0 1 2 3'
		variable = momentum_x
	[../]
	[./euler_far_field2]
		type = IsoVortexBC
		boundary = '0 1 2 3'
		variable = momentum_y
	[../]
	[./euler_far_field3]
		type = IsoVortexBC
		boundary = '0 1 2 3'
		variable = momentum_z
	[../]
	[./euler_far_field4]
		type = IsoVortexBC
		boundary = '0 1 2 3'
		variable = rhoe
	[../]	

[]

# 非线性系统求解
[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'bjacobi'
 	scheme = 'bdf2'
  dt = 0.001
  num_steps = 1
  
  # 线性迭代步的残差下降（相对）量级
 	l_tol = 1e-3
	# l_abs_step_tol = -1e-04
  # 最大线性迭代步	
 	l_max_its = 100
 	
 	# 最大非线性迭代步
 	nl_max_its = 10
 	# 非线性迭代的残值下降（相对）量级
 	nl_rel_tol = 1e-02
 	# 非线性迭代绝对残值
 	#nl_abs_tol = 1e-05
  	
	#abort_on_solve_fail = true	
  #end_time = 0.1
  
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
 	#use_displaced = true
	[./exodus]
		type = Exodus
		output_initial = true
		
		interval = 1 					#间隔
		oversample = true
		refinements = 0
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



