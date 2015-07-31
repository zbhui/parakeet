# 全局变量
[GlobalParams]
 	order = FIRST
 	family = MONOMIAL
  	
  mach = 0.1
  reynolds = 40.0
  	
	variable = rho
  variables = 'rho momentum_x momentum_y momentum_z rhoe'
	
[]

# 网格
[Mesh]
  type = FileMesh
  file = grids/cylinder_tri.msh
  dim = 2
  
  block_id = 10
  block_name = 'fluid'
  
  boundary_id = '8 9'
  boundary_name = 'far_field wall'
#  second_order = 0
  uniform_refine = 0
[]

# 变量
[Variables]
	active = 'rho momentum_x momentum_y momentum_z rhoe'

	[./rho]
		[./InitialCondition] #初始条件
			type = CFDPassFlowIC
		[../]
 	[../]
  
 	[./momentum_x]
		[./InitialCondition] #初始条件
			type = CFDPassFlowIC
		[../]
 	[../]
  	
 	[./momentum_y]
		[./InitialCondition] #初始条件
			type = CFDPassFlowIC
		[../]
 	[../]
  	
 	[./momentum_z]
		[./InitialCondition] #初始条件
			type = CFDPassFlowIC
		[../]
 	[../]
  	
  [./rhoe]
		[./InitialCondition] #初始条件
			type = CFDPassFlowIC
		[../]
  [../]	
		
[]

  
# 体积分
[Kernels]
	[./mass_time]
		type =TimeDerivative
		variable = rho
	[../]	
	
	[./x-momentumum_time]
		type = TimeDerivative
		variable = momentum_x
	[../]
		
	[./y-momentumum_time]
		type = TimeDerivative
		variable = momentum_y
	[../]
	
	[./z-momentumum_time]
		type = TimeDerivative
		variable = momentum_z
	[../]	
		
	[./total-energy_time]
		type = TimeDerivative
		variable = rhoe
	[../]
	
	[./multi_kernel]
		type = NavierStokesCellKernel
	[../]		
[]

#面积分
[DGKernels]
	[./multi_dg_kernel]
		type = NavierStokesFaceKernel
	[../]
[]

# 边界条件
[BCs]
	[./euler_far_field]
		type = NavierStokesBC
		bc_type = far-field
		boundary = far_field 
	[../]
		

	[./euler_wall]
		type = NavierStokesBC
		bc_type = wall
		boundary = wall 
	[../]

[]

# 非线性系统求解
[Executioner]
  type = Transient
  solve_type = PJFNK
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'jacobi'
	 	#scheme = 'rk-2'
  	dt = 0.1
  	num_steps = 1000
  
    # 线性迭代步的残差下降（相对）量级
 	l_tol = 1e-01
 # l_abs_step_tol = -1e-04
   # 最大线性迭代步	
 	l_max_its = 10
 	
 	# 最大非线性迭代步
 	nl_max_its = 20
 	# 非线性迭代的残值下降（相对）量级
  	nl_rel_tol = 1e-01
  	# 非线性迭代绝对残值
  	#nl_abs_tol = 1e-03
  	
	 #abort_on_solve_fail = true	
  #end_time = 0.1
  
  	[./Adaptivity]
  	
 	[../]
[]

[Functions]
  [./exact_rho]
    type = IsoVortexExact
  [../]
[]

[Postprocessors]
  #[./h]
   # type = AverageElementSize
    #variable = rho
 # [../]

 # [./dofs]
  #  type = NumDOFs
 # [../]

  [./l2_err]
    type = ElementL2Error
    variable = rho
    function = exact_rho
  [../]
  
  #[./nodes]
  #  type = NumNodes
  #[../]

  #[./elements]
   # type = NumElems
  #[../]

 # [./residuals]
 #   type = Residual
 # [../]
  
 #[./integral_left]
 #  type = ElementIntegralVariablePostprocessor
 #  variable = rho
 #[../]  
[]

# 输出和后处理
[Outputs]
	[./console]
		type = Console	
		perf_log = true
		linear_residuals = true
	  nonlinear_residuals =  true	
	[../]

	[./exodus]
		type = Exodus
		output_initial = true
		
		interval = 1  					#间隔
		oversample = true
		refinements = 0
	[../]
	
[]



