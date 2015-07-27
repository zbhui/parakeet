# 全局变量
[GlobalParams]
 	order = THIRD
 	family = MONOMIAL
	  	
  mach = 0.2
  reynolds = 100.0
  	
	variable = rho
 	variables = 'rho momentum_x momentum_y momentum_z rhoe'
[]

# 网格
[Mesh]
  type = GeneratedMesh
  dim = 2
  
  nx = 20
  ny = 10  
  
  xmin = 0
  xmax = 4

  ymin = 0
  ymax = 2

	#uniform_refine = 0
	#second_order = false
	#elem_type = TRI6
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
	[./multi_bc]
		type = CouetteFlowBC
		boundary = 'left right bottom top'
	[../]			
[]

# 非线性系统求解
[Executioner]
  type = Transient
  solve_type = PJFNK
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'bjacobi'
  dt = 1
  num_steps = 10000
  
 	l_tol = 1e-01
 	l_max_its = 10
 	
 	nl_max_its = 10
  nl_rel_tol = 1e-01
  #nl_abs_tol = 1e-05
  

[]

[Functions]
  [./exact_rho]
    type = CouetteFlowExact
  [../]
[]


# 输出和后处理
[Outputs]
  #csv = true
 	
	[./exodus]
		type = Exodus
		output_initial = true	
		interval = 1 					#间隔
	[../]
	
	[./console]
		type = Console	
		perf_log = true
		linear_residuals = true
	  	nonlinear_residuals =  true	
		interval = 1 
	[../]
[]



