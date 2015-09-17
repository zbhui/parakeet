[Mesh]
  file = grids/cylinder2.msh
  
  block_id = 10
  block_name = 'fluid'
  
  boundary_id = '8 9'
  boundary_name = 'far_field wall'
  uniform_refine = 2
[]

[Problem]
  type = EulerProblem
 	mach = 3
 	reynolds = 40.0
  [./Variables]
    order = SECOND
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
	  [./euler_far_field]
		  type = FarFieldPressure
		  boundary = far_field 
	  [../]

	  [./euler_wall]
		  type = SlipWall
		  boundary = wall 
	  [../]
  [../]
[]

[ICs]
  type = CFDInitialCondition
[]


[Preconditioning]
  [./SMP]
    type = SMP
    full = true
    #off_diag_row ='rhoe '
    #off_diag_column = 'rhoe '
  [../]
[]

[Adaptivity]
  [./Indicators]
    [./error]
      type = VariableJumpIndicator
      variable = density
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

[Executioner]
  no_fe_reinit = true
  type = Transient
  solve_type = newton
  num_steps = 100
  l_tol = 1e-01
  #l_abs_step_tol = -1e-04
  l_max_its = 30
 	
  nl_max_its = 10
  nl_rel_tol = 1e-01

    petsc_options_iname = '-ksp_type  -pc_type -snes_lag_jacobian -snes_lag_preconditioner'
    petsc_options_value = 'gmres       bjacobi 1 1'
  [./TimeStepper]
    type = RatioTimeStepper
    dt = 0.010
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
		output_on = 'linear nonlinear'
	[../]

	[./exodus]
		type = Exodus
    output_on = 'initial timestep_end'
    interval = 1
	[../]
	
[]



