[Mesh]
  file = ./grids/forward_step.e
[]

[Problem]
  type = ForwardStepProblem
  jacobian_delay = 1
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
	[./inner]
	  type = CFDBC
	  boundary = 'inner '
	[../]
      [./wall]
         type = Symmetric
         boundary = 'wall'
      [../]
      [./out]
          type = FarFieldRiemann
          boundary = 'out'
      [../]
  [../]
[]

[ICs]
  type = CFDInitialCondition 
  # constant_ic = true
[]

[Adaptivity]
  [./Indicators]
	[./error]
	  type = FluxJumpIndicator
	  variables = 'density momx momy momz rhoe'
	  variable = density
	  scale = 0.5
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
    scheme = bdf2
    num_steps = 2000
    l_tol = 1e-01
    l_max_its = 10
  	
    nl_max_its = 5
    nl_rel_tol = 1e-02
    end_time = 4

    [./TimeStepper]
      type = RatioTimeStepper
      dt = 0.0010
      ratio = 2
      step = 2
      max_dt = 0.01  
  [../]
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

