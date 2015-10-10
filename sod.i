[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 100
[]

[Problem]
  type = Riemann1DProblem
  sub_type = sod
  [./Variables]
    order = SECOND
    family = MONOMIAL
    variables = 'density momx momy momz rhoe'
  [../]
  [./AuxVariables]
    [./axu]
      order = THIRD
      family = MONOMIAL
      type = NSAuxVariable
      variables = 'density momx momy momz rhoe'
      aux_variables = 'p mach '
      execute_on = 'initial timestep_end'
    [../]
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
      boundary = '0 1'
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
      scale = 0.2
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
  dt = 0.001
  scheme = bdf2
  num_steps = 2000
  l_tol = 1e-01
  l_max_its = 10
 	
  nl_max_its = 10
  nl_rel_tol = 1e-03
  end_time = 0.2
  
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

