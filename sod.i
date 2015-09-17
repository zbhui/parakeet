[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 100
[]

[Problem]
  type = SodProblem
   mach = 0.85
   reynolds = 40.0
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
  dt = 0.002
  num_steps = 1000
  scheme = bdf2
  l_tol = 1e-01
  l_max_its = 100
 	
  nl_max_its = 100
  nl_rel_tol = 1e-03
[]

[Outputs]
  [./exodus]
    type = Exodus
    interval = 1 					
  [../]
	
  [./console]
    type = Console	
    perf_log = true
    output_on = 'linear nonlinear'
  [../]
[]

