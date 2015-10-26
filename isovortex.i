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
  flux_type = HLLC-PV
  [./Kernels]
    type = CFDCellKernel
    #indicator = error
  [../]

  [./DGKernels]
    type = CFDFaceKernel
  [../]
[]

[Variables]
    type =  CFDInitialCondition
    order = FIRST
    family = MONOMIAL
    variables = 'rho momx momy momz rhoe'
[]

[BoundaryCondition]
    [./euler_far_field]
        type =  CFDBC
        boundary = '0 1 2 3'
    [../]
[]


[ArtificialViscosity]
  [./Indicators]
     [./error]
       type = GradientJumpIndicator
       variable = rho

     [../]
  [../]
  [./Markers]
     [./marker]
       type = ErrorMaxFractionMarker
       indicator = error
       coarsen = 0.7
       refine = 0.9
     [../]
   [../]
[]

 [Preconditioning]
   [./SMP]
     type = SMP
     full = true
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
    type = ProblemElementalL2Error
  [../]
[]

[Outputs]
 	csv = true
	gnuplot = true
	[./exodus]
		type = Exodus
		execute_on = 'initial timestep_end'
	[../]

	[./console]
		type = Console
		perf_log = true
		execute_on = 'linear nonlinear'
	[../]
[]
