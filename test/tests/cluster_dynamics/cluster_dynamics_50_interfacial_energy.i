[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1
[]

[Variables]
  [clusters]
    family = LAGRANGE
    order = FIRST
    components = 1000
  []
[]

[Functions]
  [index_fn]
    type = PiecewiseConstant
    x = '0 0.9'
    y = '0.0134 0'
    direction = left
  []
[]

[ICs]
  [u_ic]
    type = ArrayFunctionIndexIC
    variable = clusters
    function = index_fn
  []
[]

[NodalKernels]
  [clusters_dot]
    type = ArrayTimeDerivativeNodalKernel
    variable = clusters
  []
  [clusters_dyn]
    type = ClusterDynamicsNodalKernel
    variable = clusters
    generation = 0.0
    sink = 0.0

    rate_model = interfacial_energy
    diffusivity_model = arrhenius

    temperature = 563.15
    D0 = 6.3e-5
    Q_eV = 2.29
    radiation_enhanced_factor = 2.5e7
    sigma = 0.37
    atomic_volume = 1.1782924e-29
    Omega_kB_K = 6255.0
    DeltaS_kB = 0.866
  []
[]

[Postprocessors]
  [avg_cluster_radius]
    type = ClusterAverageRadius
    clusters = clusters
    r1 = 1.41e-10
    n_minimum = 11
  []
  [monomer_concentration]
    type = ClusterSizeConcentration
    clusters = clusters
    n_size = 1
  []
  [total_cluster_density]
    type = ClusterTotalDensity
    clusters = clusters
    n_minimum = 11
    atomic_volume = 1.1782924e-29
  []
[]

# [Executioner]
#   type = Transient
#   solve_type = PJFNK
#   petsc_options_iname='-pc_type'
#   petsc_options_value='lu'
#   # nl_rel_tol = 1e-8
#   # nl_abs_tol = 1e-9
#   # dtmin = 1e-6
#   # dtmax = 100
  
#   scheme = bdf2

#   l_tol = 1.0e-4
#   l_max_its = 30
#   nl_max_its = 20
#   nl_rel_tol = 1.0e-8
#   nl_abs_tol = 1.0e-16
#   end_time = 100000
#   automatic_scaling=true
#   compute_scaling_once = false


#   [TimeStepper]
#     type = IterationAdaptiveDT
#     dt = 0.000001
#     optimal_iterations = 6
#     iteration_window = 2
#     growth_factor = 1.1
#     cutback_factor = 0.5
#   []
# []

[Executioner]
  type = Transient
  solve_type = PJFNK
  #solve_type = NEWTON

  scheme = bdf2

  start_time = 0
  end_time = 100000

  l_tol = 1.0e-4
  l_max_its = 30
  nl_max_its = 20
  nl_rel_tol = 1.0e-6
  nl_abs_tol = 1.0e-12

  [./TimeStepper]
    type = IterationAdaptiveDT

    dt = 1e-4
    optimal_iterations = 12
    iteration_window = 2
  [../]
[]



# [Executioner]
#   type = Transient
#   solve_type = PJFNK
#   scheme = bdf2

#   l_tol = 1.0e-4
#   l_max_its = 30
#   nl_max_its = 20
#   nl_rel_tol = 1.0e-6
#   nl_abs_tol = 1.0e-9

#   dt = 1
#   end_time = 100000
# []


[Outputs]
  csv = true
[]
