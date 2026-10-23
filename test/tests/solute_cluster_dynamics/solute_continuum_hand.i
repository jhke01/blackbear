# The problem of solute_continuum_action.i built from ContinuumClusterDynamicsNodalKernel, giving the reference solution that the
# [SoluteClusterDynamics] action reproduces with tail = fokker_planck.

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1
[]
[Problem]
  use_hash_table_matrix_assembly = true
  restore_original_nonzero_pattern = false
[]
[Variables]
  [SoluteClusterDynamics_clusters]
    family = LAGRANGE
    order = FIRST
    components = 132
    array = true
  []
[]
[Functions]
  [SoluteClusterDynamics_index_fn]
    type = PiecewiseConstant
    x = '0 0.9'
    y = '0.0134 0'
    direction = left
  []
[]
[ICs]
  [SoluteClusterDynamics_ic]
    type = ArrayFunctionIndexIC
    variable = SoluteClusterDynamics_clusters
    function = SoluteClusterDynamics_index_fn
  []
[]
[NodalKernels]
  [SoluteClusterDynamics_dot]
    type = ArrayTimeDerivativeNodalKernel
    variable = SoluteClusterDynamics_clusters
  []
  [SoluteClusterDynamics_dyn]
    type = ContinuumClusterDynamicsNodalKernel
    variable = SoluteClusterDynamics_clusters
    generation = 0
    sink = 0
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
    num_cluster_sizes = 2000
    switch_size = 20
  []
[]
[Postprocessors]
  [SoluteClusterDynamics_monomer_concentration]
    type = ClusterSizeConcentration
    clusters = SoluteClusterDynamics_clusters
    n_size = 1
  []
  [SoluteClusterDynamics_total_density]
    type = ClusterContinuumMoment
    clusters = SoluteClusterDynamics_clusters
    quantity = density
    num_cluster_sizes = 2000
    switch_size = 20
    n_minimum = 1
    r1 = 1.4116319815634138e-10
    radius_exponent = 0.3333333333333333
    atomic_volume = 1.1782924e-29
  []
  [SoluteClusterDynamics_resolvable_density]
    type = ClusterContinuumMoment
    clusters = SoluteClusterDynamics_clusters
    quantity = density
    num_cluster_sizes = 2000
    switch_size = 20
    n_minimum = 45
    r1 = 1.4116319815634138e-10
    radius_exponent = 0.3333333333333333
    atomic_volume = 1.1782924e-29
  []
  [SoluteClusterDynamics_resolvable_avg_radius]
    type = ClusterContinuumMoment
    clusters = SoluteClusterDynamics_clusters
    quantity = average_radius
    num_cluster_sizes = 2000
    switch_size = 20
    n_minimum = 45
    r1 = 1.4116319815634138e-10
    radius_exponent = 0.3333333333333333
    atomic_volume = 1.1782924e-29
  []
[]
[Preconditioning]
  [smp]
    type = SMP
    full = false
  []
[]
[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_ordering_type'
  petsc_options_value = 'lu       nd'
  line_search = none
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-14
  dt = 0.05
  num_steps = 20
[]
[Outputs]
  csv = true
[]
