# The problem of defect_dynamics_fokker_planck_action.i built from ContinuumLoopNodalKernel, CascadeSpectrumNodalKernel, DislocationRecoveryNodalKernel and PointDefectRecombinationNodalKernel, giving the reference solution that the [DefectDynamics] action reproduces with tail = fokker_planck.
[Mesh]
  type = SingleNodeMesh
  position = '0 0 0'
[]
[Problem]
  use_hash_table_matrix_assembly = true
  restore_original_nonzero_pattern = false
[]
[Variables]
  [vacancy_clusters]
    family = LAGRANGE
    order = FIRST
    components = 132
    array = true
  []
  [interstitial_clusters]
    family = LAGRANGE
    order = FIRST
    components = 132
    array = true
  []
[]
[Functions]
  [vacancy_index_fn]
    type = PiecewiseConstant
    x = '0 0.9'
    y = '2.4961095795652884e-14 0'
    direction = left
  []
  [interstitial_index_fn]
    type = PiecewiseConstant
    x = '0 0.9'
    y = '2.205126527820641e-34 0'
    direction = left
  []
[]
[ICs]
  [vacancy_ic]
    type = ArrayFunctionIndexIC
    variable = vacancy_clusters
    function = vacancy_index_fn
  []
  [interstitial_ic]
    type = ArrayFunctionIndexIC
    variable = interstitial_clusters
    function = interstitial_index_fn
  []
[]
[NodalKernels]
  [vacancy_dot]
    type = ArrayTimeDerivativeNodalKernel
    variable = vacancy_clusters
  []
  [vacancy_dyn]
    type = ContinuumLoopNodalKernel
    variable = vacancy_clusters
    generation = 0.267
    temperature = 603.15
    atomic_volume = 1.1664e-29
    lattice_parameter = 3.6e-10
    D0 = 6e-05
    Em_eV = 1.35
    Z = 1.0
    Z_loop = 33.0
    alpha = 0.65
    formation_energy_eV = 1.7
    energy_exponent = 0.8
    energy_coefficient_eV = 1.619212219414472
    energy_constant_eV = 0.08078778058552794
    num_cluster_sizes = 2000
    switch_size = 20
    opposite_coupled_var = interstitial_clusters
    opposite_D0 = 1e-07
    opposite_Em_eV = 0.45
    opposite_Z = 1.1
    opposite_Z_loop = 63.0
    opposite_alpha = 0.8
    partner_Z_loop = 33.0
    partner_alpha = 0.65
  []
  [vacancy_cascade]
    type = CascadeSpectrumNodalKernel
    variable = vacancy_clusters
    generation = 0.3
    f2 = 0.06
    f3 = 0.03
    f4 = 0.02
  []
  [vacancy_recovery]
    type = DislocationRecoveryNodalKernel
    variable = vacancy_clusters
    D0 = 6e-05
    Em_eV = 1.35
    temperature = 603.15
    Z = 1.0
    rho0 = 1e10
    K = 16342508.0
    burgers_vector = 2.0784609690826528e-10
  []
  [vacancy_recomb]
    type = PointDefectRecombinationNodalKernel
    variable = vacancy_clusters
    coupled_var = interstitial_clusters
    recombination_rate = 13103131784.220617
  []
  [interstitial_dot]
    type = ArrayTimeDerivativeNodalKernel
    variable = interstitial_clusters
  []
  [interstitial_dyn]
    type = ContinuumLoopNodalKernel
    variable = interstitial_clusters
    generation = 0.072
    temperature = 603.15
    atomic_volume = 1.1664e-29
    lattice_parameter = 3.6e-10
    D0 = 1e-07
    Em_eV = 0.45
    Z = 1.1
    Z_loop = 63.0
    alpha = 0.8
    formation_energy_eV = 4.1
    energy_exponent = 0.8
    energy_coefficient_eV = 4.722702306625544
    energy_constant_eV = -0.6227023066255439
    num_cluster_sizes = 2000
    switch_size = 20
    opposite_coupled_var = vacancy_clusters
    opposite_D0 = 6e-05
    opposite_Em_eV = 1.35
    opposite_Z = 1.0
    opposite_Z_loop = 33.0
    opposite_alpha = 0.65
    partner_Z_loop = 63.0
    partner_alpha = 0.8
  []
  [interstitial_cascade]
    type = CascadeSpectrumNodalKernel
    variable = interstitial_clusters
    generation = 0.3
    f2 = 0.5
    f3 = 0.2
    f4 = 0.06
  []
  [interstitial_recovery]
    type = DislocationRecoveryNodalKernel
    variable = interstitial_clusters
    D0 = 1e-07
    Em_eV = 0.45
    temperature = 603.15
    Z = 1.1
    rho0 = 1e10
    K = 16342508.0
    burgers_vector = 2.0784609690826528e-10
  []
  [interstitial_recomb]
    type = PointDefectRecombinationNodalKernel
    variable = interstitial_clusters
    coupled_var = vacancy_clusters
    recombination_rate = 13103131784.220617
  []
[]
[Postprocessors]
  [vacancy_monomer_concentration]
    type = ClusterSizeConcentration
    clusters = vacancy_clusters
    n_size = 1
  []
  [vacancy_resolvable_density]
    type = ClusterContinuumMoment
    clusters = vacancy_clusters
    quantity = density
    num_cluster_sizes = 2000
    switch_size = 20
    n_minimum = 32
    atomic_volume = 1.1664e-29
    r1 = 1.3365274487142442e-10
    radius_exponent = 0.5
  []
  [vacancy_resolvable_avg_radius]
    type = ClusterContinuumMoment
    clusters = vacancy_clusters
    quantity = average_radius
    num_cluster_sizes = 2000
    switch_size = 20
    n_minimum = 32
    atomic_volume = 1.1664e-29
    r1 = 1.3365274487142442e-10
    radius_exponent = 0.5
  []
  [interstitial_monomer_concentration]
    type = ClusterSizeConcentration
    clusters = interstitial_clusters
    n_size = 1
  []
  [interstitial_resolvable_density]
    type = ClusterContinuumMoment
    clusters = interstitial_clusters
    quantity = density
    num_cluster_sizes = 2000
    switch_size = 20
    n_minimum = 32
    atomic_volume = 1.1664e-29
    r1 = 1.3365274487142442e-10
    radius_exponent = 0.5
  []
  [interstitial_resolvable_avg_radius]
    type = ClusterContinuumMoment
    clusters = interstitial_clusters
    quantity = average_radius
    num_cluster_sizes = 2000
    switch_size = 20
    n_minimum = 32
    atomic_volume = 1.1664e-29
    r1 = 1.3365274487142442e-10
    radius_exponent = 0.5
  []
[]
[Preconditioning]
  [smp]
    type = SMP
    full = false
    off_diag_row = 'vacancy_clusters_0 interstitial_clusters_0'
    off_diag_column = 'interstitial_clusters_0 vacancy_clusters_0'
  []
[]
[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_ordering_type'
  petsc_options_value = 'lu       nd'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-20
  dt = 0.01
  num_steps = 10
[]
[Outputs]
  csv = true
[]
