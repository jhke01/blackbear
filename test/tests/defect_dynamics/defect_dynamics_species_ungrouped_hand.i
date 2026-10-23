# Hand-built twin of defect_dynamics_species_ungrouped_action.i: the loop arrays have use_grouping = false
# (20 sizes, one unknown each) wherever they are described, and the void arrays stay grouped.
# The run uses a high dose rate and a long step (dpa_rate 1.5, dt 1e-2, 6 steps) so that clusters well
# beyond the explicit sizes are populated. With the original 3 steps of 1e-5 s only the monomer and dimer
# were non-negligible, so changing a bin width or explicit_max left the output unchanged. Now a change of
# the void bin width moves the result by about 0.3% and of explicit_max by about 0.04%, against a
# comparison tolerance of 5.5e-6.

[Mesh]
  type = SingleNodeMesh
  position = '0 0 0'
[]

[Problem]
  use_hash_table_matrix_assembly = true
  restore_original_nonzero_pattern = false
[]

[GroupedVariables]
  [vacancy_clusters]
    family = LAGRANGE
    order = FIRST
    num_cluster_sizes = 40
    explicit_max = 6
    group_spacing = linear
    group_linear_width = 4
    group_geometric_ratio = 1.5
  []
  [interstitial_clusters]
    family = LAGRANGE
    order = FIRST
    num_cluster_sizes = 20
    explicit_max = 20
    group_spacing = linear
    group_linear_width = 4
    group_geometric_ratio = 1.5
  []
[]

[Functions]
  [vacancy_index_fn]
    type = PiecewiseConstant
    x = '0 0.9'
    y = '4.465175149565e-15 0'
    direction = left
  []
  [interstitial_index_fn]
    type = PiecewiseConstant
    x = '0 0.9'
    y = '3.473920870009918e-36 0'
    direction = left
  []
[]

[ICs]
  [vacancy_ic]
    type = GroupedArrayFunctionIndexIC
    variable = vacancy_clusters
    function = vacancy_index_fn
    num_cluster_sizes = 40
    use_grouping = true
    explicit_max = 6
    group_spacing = linear
    group_linear_width = 4
    group_geometric_ratio = 1.5
  []
  [interstitial_ic]
    type = GroupedArrayFunctionIndexIC
    variable = interstitial_clusters
    function = interstitial_index_fn
    num_cluster_sizes = 20
    use_grouping = false
    explicit_max = 20
    group_spacing = linear
    group_linear_width = 4
    group_geometric_ratio = 1.5
  []
[]

[NodalKernels]
  [vacancy_dot]
    type = ArrayTimeDerivativeNodalKernel
    variable = vacancy_clusters
  []
  [vacancy_dyn]
    type = VacancyVoidGroupedNodalKernel
    variable = vacancy_clusters

    generation = 0.15
    sink = 2.2521186324898173

    temperature = 573
    D0_v = 1.19073105525e-8
    Em_v_eV = 0.7644
    atomic_volume = 1.1664e-29
    lattice_parameter = 0.36e-9

    formation_energy_eV = 1.7
    void_energy_exponent = 0.67
    void_energy_coefficient_eV = 2.030206194432379
    void_energy_constant_eV = -0.33020619443237886

    D0_i = 8.781e-8
    Em_i_eV = 0.4684

    num_cluster_sizes = 40
    use_grouping = true
    explicit_max = 6
    group_spacing = linear
    group_linear_width = 4
    group_geometric_ratio = 1.5
    enforce_group_nonnegative = true
    opposite_coupled_var = interstitial_clusters
    opposite_z_v = 1.0
    opposite_z_loop_v = 33.0
    opposite_alpha_loop_v = 0.65
  
    opposite_num_cluster_sizes = 20
    opposite_use_grouping = false
    opposite_explicit_max = 20
    opposite_group_spacing = linear
    opposite_group_linear_width = 4
    opposite_group_geometric_ratio = 1.5
    opposite_enforce_group_nonnegative = true
  []

  [interstitial_dot]
    type = ArrayTimeDerivativeNodalKernel
    variable = interstitial_clusters
  []
  [interstitial_dyn]
    type = InterstitialLoopGroupedNodalKernel
    variable = interstitial_clusters

    generation = 0.15
    sink = 7330.943459224491

    temperature = 573
    D0_i = 8.781e-8
    Em_i_eV = 0.4684
    atomic_volume = 1.1664e-29
    lattice_parameter = 0.36e-9

    z_i = 1.10
    z_loop_i = 63.0
    alpha_loop_i = 0.8

    formation_energy_eV = 4.1
    loop_energy_exponent = 0.80
    loop_energy_coefficient_eV = 4.722702306625544
    loop_energy_constant_eV = -0.6227023066255439

    D0_v = 1.19073105525e-8
    Em_v_eV = 0.7644
    z_v = 1.0
    z_loop_v = 33.0
    alpha_loop_v = 0.65

    num_cluster_sizes = 20
    use_grouping = false
    explicit_max = 20
    group_spacing = linear
    group_linear_width = 4
    group_geometric_ratio = 1.5
    enforce_group_nonnegative = true
    opposite_coupled_var = vacancy_clusters
  
    opposite_num_cluster_sizes = 40
    opposite_use_grouping = true
    opposite_explicit_max = 6
    opposite_group_spacing = linear
    opposite_group_linear_width = 4
    opposite_group_geometric_ratio = 1.5
    opposite_enforce_group_nonnegative = true
  []
[]

[Postprocessors]
  [vacancy_monomer_concentration]
    type = ClusterSizeConcentrationGrouped
    clusters = vacancy_clusters
    n_size = 1
  
    num_cluster_sizes = 40
    use_grouping = true
    explicit_max = 6
    group_spacing = linear
    group_linear_width = 4
    group_geometric_ratio = 1.5
    enforce_group_nonnegative = true
  []
  [interstitial_monomer_concentration]
    type = ClusterSizeConcentrationGrouped
    clusters = interstitial_clusters
    n_size = 1
  
    num_cluster_sizes = 20
    use_grouping = false
    explicit_max = 20
    group_spacing = linear
    group_linear_width = 4
    group_geometric_ratio = 1.5
    enforce_group_nonnegative = true
  []
  [vacancy_resolvable_density]
    type = ClusterTotalDensityGrouped
    clusters = vacancy_clusters
    n_minimum = 2
    atomic_volume = 1.1664e-29
    num_cluster_sizes = 40
    use_grouping = true
    explicit_max = 6
    group_spacing = linear
    group_linear_width = 4
    group_geometric_ratio = 1.5
    enforce_group_nonnegative = true
  []
  [vacancy_resolvable_avg_radius]
    type = ClusterAverageRadiusGrouped
    clusters = vacancy_clusters
    r1 = 1.4068667552341914e-10
    radius_exponent = 0.3333333333333333
    n_minimum = 2
    num_cluster_sizes = 40
    use_grouping = true
    explicit_max = 6
    group_spacing = linear
    group_linear_width = 4
    group_geometric_ratio = 1.5
    enforce_group_nonnegative = true
  []
  [interstitial_resolvable_density]
    type = ClusterTotalDensityGrouped
    clusters = interstitial_clusters
    n_minimum = 2
    atomic_volume = 1.1664e-29
    num_cluster_sizes = 20
    use_grouping = false
    explicit_max = 20
    group_spacing = linear
    group_linear_width = 4
    group_geometric_ratio = 1.5
    enforce_group_nonnegative = true
  []
  [interstitial_resolvable_avg_radius]
    type = ClusterAverageRadiusGrouped
    clusters = interstitial_clusters
    r1 = 1.3365274487142442e-10
    radius_exponent = 0.5
    n_minimum = 2
    num_cluster_sizes = 20
    use_grouping = false
    explicit_max = 20
    group_spacing = linear
    group_linear_width = 4
    group_geometric_ratio = 1.5
    enforce_group_nonnegative = true
  []
[]

[Preconditioning]
  [smp]
    # See point_defect_recombination_small.i for why full=true is used here (SMP's
    # off_diag_row/off_diag_column do not resolve array-variable logical names). Affordable at
    # this test's small component count.
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_ordering_type'
  petsc_options_value = 'lu       nd'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-20
  dt = 1e-2
  num_steps = 6
[]

[Outputs]
  csv = true
[]
