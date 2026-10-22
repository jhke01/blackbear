# Defect conservation with two loop species. Each absorption event of one species' monomer on the
# other species' loop removes one defect of EACH type, so with no generation, no sinks and no
# recombination the difference (total vacancy atoms - total interstitial atoms) must stay constant.
# That only holds if both sides of one event use one rate: the interstitial-loop kernel shrinks a
# loop with its own vacancy-capture coefficient (Piece A), and the vacancy-loop kernel removes the
# vacancy monomer with opposite_geometry = loop (Piece B). With the void formula in Piece B the
# difference drifts by tens of percent of the held atoms; with opposite_geometry = loop it stays
# constant to about 1e-8. The held_difference column of the gold file shows the constant.

[Mesh]
  type = SingleNodeMesh
  position = '0 0 0'
[]
[Problem]
  use_hash_table_matrix_assembly = true
  restore_original_nonzero_pattern = false
[]
[GlobalParams]
  num_cluster_sizes = 30
  use_grouping = true
  explicit_max = 30
  group_spacing = linear
  group_linear_width = 4
  group_geometric_ratio = 1.5
  enforce_group_nonnegative = true
  group_nonnegative_tolerance_factor = 0.0
  opposite_num_cluster_sizes = 30
  opposite_use_grouping = true
  opposite_explicit_max = 30
  opposite_group_spacing = linear
  opposite_group_linear_width = 4
  opposite_group_geometric_ratio = 1.5
  opposite_enforce_group_nonnegative = true
  opposite_group_nonnegative_tolerance_factor = 0.0
[]
[GroupedVariables]
  [vacancy_clusters]
    family = LAGRANGE
    order = FIRST
  []
  [interstitial_clusters]
    family = LAGRANGE
    order = FIRST
  []
[]
[Functions]
  [vacancy_index_fn]
    type = ParsedFunction
    expression = 'if(x < 15.5, 1.5e-9 * exp(-0.2 * x), 1e-11)'
  []
  [interstitial_index_fn]
    type = ParsedFunction
    expression = 'if(x < 15.5, 1.2e-9 * exp(-0.15 * x), 1e-11)'
  []
[]
[ICs]
  [v_ic]
    type = GroupedArrayFunctionIndexIC
    variable = vacancy_clusters
    function = vacancy_index_fn
  []
  [i_ic]
    type = GroupedArrayFunctionIndexIC
    variable = interstitial_clusters
    function = interstitial_index_fn
  []
[]
[NodalKernels]
  [vac_dot]
    type = ArrayTimeDerivativeNodalKernel
    variable = vacancy_clusters
  []
  [vac_dyn]
    type = InterstitialLoopGroupedNodalKernel
    variable = vacancy_clusters
    generation = 0.0
    sink = 0.0
    temperature = 573
    D0_i = 1.19073105525e-08
    Em_i_eV = 0.7644
    atomic_volume = 1.1664e-29
    lattice_parameter = 3.6e-10
    z_i = 1.0
    z_loop_i = 33.0
    alpha_loop_i = 0.65
    formation_energy_eV = 1.7
    loop_energy_exponent = 0.65
    loop_energy_coefficient_eV = 2.1083398701275446
    loop_energy_constant_eV = -0.4083398701275447
    D0_v = 8.781e-08
    Em_v_eV = 0.4684
    z_v = 1.1
    z_loop_v = 63.0
    alpha_loop_v = 0.8
    opposite_coupled_var = interstitial_clusters
    opposite_geometry = loop
    opposite_z_i = 1.0
    opposite_z_loop_i = 33.0
    opposite_alpha_loop_i = 0.65
  []
  [int_dot]
    type = ArrayTimeDerivativeNodalKernel
    variable = interstitial_clusters
  []
  [int_dyn]
    type = InterstitialLoopGroupedNodalKernel
    variable = interstitial_clusters
    generation = 0.0
    sink = 0.0
    temperature = 573
    D0_i = 8.781e-08
    Em_i_eV = 0.4684
    atomic_volume = 1.1664e-29
    lattice_parameter = 3.6e-10
    z_i = 1.1
    z_loop_i = 63.0
    alpha_loop_i = 0.8
    formation_energy_eV = 4.1
    loop_energy_exponent = 0.8
    loop_energy_coefficient_eV = 4.722702306625544
    loop_energy_constant_eV = -0.6227023066255439
    D0_v = 1.19073105525e-08
    Em_v_eV = 0.7644
    z_v = 1.0
    z_loop_v = 33.0
    alpha_loop_v = 0.65
    opposite_coupled_var = vacancy_clusters
    opposite_geometry = loop
    opposite_z_i = 1.1
    opposite_z_loop_i = 63.0
    opposite_alpha_loop_i = 0.8
  []
[]
[Postprocessors]
  [held_vacancy]
    type = ClusterTotalContentGrouped
    clusters = vacancy_clusters
    num_cluster_sizes = 30
    explicit_max = 30
    group_spacing = linear
    group_linear_width = 4
    group_geometric_ratio = 1.5
  []
  [held_interstitial]
    type = ClusterTotalContentGrouped
    clusters = interstitial_clusters
    num_cluster_sizes = 30
    explicit_max = 30
    group_spacing = linear
    group_linear_width = 4
    group_geometric_ratio = 1.5
  []
  [held_difference]
    type = ParsedPostprocessor
    pp_names = 'held_vacancy held_interstitial'
    expression = 'held_vacancy - held_interstitial'
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
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-20
  dt = 1e-4
  num_steps = 20
[]
[Outputs]
  csv = true
[]
