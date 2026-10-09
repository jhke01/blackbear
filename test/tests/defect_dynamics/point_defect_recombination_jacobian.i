# Jacobian test for PointDefectRecombinationNodalKernel.
# This test verifies that the hand-coded on- and off-diagonal Jacobian entries
# for the recombination coupling term match the finite-difference Jacobian
# computed by PETSc. Uses the same smooth, well-conditioned decaying initial
# profile and reduced test-only D0_i as interstitial_loop_jacobian_grouped.i,
# for the same reason given there (the physical D0_i makes beta(1) ~4.6e9,
# far outside what PETSc's -snes_test_jacobian can resolve with its single
# global finite-difference step length).

[Mesh]
  type = SingleNodeMesh
  position = '1 1 1'
[]

[Problem]
  use_hash_table_matrix_assembly = true
  restore_original_nonzero_pattern = false
[]

[GlobalParams]
  num_cluster_sizes = 30
  use_grouping = true
  explicit_max = 6
  group_spacing = linear
  group_linear_width = 4
  group_geometric_ratio = 1.5
  enforce_group_nonnegative = true
  group_nonnegative_tolerance_factor = 0.0
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
  [vacancy_ic]
    type = GroupedArrayFunctionIndexIC
    variable = vacancy_clusters
    function = vacancy_index_fn
  []
  [interstitial_ic]
    type = GroupedArrayFunctionIndexIC
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
    type = VacancyVoidGroupedNodalKernel
    variable = vacancy_clusters

    generation = 1.5e-4
    sink = 0.07330943459225948

    temperature = 573
    D0_v = 1.19073105525e-13
    Em_v_eV = 0.7644
    atomic_volume = 1.1664e-29
    lattice_parameter = 0.36e-9

    formation_energy_eV = 1.7
    void_energy_exponent = 0.67
    void_energy_coefficient_eV = 2.030206194432379
    void_energy_constant_eV = -0.33020619443237886
  []
  [vacancy_recombination]
    type = PointDefectRecombinationNodalKernel
    variable = vacancy_clusters
    coupled_var = interstitial_clusters
    recombination_rate = 64642.59505799026
  []

  [interstitial_dot]
    type = ArrayTimeDerivativeNodalKernel
    variable = interstitial_clusters
  []
  [interstitial_dyn]
    type = InterstitialLoopGroupedNodalKernel
    variable = interstitial_clusters

    generation = 1.5e-4
    sink = 0.07330943459225948

    temperature = 573
    D0_i = 8.781e-13
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
  []
  [interstitial_recombination]
    type = PointDefectRecombinationNodalKernel
    variable = interstitial_clusters
    coupled_var = vacancy_clusters
    recombination_rate = 64642.59505799026
  []
[]

[Preconditioning]
  [smp]
    # See point_defect_recombination_small.i for why full=true is used here (SMP's
    # off_diag_row/off_diag_column and coupled_groups do not resolve array-variable logical
    # names). Affordable at this test's small component count.
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_ordering_type'
  petsc_options_value = 'lu       nd'
  dt = 1e-5
  num_steps = 1
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-20
[]

[Outputs]
  console = true
[]
