# Jacobian test for the new opposite-species absorption terms in
# VacancyVoidGroupedNodalKernel and InterstitialLoopGroupedNodalKernel. This
# verifies that the hand-coded on- and off-diagonal Jacobian entries for the
# opposite-absorption coupling (an interstitial monomer shrinking a void, a
# vacancy monomer shrinking a loop, and the matching monomer-pool drain terms)
# match the finite-difference Jacobian computed by PETSc. Uses the same
# smooth, well-conditioned decaying initial profile and reduced test-only
# D0_v/D0_i as point_defect_recombination_jacobian.i, for the same reason
# given there (the physical diffusivities make beta(1)/beta_opposite(1)
# ~1e9-scale, far outside what PETSc's -snes_test_jacobian can resolve with
# its single global finite-difference step length). This was independently
# cross-checked against a from-scratch Python re-derivation of the exact
# residual/Jacobian before trusting this PETSc-reported ratio.
#
# The Jacobian is compared WITHOUT a time-derivative kernel and with the nonnegativity limiter off. A
# 1/dt mass term makes ||J|| about 6e5 against a physics part of about 0.5, so the tolerance below would
# accept a missing or sign-flipped off-diagonal block; the limiter's kinks are not differentiable and
# would hide the rest. The profile amplitude is large enough that the physics terms are well above the
# finite-difference noise. With SMP full = false (no cross-species block) this input gives a ratio of
# 0.1 to 0.7, against about 2e-8 here.

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
  enforce_group_nonnegative = false
  group_nonnegative_tolerance_factor = 0.0

  opposite_num_cluster_sizes = 30
  opposite_use_grouping = true
  opposite_explicit_max = 6
  opposite_group_spacing = linear
  opposite_group_linear_width = 4
  opposite_group_geometric_ratio = 1.5
  opposite_enforce_group_nonnegative = false
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
    expression = 'if(x < 15.5, 1 * exp(-0.2 * x), 1e-2)'
  []
  [interstitial_index_fn]
    type = ParsedFunction
    expression = 'if(x < 15.5, 1 * exp(-0.15 * x), 1e-2)'
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

    D0_i = 8.781e-13
    Em_i_eV = 0.4684

    opposite_coupled_var = interstitial_clusters
    opposite_z_v = 1.0
    opposite_z_loop_v = 33.0
    opposite_alpha_loop_v = 0.65
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

    D0_v = 1.19073105525e-13
    Em_v_eV = 0.7644
    z_v = 1.0
    z_loop_v = 33.0
    alpha_loop_v = 0.65

    opposite_coupled_var = vacancy_clusters
  []
[]

[Preconditioning]
  [smp]
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
