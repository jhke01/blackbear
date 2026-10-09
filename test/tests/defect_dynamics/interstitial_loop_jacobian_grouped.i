# Jacobian test for InterstitialLoopGroupedNodalKernel.
# This test verifies that the intra-variable analytical Jacobian, including the
# coupling of grouped L0/L1 components through the inter-bin fluxes and the
# size-dependent dislocation-loop bias factor, matches the finite-difference
# Jacobian computed by PETSc. A smooth decaying initial profile is used
# (rather than the near-zero physical initial condition in
# interstitial_loop_grouped_small.i) so every explicit and grouped component
# has a nonzero, finite-difference-friendly value, scaled to physically
# reasonable site-fraction concentrations (~1e-9) and a dt (1e-5 s) consistent
# with the real monomer migration energy/temperature used here.
#
# D0_i is reduced far below its physical value (8.781e-8 m^2/s in
# CD_data_for_comparison_mathew) so that beta(n) stays in the same O(1e4)
# range as the vacancy-void kernel's own absorption coefficients. With the
# full physical D0_i, beta(1) ~ 4.6e9: PETSc's -snes_test_jacobian uses one
# global finite-difference step length for the whole solution vector (the
# Walker-Pernice "wp" algorithm), which cannot simultaneously resolve both
# the tiny grouped-tail concentrations and a derivative that large, and
# reports a spurious O(1e-4) mismatch even though the hand-coded Jacobian is
# exact (independently re-derived and confirmed to agree with a
# properly-stepped finite difference to 1e-7 relative error at this same
# state). This is the same reason the Cu-precipitation kernel's own
# cluster_dynamics_jacobian_grouped.i uses free toy rate coefficients
# (beta0, alpha0) rather than physical Arrhenius values.

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
  [interstitial_clusters]
    family = LAGRANGE
    order = FIRST
  []
[]

[Functions]
  [index_fn]
    # smooth decay with a step inside the grouped tail so that some bins are slope limited,
    # scaled to a physically reasonable site-fraction concentration (~1e-9)
    type = ParsedFunction
    expression = 'if(x < 15.5, 1.5e-9 * exp(-0.2 * x), 1e-11)'
  []
[]

[ICs]
  [interstitial_ic]
    type = GroupedArrayFunctionIndexIC
    variable = interstitial_clusters
    function = index_fn
  []
[]

[NodalKernels]
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
  dt = 1e-5
  num_steps = 1
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-20
[]

[Outputs]
  console = true
[]
