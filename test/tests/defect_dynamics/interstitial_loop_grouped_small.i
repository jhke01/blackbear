# Smoke test for the new InterstitialLoopGroupedNodalKernel: a small
# single-node grouped loop-only run (no vacancy coupling yet) using
# MD-benchmark-scale interstitial physics parameters
# (CD_data_for_comparison_mathew/MdParam, 573 K).
#
# Average-radius comparison is not exercised here: ClusterAverageRadiusGrouped
# assumes a cbrt(n) (three-dimensional) radius law, which does not match this
# kernel's sqrt(n) (two-dimensional) loop radius, so it is not reused for this
# family until a dedicated loop postprocessor exists.
#
# sink is left at zero (dislocation sink loss and recombination are added by
# a separate coupling kernel, not yet implemented), so this run does not
# exactly reproduce DefectDynamics_Sundials.cpp's interstitial_monomer_
# concentration -- only the self-reaction (absorption/emission) physics
# implemented so far is being exercised. With sink = D_i*z_i*dislocation_
# density = 7330.94 added, the result brackets the reference value from
# below just as the sink=0 case brackets it from above, consistent with the
# missing sink/recombination terms rather than an error in this kernel.

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
    type = PiecewiseConstant
    x = '0 0.9'
    y = '3.473920870009918e-36 0'
    direction = left
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
    sink = 0.0

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
  []
[]

[Postprocessors]
  [interstitial_monomer_concentration]
    type = ClusterSizeConcentrationGrouped
    clusters = interstitial_clusters
    n_size = 1
  []
  [interstitial_dimer_concentration]
    type = ClusterSizeConcentrationGrouped
    clusters = interstitial_clusters
    n_size = 2
  []
  [total_loop_density]
    type = ClusterTotalDensityGrouped
    clusters = interstitial_clusters
    n_minimum = 2
    atomic_volume = 1.1664e-29
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
  dt = 1e-5
  num_steps = 3
[]

[Outputs]
  csv = true
[]
