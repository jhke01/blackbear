# Exercises [DefectDynamics]'s output_size_distribution/size_distribution_times: the same small
# two-loop-species grid as defect_dynamics_two_loop_species_action.i (see that file for the physics
# background), with vacancy additionally requesting its full per-size concentration distribution
# at one snapshot time via a ClusterSizeDistributionGrouped vector postprocessor. Confirms the
# dedicated sync_times-only CSV output writes one row per physical cluster size (1 to
# num_cluster_sizes = 30), with explicit sizes (1-6) and grouped-bin-reconstructed sizes (7-30)
# both present, while interstitial (output_size_distribution left at its false default) gets no
# such output and the regular [Outputs] csv = true scalar history is unaffected.

[Mesh]
  type = SingleNodeMesh
  position = '0 0 0'
[]

[Problem]
  use_hash_table_matrix_assembly = true
  restore_original_nonzero_pattern = false
[]

[DefectDynamics]
  temperature = 573
  atomic_volume = 1.1664e-29
  lattice_parameter = 0.36e-9
  dpa_rate = 1.5e-3
  cascade_efficiency = 0.1
  dislocation_density = 1e15
  recombination_radius = 9e-10

  size_distribution_times = '2e-5'

  [Grid]
    explicit_max = 6
    spacing = linear
    linear_width = 4
    geometric_ratio = 1.5
  []

  [vacancy]
    geometry = loop
    D0 = 1.19073105525e-8
    Em = 0.7644
    Ef = 1.7
    Eb2 = 0.5
    exponent = 0.65
    Z = 1.0
    Z_loop = 33.0
    alpha = 0.65
    Z_loop_opposite = 63.0
    alpha_opposite = 0.8
    num_cluster_sizes = 30
    output_size_distribution = true
  []

  [interstitial]
    geometry = loop
    D0 = 8.781e-8
    Em = 0.4684
    Ef = 4.1
    Eb2 = 0.6
    exponent = 0.8
    Z = 1.1
    Z_loop = 63.0
    alpha = 0.8
    Z_loop_opposite = 33.0
    alpha_opposite = 0.65
    num_cluster_sizes = 30
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
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-20
  dt = 1e-5
  num_steps = 3
[]

[Outputs]
  csv = true
[]
