# Two species with about 430 array components each, so about 860 in all, and the pair of their arrays
# named in the preconditioner. MOOSE then allocates local block storage for every pair of components of
# the whole system (about 1.1 GB by the estimate in ClusterDynamics::checkSolverSetup), and the kernels
# warn about it before the first solve. The test checks the warning; the run stops at the warning
# because the harness treats warnings as errors, so nothing is allocated.

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
  dpa_rate = 1.5
  cascade_efficiency = 0.1
  dislocation_density = 1e15

  [vacancy]
    geometry = void
    D0 = 1.19073105525e-8
    Em = 0.7644
    Ef = 1.7
    Eb2 = 0.5
    exponent = 0.67
    Z = 1.0
    num_cluster_sizes = 2000
    resolvable_radius = 1.772539e-10
    explicit_max = 400
    group_spacing = linear
    group_linear_width = 100
    group_geometric_ratio = 1.5
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
    num_cluster_sizes = 2000
    resolvable_radius = 1.890133e-10
    explicit_max = 400
    group_spacing = geometric
    group_linear_width = 100
    group_geometric_ratio = 1.4
  []
[]

[Preconditioning]
  [smp]
    # See point_defect_recombination_small.i for why full=true is used here (SMP's
    # off_diag_row/off_diag_column do not resolve array-variable logical names). Affordable at
    # this test's small component count.
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
  dt = 1e-2
  num_steps = 6
[]

[Outputs]
  csv = true
[]
