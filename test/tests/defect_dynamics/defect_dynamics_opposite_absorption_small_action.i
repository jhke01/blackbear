# [DefectDynamics] action equivalent of defect_dynamics_opposite_absorption_small.i: a void and a
# loop species coupled only through opposite-species absorption (no recombination_radius given,
# so no PointDefectRecombinationNodalKernel is added -- matching the hand-written test's
# deliberate omission, see its own header comment). Every "opposite_*"/mirrored parameter the two
# low-level kernels need (D0_i/Em_i_eV on the void, D0_v/Em_v_eV/z_v/z_loop_v/alpha_loop_v on the
# loop, both grids) is derived from the two species blocks below instead of hand-typed.

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

  [vacancy]
    geometry = void
    D0 = 1.19073105525e-8
    Em = 0.7644
    Ef = 1.7
    Eb2 = 0.5
    exponent = 0.67
    Z = 1.0
    num_cluster_sizes = 30
    explicit_max = 6
    group_spacing = linear
    group_linear_width = 4
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
    num_cluster_sizes = 30
    explicit_max = 6
    group_spacing = linear
    group_linear_width = 4
    group_geometric_ratio = 1.5
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
  dt = 1e-5
  num_steps = 3
[]

[Outputs]
  csv = true
[]
