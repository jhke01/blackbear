# [DefectDynamics] action equivalent of defect_dynamics_opposite_absorption_small_action.i, but
# with opposite_absorption = false: the two species run side by side with no coupling between
# them at all (no opposite-species absorption, and no recombination since recombination_radius
# is also left at zero). Each species therefore reduces exactly to the self-reaction-only
# behavior the low-level kernels already have when opposite_coupled_var is left unset -- the
# same behavior a single-species-only [DefectDynamics] block gets, just run for two species in
# one input. Verified directly against vacancy_void_grouped_small_action.i/
# interstitial_loop_grouped_small_action.i run with this same dislocation_density: bit-identical
# monomer concentrations.

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
  opposite_absorption = false

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
  dt = 1e-5
  num_steps = 3
[]

[Outputs]
  csv = true
[]
