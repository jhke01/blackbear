# [DefectDynamics] action equivalent of vacancy_void_grouped_small_action.i, but with the
# species' maximum cluster size given as max_cluster_radius_m instead of num_cluster_sizes. The
# radius below is exactly r_30 = r1*30^(1/3) for this input's own atomic_volume (void geometry:
# r1 = (3*V_at/(4*pi))^(1/3)), so the radius-to-n inversion (ceil((R/r1)^3)) resolves to exactly
# num_cluster_sizes = 30 -- the same value vacancy_void_grouped_small_action.i sets directly --
# making this run's output bit-identical to that test's, proving the two mutually exclusive ways
# to specify a species' maximum cluster size converge on the same physics.

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
  dislocation_density = 0

  [vacancy]
    geometry = void
    D0 = 1.19073105525e-8
    Em = 0.7644
    Ef = 1.7
    Eb2 = 0.5
    exponent = 0.67
    max_cluster_radius_m = 4.3714621134095103e-10
    resolvable_radius = 1.6e-10
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
