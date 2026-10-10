# [DefectDynamics] action equivalent of vacancy_void_grouped_small.i: same physics (a single
# void species, no interstitial coupling), but the array variable, Function/IC pair, time
# derivative, and VacancyVoidGroupedNodalKernel's capillary coefficients/initial condition are all
# generated from the physical parameters below instead of hand-typed. The monomer concentration
# this produces is bit-identical to vacancy_void_grouped_small_out.csv's (verified directly); the
# TEM-visible postprocessors below exercise the action's radius-to-n_minimum derivation, which
# vacancy_void_grouped_small.i's own ClusterAverageRadiusGrouped/ClusterTotalDensityGrouped
# postprocessors (hand-given r1/n_minimum=2) do not test.

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
  tem_visible_radius = 1.6e-10

  [Grid]
    num_cluster_sizes = 30
    explicit_max = 6
    spacing = linear
    linear_width = 4
    geometric_ratio = 1.5
  []

  [vacancy]
    geometry = void
    D0 = 1.19073105525e-8
    Em = 0.7644
    Ef = 1.7
    Eb2 = 0.5
    exponent = 0.67
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
