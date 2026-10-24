# [DefectDynamics] action with the cellular network dislocation density of Rest (2006) in [Recovery]
# (model = cellular_network): the network density rho is advanced at the end of each step from the
# interstitial monomer concentration, and both species' dislocation sink D*Z*rho uses it, read
# through dislocation_density_postprocessor. Compare defect_dynamics_recovery_action.i, which uses
# Pokor's closed-form recovery. The steps are long enough for rho to change from its initial value.

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

  # Rest (2006) cellular network dislocation density in place of Pokor's recovery law. The bias is a
  # round number, not a measured one.
  [Recovery]
    model = cellular_network
    rho0 = 1e13
    bias = 0.1
    ca_crho = 2.0
    poisson_ratio = 0.3
    grain_size = 4e-5
    interstitial_species = interstitial
    burgers_vector = 2.0784609690826528e-10
  []

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
  dt = 0.02
  num_steps = 5
[]

[Outputs]
  csv = true
[]
