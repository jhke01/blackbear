# Phase 3c combined check: [Recovery] and cascade seeding active together, confirming the two
# modules compose without surprises. Combines defect_dynamics_recovery_action.i's [Recovery]
# block with defect_dynamics_cascade_action.i's per-species cascade_f2/f3/f4 values on the same
# two-species small-grid setup both files already use. The two modules modify different terms of
# the monomer-row equation (Recovery: sink; Cascade: generation) and touch different array
# components for their own dedicated kernels (Recovery: component 0, the monomer; Cascade:
# components 1-3, sizes 2-4), so no interaction is expected -- this test is the regression gate
# for that claim, including that it runs with zero warnings (dislocation_density is left unset,
# matching defect_dynamics_recovery_action.i's own convention, so the dislocation_density-vs-
# [Recovery] warning does not fire).

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

  [Grid]
    explicit_max = 6
    spacing = linear
    linear_width = 4
    geometric_ratio = 1.5
  []

  [Recovery]
    rho0 = 1e15
    K = 5e13
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
    cascade_f2 = 0.2
    cascade_f3 = 0.05
    cascade_f4 = 0.02
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
    cascade_f2 = 0.05
    cascade_f3 = 0.05
    cascade_f4 = 0.1
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
