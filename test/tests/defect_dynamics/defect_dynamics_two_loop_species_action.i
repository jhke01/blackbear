# [DefectDynamics] action equivalent of a hand-built two-geometry=loop [NodalKernels] block (see
# blackbear_runs/cdDefectDynamics_pokorTable3_test/sa304_osiris_vacancy_as_loop.i for the full-scale
# validated precedent this reproduces the KERNEL WIRING of): vacancy modeled as geometry = loop
# (following Pokor et al. 2004's own Eqs 4a-4f, which treat vacancy and interstitial clusters with
# the same 2-D loop absorption form), coupled to an interstitial geometry = loop species through
# both opposite-species absorption (default on) and point-defect recombination
# (recombination_radius > 0). Confirms the [DefectDynamics] action can wire two loop species the
# same way it already wires one void + one loop, including correctly resolving each kernel's
# opposite_* parameters to "the other species" instead of silently defaulting to a null void
# species (the latent bug this change fixes in actAddNodalKernels()).

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
