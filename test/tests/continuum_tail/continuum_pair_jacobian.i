# Jacobian test for ContinuumLoopNodalKernel with two coupled families (vacancy and interstitial loops),
# each discrete up to size 20 and a Fokker-Planck continuum above it (132 components, largest size 2000),
# with the minmod limiter on. It covers the shrinking of each family by the other's monomer and the loss
# of both monomers, so the off-diagonal blocks (the partner's monomer in every flux, the partner's
# clusters in the monomer row) are compared with a finite difference. No time-derivative kernel is
# present: its 1/dt term would hide an error in the other entries.

[Mesh]
  type = SingleNodeMesh
  position = '1 1 1'
[]
[Problem]
  use_hash_table_matrix_assembly = true
  restore_original_nonzero_pattern = false
[]
[Variables]
  [v]
    family = LAGRANGE
    order = FIRST
    components = 132
    array = true
  []
  [i]
    family = LAGRANGE
    order = FIRST
    components = 132
    array = true
  []
[]
[Functions]
  [v_fn]
    type = ParsedFunction
    expression = '1 * exp(-0.05 * x) * (0.6 + 0.4 * sin(0.9 * x))'
  []
  [i_fn]
    type = ParsedFunction
    expression = '2 * exp(-0.04 * x) * (0.6 + 0.4 * cos(0.7 * x))'
  []
[]
[ICs]
  [v_ic]
    type = ArrayFunctionIndexIC
    variable = v
    function = v_fn
  []
  [i_ic]
    type = ArrayFunctionIndexIC
    variable = i
    function = i_fn
  []
[]
[NodalKernels]
  [vacancy_dyn]
    type = ContinuumLoopNodalKernel
    variable = v
    generation = 1.1e-4
    sink = 0.37
    temperature = 603.15
    atomic_volume = 1.1664e-29
    lattice_parameter = 0.36e-9
    D0 = 6e-05
    Em_eV = 1.35
    Z = 1.0
    Z_loop = 33.0
    alpha = 0.65
    formation_energy_eV = 1.7
    energy_exponent = 0.8
    energy_coefficient_eV = 1.619212219414472
    energy_constant_eV = 0.08078778058552794
    num_cluster_sizes = 2000
    switch_size = 20
    use_limiter = true
    opposite_coupled_var = i
    opposite_D0 = 1e-07
    opposite_Em_eV = 0.45
    opposite_Z = 1.1
    opposite_Z_loop = 63.0
    opposite_alpha = 0.8
  []
  [interstitial_dyn]
    type = ContinuumLoopNodalKernel
    variable = i
    generation = 1.1e-4
    sink = 0.37
    temperature = 603.15
    atomic_volume = 1.1664e-29
    lattice_parameter = 0.36e-9
    D0 = 1e-07
    Em_eV = 0.45
    Z = 1.1
    Z_loop = 63.0
    alpha = 0.8
    formation_energy_eV = 4.1
    energy_exponent = 0.8
    energy_coefficient_eV = 4.722702306625544
    energy_constant_eV = -0.6227023066255439
    num_cluster_sizes = 2000
    switch_size = 20
    use_limiter = true
    opposite_coupled_var = v
    opposite_D0 = 6e-05
    opposite_Em_eV = 1.35
    opposite_Z = 1.0
    opposite_Z_loop = 33.0
    opposite_alpha = 0.65
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
  dt = 1
  num_steps = 1
[]
[Outputs]
  csv = true
[]
