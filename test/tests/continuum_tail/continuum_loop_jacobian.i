# Jacobian test for ContinuumLoopNodalKernel: a loop family that is discrete up to size 20 and a
# Fokker-Planck continuum above it (132 components, largest size 2000), with the minmod limiter on.
# The hand-coded Jacobian differentiates through the limited slopes; PETSc compares it with a finite
# difference. The profile is rough (a sine on a decay) so that the limiter changes branch from cell to
# cell, and its amplitude is large enough for the finite-difference step to resolve the physics terms.
# No time-derivative kernel is present: its 1/dt term would hide an error in the other entries.

[Mesh]
  type = SingleNodeMesh
  position = '1 1 1'
[]
[Problem]
  use_hash_table_matrix_assembly = true
  restore_original_nonzero_pattern = false
[]
[Variables]
  [u]
    family = LAGRANGE
    order = FIRST
    components = 132
    array = true
  []
[]
[Functions]
  [ic_fn]
    type = ParsedFunction
    expression = '1 * exp(-0.05 * x) * (0.6 + 0.4 * sin(0.9 * x))'
  []
[]
[ICs]
  [u_ic]
    type = ArrayFunctionIndexIC
    variable = u
    function = ic_fn
  []
[]
[NodalKernels]
  [dyn]
    type = ContinuumLoopNodalKernel
    variable = u
    generation = 8.699999999999998e-08
    sink = 3.1468192651509726e-06
    temperature = 603.15
    atomic_volume = 1.1664e-29
    lattice_parameter = 0.36e-9
    burgers_vector = 2.078460969082653e-10
    D0 = 6e-5
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
  dt = 1
  num_steps = 1
[]
[Outputs]
  csv = true
[]
