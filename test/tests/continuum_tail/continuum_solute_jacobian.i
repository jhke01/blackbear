# Jacobian test for ContinuumClusterDynamicsNodalKernel (solute clustering, interfacial-energy model): discrete up to size 20, a
# Fokker-Planck continuum up to size 2000 (132 components), minmod limiter on, with generation and a sink so that the monomer row
# is exercised. No time-derivative kernel, whose 1/dt term would hide an error in the other entries.

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
    expression = '0.1 * exp(-0.05 * x) * (0.6 + 0.4 * sin(0.9 * x))'
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
    type = ContinuumClusterDynamicsNodalKernel
    variable = u
    generation = 1e-6
    sink = 1e-3
    rate_model = interfacial_energy
    diffusivity_model = arrhenius
    temperature = 563.15
    D0 = 6.3e-5
    Q_eV = 2.29
    sigma = 0.37
    atomic_volume = 1.1782924e-29
    Omega_kB_K = 6255.0
    DeltaS_kB = 0.866
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
