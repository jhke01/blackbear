# Jacobian test for ContinuumClusterDynamicsNodalKernel with vacancy_coupled_var: the monomer diffusivity is scaled by the vacancy
# monomer over its thermal value, so the off-diagonal column of the vacancy monomer is compared with a finite difference too.

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
  [vacancy]
    family = LAGRANGE
    order = FIRST
    components = 1
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
  [vacancy_ic]
    type = ArrayConstantIC
    variable = vacancy
    value = '2e-12'
  []
[]
[NodalKernels]
  [dyn]
    type = ContinuumClusterDynamicsNodalKernel
    variable = u
    generation = 1e-6
    sink = 1e-3
    vacancy_coupled_var = vacancy
    vacancy_formation_energy_eV = 1.6
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
  [vacancy_dot]
    type = ArrayTimeDerivativeNodalKernel
    variable = vacancy
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
