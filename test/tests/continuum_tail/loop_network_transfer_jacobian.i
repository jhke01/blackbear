# Jacobian test for LoopNetworkTransferNodalKernel: loops of every size from 3 up (discrete sizes
# and continuum cells) leave at the rate k = coefficient * c_1 * sqrt(rho), which depends on the
# monomer (component 0) of the same variable, so the Jacobian has a column in component 0 for every row
# and a diagonal. PETSc compares it with a finite difference. The profile is rough and the monomer is
# large enough for the finite-difference step to resolve the dependence of k on it. No time-derivative
# kernel is present: its 1/dt term would hide an error in the other entries.

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
[Postprocessors]
  [rho]
    type = ConstantPostprocessor
    value = 3.0e14
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
[NodalKernels]
  [transfer]
    type = LoopNetworkTransferNodalKernel
    variable = u
    dislocation_density_postprocessor = rho
    D0 = 1e-7
    Em_eV = 0.45
    temperature = 603.15
    burgers_vector = 2.078460969082653e-10
    bias = 0.1
    ca_crho = 2.0
    poisson_ratio = 0.3
    minimum_size = 3
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
