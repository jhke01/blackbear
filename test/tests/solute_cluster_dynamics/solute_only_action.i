# Smoke test for [SoluteClusterDynamics] with no [DefectDynamics] block present (the solute-only
# case, no point-defect coupling): the array variable, Function/IC pair, time derivative, and
# ClusterDynamicsGroupedNodalKernel instance are all generated from the physical parameters below
# instead of hand-typed via [GroupedVariables]/[NodalKernels]. Same rate_model = simple/beta0/
# alpha0/generation/sink physics and grouping layout as cluster_dynamics_grouped_small.i, but with
# this action's own IC convention (monomer = initial_concentration, every other size starts at
# zero) rather than that test's hand-typed index function, so this is its own gold rather than a
# bit-identical match to that file's.

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1
[]

[Problem]
  use_hash_table_matrix_assembly = true
  restore_original_nonzero_pattern = false
[]

[SoluteClusterDynamics]
  initial_concentration = 1
  generation = 0.03
  sink = 0.12
  beta0 = 0.14
  alpha0 = 0.05

  num_cluster_sizes = 50
  use_grouping = true
  explicit_max = 10
  group_spacing = linear
  group_linear_width = 5
  group_geometric_ratio = 1.5
  enforce_group_nonnegative = true
  group_nonnegative_tolerance_factor = 0.0
[]

[Preconditioning]
  [smp]
    type = SMP
    full = false
  []
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_factor_mat_ordering_type'
  petsc_options_value = 'ilu      rcm'
  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-12
  dt = 0.1
  num_steps = 5
[]

[Outputs]
  csv = true
[]
