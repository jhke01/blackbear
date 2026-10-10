[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1
[]

[Problem]
  # Required for the intra-variable Jacobian; keep the fixed sparsity after the first assembly
  use_hash_table_matrix_assembly = true
  restore_original_nonzero_pattern = false
[]

[GlobalParams]
  num_cluster_sizes = 50
  use_grouping = true
  explicit_max = 10
  group_spacing = linear
  group_linear_width = 5
  group_geometric_ratio = 1.5
  enforce_group_nonnegative = true
  group_nonnegative_tolerance_factor = 0.0
[]

[GroupedVariables]
  [clusters]
    family = LAGRANGE
    order = FIRST
  []
[]

[Functions]
  [index_fn]
    type = PiecewiseConstant
    x = '0 0.9'
    y = '1 0.02'
    direction = left
  []
[]

[ICs]
  [clusters_ic]
    type = GroupedArrayFunctionIndexIC
    variable = clusters
    function = index_fn
  []
[]

[NodalKernels]
  [clusters_dot]
    type = ArrayTimeDerivativeNodalKernel
    variable = clusters
  []
  [clusters_dyn]
    type = ADClusterDynamicsGroupedNodalKernel
    variable = clusters
    generation = 0.03
    sink = 0.12
    beta0 = 0.14
    alpha0 = 0.05
  []
[]

[Postprocessors]
  # Golubov et al. (2001) criterion (b): d/dt[sum_{n>=2} C_n] = J(1,t), the flux of monomers
  # combining into dimers. Golubov's N(t) counts clusters only (size >= 2), excluding free
  # monomers, so total_cluster_density below uses n_minimum = 2. For rate_model = simple,
  # beta(n) = beta0 * n^(1/3) and alpha(n) = alpha0 * n^(1/3) (see
  # GenericClusterDynamicsNodalKernel.md), so J(1,t) = beta(1)*C_1^2 - alpha(2)*C_2
  # = beta0*C_1^2 - alpha0*2^(1/3)*C_2. This is independent of this model's added generation/sink
  # terms, which act only on the monomer component and therefore drop out of a sum that excludes
  # it -- confirmed directly against GenericClusterDynamicsGroupedNodalKernel.C, where _generation
  # and _sink only appear in the monomer residual row.
  [total_cluster_density]
    type = ClusterTotalDensityGrouped
    clusters = clusters
    n_minimum = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [monomer_concentration]
    type = ClusterSizeConcentrationGrouped
    clusters = clusters
    n_size = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [dimer_concentration]
    type = ClusterSizeConcentrationGrouped
    clusters = clusters
    n_size = 2
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [density_rate]
    type = ChangeOverTimePostprocessor
    postprocessor = total_cluster_density
    divide_by_dt = true
    execute_on = 'TIMESTEP_END'
  []
  [expected_density_rate]
    type = ParsedPostprocessor
    pp_names = 'monomer_concentration dimer_concentration'
    expression = '0.14 * monomer_concentration^2 - 0.05 * 2^(1/3) * dimer_concentration'
    execute_on = 'TIMESTEP_END'
  []
  [density_rate_error]
    type = DifferencePostprocessor
    value1 = density_rate
    value2 = expected_density_rate
    execute_on = 'TIMESTEP_END'
  []

  # Golubov et al. (2001) criterion (c): d/dt[sum_n n*C_n] = J(1,t) + sum_n J(n,t), the total
  # point-defect content conservation law. Re-deriving this model's flux-form residuals directly
  # (rather than Golubov's unmodified equations, which have no generation/sink terms) shows every
  # internal absorption/emission flux cancels telescopically in this sum, leaving only the
  # monomer's own source and sink: d(total_cluster_content)/dt = G_1 - k_s*C_1.
  [total_cluster_content]
    type = ClusterTotalContentGrouped
    clusters = clusters
    n_minimum = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [content_rate]
    type = ChangeOverTimePostprocessor
    postprocessor = total_cluster_content
    divide_by_dt = true
    execute_on = 'TIMESTEP_END'
  []
  [expected_content_rate]
    type = ParsedPostprocessor
    pp_names = 'monomer_concentration'
    expression = '0.03 - 0.12 * monomer_concentration'
    execute_on = 'TIMESTEP_END'
  []
  [content_rate_error]
    type = DifferencePostprocessor
    value1 = content_rate
    value2 = expected_content_rate
    execute_on = 'TIMESTEP_END'
  []
[]

[Preconditioning]
  # Avoid the default full coupling matrix, which is O(N^2) in the number of array components
  [smp]
    type = SMP
    full = false
  []
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  # Reorder so that factoring the dense monomer row and column does not cost O(N^2)
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
