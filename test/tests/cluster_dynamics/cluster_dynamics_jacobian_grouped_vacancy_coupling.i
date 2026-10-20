# Jacobian test for the new off-diagonal block in ClusterDynamicsGroupedNodalKernel's
# vacancy_coupled_var coupling (GenericClusterDynamicsGroupedNodalKernelTempl<false>::
# computeQpOffDiagJacobian). vacancy_clusters is a plain (ungrouped) array variable -- only its
# component 0 is ever read, and ClusterGroupingLayout::monomerComponent() is always 0 regardless
# of grouping scheme, so no grouping layout is needed for the coupled variable itself. Only
# ArrayTimeDerivativeNodalKernel acts on it, so this is purely a Jacobian structure test. Its IC
# is set away from the thermal-equilibrium value (S_v = C_v(1)/C_v_thermal != 1) so the coupling
# derivative terms (betaPrime/alphaPrime) are genuinely nonzero and exercised by
# -snes_test_jacobian.

[Mesh]
  type = SingleNodeMesh
  position = '1 1 1'
[]

[Problem]
  use_hash_table_matrix_assembly = true
  restore_original_nonzero_pattern = false
[]

[GlobalParams]
  num_cluster_sizes = 30
  use_grouping = true
  explicit_max = 6
  group_spacing = linear
  group_linear_width = 4
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

[Variables]
  [vacancy_clusters]
    family = LAGRANGE
    order = FIRST
    components = 1
    array = true
  []
[]

[Functions]
  [index_fn]
    # smooth decay with a step inside the grouped tail so that some bins are slope limited
    type = ParsedFunction
    expression = 'if(x < 15.5, 1.5 * exp(-0.2 * x), 1e-6)'
  []
[]

[ICs]
  [clusters_ic]
    type = GroupedArrayFunctionIndexIC
    variable = clusters
    function = index_fn
  []
  [vacancy_ic]
    type = ArrayConstantIC
    variable = vacancy_clusters
    # 3x the thermal-equilibrium concentration at Ef = 1.6 eV, T = 600 K, so S_v != 1
    value = '4.362874814683014e-13'
  []
[]

[NodalKernels]
  [clusters_dot]
    type = ArrayTimeDerivativeNodalKernel
    variable = clusters
  []
  [clusters_dyn]
    type = ClusterDynamicsGroupedNodalKernel
    variable = clusters
    generation = 0.05
    sink = 0.15
    rate_model = interfacial_energy
    diffusivity_model = constant
    monomer_diffusivity = 1.0e-21
    temperature = 600
    atomic_volume = 1.1782924e-29
    binding_energy_model = formation_energy_table
    formation_energy_table_eV = '0 0.05 0.12 0.2 0.3 0.41 0.53 0.66 0.8 0.95 1.11 1.28 1.46 1.65
                                 1.85 2.06 2.28 2.51 2.75 3.0 3.26 3.53 3.81 4.1 4.4 4.71 5.03
                                 5.36 5.7 6.05'
    vacancy_coupled_var = vacancy_clusters
    vacancy_formation_energy_eV = 1.6
  []
  # No kernel besides the time derivative acts on vacancy_clusters.
  [vacancy_dot]
    type = ArrayTimeDerivativeNodalKernel
    variable = vacancy_clusters
  []
[]

[Preconditioning]
  [smp]
    # SMP's off_diag_row/off_diag_column and coupled_groups do not resolve array-variable logical
    # names (only their per-component names) -- a MOOSE-side limitation unrelated to this kernel.
    # Affordable at this test's small component count.
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_ordering_type'
  petsc_options_value = 'lu       nd'
  dt = 0.01
  num_steps = 1
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
[]

[Outputs]
  console = true
[]
