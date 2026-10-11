# Confirms CascadeSpectrumNodalKernel conserves total point-defect production: with the full
# capillary absorption/emission active (unlike cascade_spectrum_small.i's isolated setup), the
# weighted sum sum_n(n*C_n) must still grow at exactly generation = 1.5e-4 per second, since
# self-absorption/emission only move point defects between adjacent cluster sizes (a
# conservative birth-death chain) and never create or destroy them -- the ONLY net point-defect
# source in this system is the split generation term (monomer share + cascade shares), each
# divided by its own cluster size n so that n*(generation*f_n/n) = generation*f_n contributes
# exactly its own fraction of the total rate, regardless of how f2/f3/f4 split it (Pokor's own
# literal G(n) = generation*f_n form, by contrast, would double/triple/quadruple-count the
# point defects seeded at n = 2, 3, 4). generation*f2/f3/f4 are the same values
# defect_dynamics_cascade_action.i uses for the vacancy species.

[Mesh]
  type = SingleNodeMesh
  position = '0 0 0'
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
  [vacancy_clusters]
    family = LAGRANGE
    order = FIRST
  []
[]

[Functions]
  [index_fn]
    type = PiecewiseConstant
    x = '0 0.9'
    y = '0 0'
    direction = left
  []
[]

[ICs]
  [vacancy_ic]
    type = GroupedArrayFunctionIndexIC
    variable = vacancy_clusters
    function = index_fn
  []
[]

[NodalKernels]
  [vacancy_dot]
    type = ArrayTimeDerivativeNodalKernel
    variable = vacancy_clusters
  []
  [vacancy_dyn]
    type = VacancyVoidGroupedNodalKernel
    variable = vacancy_clusters

    # generation*(1 - f2 - f3 - f4): the monomer's own share after the cascade kernel below seeds
    # the rest directly, matching what DefectDynamicsAction computes for a species with these same
    # cascade_f2/f3/f4 values when a [Cascade] sub-block is present.
    generation = 9.75e-05
    sink = 0.0

    temperature = 573
    D0_v = 1.19073105525e-8
    Em_v_eV = 0.7644
    atomic_volume = 1.1664e-29
    lattice_parameter = 0.36e-9

    formation_energy_eV = 1.7
    void_energy_exponent = 0.67
    void_energy_coefficient_eV = 2.030206194432379
    void_energy_constant_eV = -0.33020619443237886
  []
  [vacancy_cascade]
    type = CascadeSpectrumNodalKernel
    variable = vacancy_clusters

    generation = 1.5e-4
    f2 = 0.2
    f3 = 0.1
    f4 = 0.05
  []
[]

[Postprocessors]
  [total_point_defects]
    type = ParsedPostprocessor
    pp_names = 'c1 c2 c3 c4'
    expression = 'c1 + 2*c2 + 3*c3 + 4*c4'
  []
  [c1]
    type = ClusterSizeConcentrationGrouped
    clusters = vacancy_clusters
    n_size = 1
  []
  [c2]
    type = ClusterSizeConcentrationGrouped
    clusters = vacancy_clusters
    n_size = 2
  []
  [c3]
    type = ClusterSizeConcentrationGrouped
    clusters = vacancy_clusters
    n_size = 3
  []
  [c4]
    type = ClusterSizeConcentrationGrouped
    clusters = vacancy_clusters
    n_size = 4
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
  dt = 1e-5
  num_steps = 3
[]

[Outputs]
  csv = true
[]
