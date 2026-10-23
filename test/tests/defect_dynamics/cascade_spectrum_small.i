# Direct numerical check of CascadeSpectrumNodalKernel's constant source terms, isolated from the
# void self-reaction kernel's own capillary absorption/emission. generation = 0 on the core
# kernel and zero initial condition everywhere keep every self-absorption term (which all
# multiply the monomer concentration) at exactly zero, since nothing but this kernel's own
# constant source feeds the monomer-free array. The remaining possible contamination is thermal
# emission FROM sizes 2/3/4 back toward the monomer -- suppressed to below floating-point
# significance here by using a tiny capillary coefficient (void_energy_coefficient_eV = 0.01, vs.
# the production value ~2) so the binding energy stays close to formation_energy_eV = 1.7 eV,
# making alpha(2) ~ 1e-10 (vs. a cascade-driven dc2/dt ~ 1.5e-5): fourteen orders of magnitude too
# small to matter. With that, each seeded component's concentration after n steps of dt is simply
# n*dt*generation*f_k/k for cluster size k -- compared directly against hand-computed values below
# (generation = 1.5e-4, f2 = 0.2, f3 = 0.1, f4 = 0.05, dt = 1e-2):
#   t=0.01: c2=1.5e-7,  c3=5e-8,   c4=1.875e-8
#   t=0.02: c2=3.0e-7,  c3=1.0e-7, c4=3.75e-8
#   t=0.03: c2=4.5e-7,  c3=1.5e-7, c4=5.625e-8
# (With the former step of 1e-5 these were 1000 times smaller, and c3 and c4 were below the CSV
# comparison floor of 1e-10, so only c2 was actually compared.)
# CascadeSpectrumNodalKernel does not touch the monomer row itself (the enclosing
# DefectDynamicsAction reduces the core kernel's own "generation" parameter instead, exercised by
# defect_dynamics_cascade_action.i), so the monomer concentration stays at (solver-converged)
# zero throughout this test.

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

    generation = 0.0
    sink = 0.0

    temperature = 573
    D0_v = 1.19073105525e-8
    Em_v_eV = 0.7644
    atomic_volume = 1.1664e-29
    lattice_parameter = 0.36e-9

    formation_energy_eV = 1.7
    void_energy_exponent = 0.67
    void_energy_coefficient_eV = 0.01
    void_energy_constant_eV = 1.69
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
  [vacancy_monomer]
    type = ClusterSizeConcentrationGrouped
    clusters = vacancy_clusters
    n_size = 1
  []
  [vacancy_c2]
    type = ClusterSizeConcentrationGrouped
    clusters = vacancy_clusters
    n_size = 2
  []
  [vacancy_c3]
    type = ClusterSizeConcentrationGrouped
    clusters = vacancy_clusters
    n_size = 3
  []
  [vacancy_c4]
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
  dt = 1e-2
  num_steps = 3
[]

[Outputs]
  csv = true
[]
