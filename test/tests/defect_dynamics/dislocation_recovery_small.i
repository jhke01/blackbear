# Direct numerical check of DislocationRecoveryNodalKernel's closed-form rho(t), isolated from
# the void self-reaction kernel's own capillary emission (generation, void_energy_* set to make
# that kernel's own contribution zero for n > 1: no monomer production, same formation_energy_eV
# at every size via a zero exponent/coefficient would be ambiguous, so instead generation = 0 and
# a short enough run that absorption/emission from the near-zero initial monomer population
# stays negligible next to the recovery sink's own loss) so the only appreciable loss mechanism
# left is this kernel's own rho(t)*D*Z term, letting the resulting monomer concentration be
# compared directly against the closed-form rho(t) evaluated by hand with the same
# rho0/K/burgers_vector (see the comment in GenericDislocationRecoveryNodalKernel.C: rho(t) =
# (rho0^(-1/2) + K*b^2*t/2)^(-2)).
#
# D0/Em_eV/temperature/Z/rho0/K/burgers_vector below exactly match dislocation_recovery_jacobian.i,
# so this test's hand-computed rho(t) values (t=1e-5: 9.997840e13 /m^3, t=2e-5: 9.995681e13,
# t=3e-5: 9.993523e13) are directly reusable to sanity-check that Jacobian test's own setup.

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
    y = '4.465175149565e-15 0'
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
    void_energy_coefficient_eV = 2.030206194432379
    void_energy_constant_eV = -0.33020619443237886
  []
  [vacancy_recovery]
    type = DislocationRecoveryNodalKernel
    variable = vacancy_clusters

    D0 = 1.19073105525e-8
    Em_eV = 0.7644
    temperature = 573
    Z = 1.0

    rho0 = 1e14
    K = 5e13
    burgers_vector = 2.0784609690826528e-10
  []
[]

[Postprocessors]
  [vacancy_monomer_concentration]
    type = ClusterSizeConcentrationGrouped
    clusters = vacancy_clusters
    n_size = 1
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
