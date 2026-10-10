# PetscJacobianTester for the new DislocationRecoveryNodalKernel: a single void species, no
# opposite-species coupling, with the core kernel's own constant sink left at zero (sink = 0.0)
# so DislocationRecoveryNodalKernel is the ONLY sink contribution being exercised. Much simpler
# than point_defect_recombination_jacobian.i -- no second variable, no opposite_coupled_var, so
# full = false SMP suffices (there is no array-to-array coupling to resolve).
#
# start_time/dt are chosen so the Jacobian is evaluated at a nonzero simulation time (not t=0):
# rho(t) depends only on _t, and at t=0 the recovery sink's value happens to equal rho0*D*Z
# regardless of whether _t is actually read correctly, so evaluating only at t=0 would not catch
# a kernel that silently ignored _t. K is chosen (not a literal copy of Pokor's own quoted
# values, whose units are inconsistent with their own governing equation -- see
# GenericDislocationRecoveryNodalKernel's own class description) so rho(t) changes measurably
# over the few residual/Jacobian evaluations PetscJacobianTester performs.

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

    generation = 1.5e-4
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
  start_time = 1e-5
  dt = 1e-5
  num_steps = 3
[]

[Outputs]
  csv = true
[]
