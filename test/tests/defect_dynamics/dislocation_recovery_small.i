# Check of DislocationRecoveryNodalKernel's time-dependent dislocation density
# rho(t) = (rho0^(-1/2) + K*b^2*t/2)^(-2) (see GenericDislocationRecoveryNodalKernel.C). The void
# kernel's own constant sink is zero, so the recovery kernel supplies the dislocation loss. A monomer
# source (generation) and a run of 2 s, long enough for rho(t) to fall by a factor of about 500, put
# the monomer concentration near 2e-5, far above the CSV comparison floor. The recovery sink is then
# visible in the result: with K = 0 (constant density) the monomer changes by 1.6%, and a 10% change
# in K by 0.017%, both far above the comparison tolerance. (An earlier version had no source and a
# 3e-5 s run, so its monomer concentration of about 4e-15 was below the comparison floor and the test
# could not fail.)
#
# D0/Em_eV/temperature/Z/rho0/K/burgers_vector below match dislocation_recovery_jacobian.i.

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
  dt = 0.25
  num_steps = 8
[]

[Outputs]
  csv = true
[]
