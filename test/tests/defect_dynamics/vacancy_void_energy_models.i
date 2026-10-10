# Exercises VacancyVoidGroupedNodalKernel's binding_energy_model parameter: binding_energy_table
# and formation_energy_table mode, using the real MD binding energies from Soneda and Diaz de la
# Rubia (1998), Table 1 (vacancy clusters n = 2..8, eV): 0.23 0.37 0.66 0.66 0.96 0.66 0.94.
# num_cluster_sizes = 30 extends well past the table's n=8 end, so the run also exercises the
# capillary-law fallback for every larger size -- the one genuinely new behavior this feature
# adds beyond the Cu-precipitation kernel's own table modes (which hard-error past the table).
#
# formation_energy_table_eV below (E_f(1)..E_f(8)) is built from the SAME Table 1 binding
# energies via E_f(n) = E_f(n-1) + E_f(1) - E_b(n), so both table modes reduce to the identical
# binding energies and the cli_args variants below are directly comparable.
#
# Capillary parameters (void_energy_coefficient_eV/void_energy_constant_eV) are still required
# (the kernel's capillary fallback is mandatory in every mode) and are set to the Soneda fit
# itself (A = 2.59, C = 0, p = 2/3), so the capillary-mode run (binding_energy_model left at its
# default) is also directly comparable to the two table-mode variants below.

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
  use_grouping = false
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

    formation_energy_eV = 1.73
    void_energy_exponent = 0.6666666666666666
    void_energy_coefficient_eV = 2.59
    void_energy_constant_eV = 0.0

    binding_energy_table_eV = '0.23 0.37 0.66 0.66 0.96 0.66 0.94'
    formation_energy_table_eV = '1.73 3.23 4.59 5.66 6.73 7.50 8.57 9.36'
  []
[]

[Postprocessors]
  [vacancy_monomer_concentration]
    type = ClusterSizeConcentrationGrouped
    clusters = vacancy_clusters
    n_size = 1
  []
  [total_void_density]
    type = ClusterTotalDensityGrouped
    clusters = vacancy_clusters
    n_minimum = 2
    atomic_volume = 1.1664e-29
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
