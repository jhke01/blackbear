# Exercises InterstitialLoopGroupedNodalKernel's binding_energy_model parameter:
# binding_energy_table and formation_energy_table mode, using the real MD binding energies from
# Soneda and Diaz de la Rubia (1998), Table 1 (interstitial clusters n = 2..10, eV):
# 0.71 1.20 2.09 1.67 2.12 1.84 3.00 2.16 1.97. num_cluster_sizes = 30 extends well past the
# table's n=10 end, so the run also exercises the capillary-law fallback for every larger size --
# the one genuinely new behavior this feature adds beyond the Cu-precipitation kernel's own table
# modes (which hard-error past the table).
#
# formation_energy_table_eV below (E_f(1)..E_f(10)) is built from the SAME Table 1 binding
# energies via E_f(n) = E_f(n-1) + E_f(1) - E_b(n), so both table modes reduce to the identical
# binding energies and the cli_args variants below are directly comparable.
#
# Capillary parameters (loop_energy_coefficient_eV/loop_energy_constant_eV) are still required
# (the kernel's capillary fallback is mandatory in every mode) and are set to the paper's own
# fitted law for interstitials, E_f(n) = 5.76*n^(2/3) eV (their Eq 3): A = 5.76, C = 0, p = 2/3.

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
  [interstitial_clusters]
    family = LAGRANGE
    order = FIRST
  []
[]

[Functions]
  [index_fn]
    type = PiecewiseConstant
    x = '0 0.9'
    y = '3.473920870009918e-36 0'
    direction = left
  []
[]

[ICs]
  [interstitial_ic]
    type = GroupedArrayFunctionIndexIC
    variable = interstitial_clusters
    function = index_fn
  []
[]

[NodalKernels]
  [interstitial_dot]
    type = ArrayTimeDerivativeNodalKernel
    variable = interstitial_clusters
  []
  [interstitial_dyn]
    type = InterstitialLoopGroupedNodalKernel
    variable = interstitial_clusters

    generation = 1.5e-4
    sink = 0.0

    temperature = 573
    D0_i = 8.781e-8
    Em_i_eV = 0.4684
    atomic_volume = 1.1664e-29
    lattice_parameter = 0.36e-9

    z_i = 1.10
    z_loop_i = 63.0
    alpha_loop_i = 0.8

    formation_energy_eV = 4.33
    loop_energy_exponent = 0.6666666666666666
    loop_energy_coefficient_eV = 5.76
    loop_energy_constant_eV = 0.0

    binding_energy_table_eV = '0.71 1.20 2.09 1.67 2.12 1.84 3.00 2.16 1.97'
    formation_energy_table_eV = '4.33 7.95 11.08 13.32 15.98 18.19 20.68 22.01 24.18 26.54'
  []
[]

[Postprocessors]
  [interstitial_monomer_concentration]
    type = ClusterSizeConcentrationGrouped
    clusters = interstitial_clusters
    n_size = 1
  []
  [total_loop_density]
    type = ClusterTotalDensityGrouped
    clusters = interstitial_clusters
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
