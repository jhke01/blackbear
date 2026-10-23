# The continuum version of solute_defect_coupled_action.i: the solute array is a Fokker-Planck continuum (tail = fokker_planck) and the
# [DefectDynamics] vacancy species is found and wired to its kernel exactly as for the grouped array.
# [SoluteClusterDynamics] + [DefectDynamics] together: the solute kernel's vacancy_coupled_var and
# vacancy_formation_energy_eV are wired automatically to the [DefectDynamics] block's vacancy species
# (identified by geometry = void), with no coupling parameter typed by hand.
#
# The vacancies are produced by irradiation here (dpa_rate, a dislocation sink), so the vacancy monomer
# concentration rises far above its thermal value and scales the solute monomer diffusivity by the
# supersaturation C_v / C_v_thermal. That makes the result depend on the coupling: the same solute
# problem without the [DefectDynamics] block gives a solute monomer concentration that is 0.1% higher
# after 4 steps, against a comparison tolerance of 5.5e-6, so the test fails if the auto-wiring is lost.
#
# The preconditioner names the pair of the solute array and the vacancy array. MOOSE calls the
# kernel's off-diagonal Jacobian only for pairs in the coupling matrix, and SMP full = false leaves
# this pair out, so without these two lines the vacancy block is never assembled and Newton needs
# smaller steps.

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1
[]

[Problem]
  use_hash_table_matrix_assembly = true
  restore_original_nonzero_pattern = false
[]

[DefectDynamics]
  temperature = 600
  atomic_volume = 1.1782924e-29
  lattice_parameter = 0.36e-9
  dpa_rate = 1e-6
  cascade_efficiency = 1
  dislocation_density = 1e14

  [vacancy]
    geometry = void
    D0 = 1.19073105525e-8
    Em = 0.7644
    Ef = 1.6
    Eb2 = 0.5
    exponent = 0.67
    num_cluster_sizes = 30
    explicit_max = 6
    group_spacing = linear
    group_linear_width = 4
    group_geometric_ratio = 1.5
  []
[]

[SoluteClusterDynamics]
  initial_concentration = 0.01
  rate_model = interfacial_energy
  diffusivity_model = constant
  monomer_diffusivity = 1.0e-21
  temperature = 600
  atomic_volume = 1.1782924e-29
  sigma = 0.37
  Omega_kB_K = 6255
  DeltaS_kB = 0.866

  num_cluster_sizes = 2000
  tail = fokker_planck
  switch_size = 20
[]

[Preconditioning]
  [smp]
    type = SMP
    full = false
    off_diag_row = 'SoluteClusterDynamics_clusters_0'
    off_diag_column = 'vacancy_clusters_0'
  []
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_factor_mat_ordering_type'
  petsc_options_value = 'ilu      rcm'
  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-12
  dt = 5e-5
  num_steps = 4
[]

[Outputs]
  csv = true
[]
