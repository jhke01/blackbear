# [SoluteClusterDynamics] + [DefectDynamics] together: confirms the solute kernel's
# vacancy_coupled_var/vacancy_formation_energy_eV are auto-wired to the [DefectDynamics] block's
# vacancy species (identified by geometry = void), with no coupling parameter typed by hand.
# vacancy's own dynamics are deliberately held fixed (generation = sink = 0, no cascade/
# recombination/opposite species), so its monomer concentration stays at the thermal-equilibrium
# initial condition [DefectDynamics] always derives (4*exp(-Ef/kT) at Ef = 1.6 eV, T = 600 K) --
# this test's own goal is confirming the auto-detected WIRING itself is correct (verified directly
# against solute_defect_coupled_equivalent.i's hand-typed vacancy_coupled_var pointing at the same
# variable/formation energy, which produces bit-identical postprocessor values); the coupling's
# own non-trivial physical effect at S_v != 1 is already covered by the kernel-level
# cluster_dynamics_jacobian_vacancy_coupling/cluster_dynamics_vacancy_coupling tests.

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
  dpa_rate = 0
  cascade_efficiency = 0
  dislocation_density = 0

  [Grid]
    explicit_max = 6
    spacing = linear
    linear_width = 4
    geometric_ratio = 1.5
  []

  [vacancy]
    geometry = void
    D0 = 1.19073105525e-8
    Em = 0.7644
    Ef = 1.6
    Eb2 = 0.5
    exponent = 0.67
    num_cluster_sizes = 30
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

  num_cluster_sizes = 50
  use_grouping = true
  explicit_max = 10
  group_spacing = linear
  group_linear_width = 5
  group_geometric_ratio = 1.5
[]

[Preconditioning]
  [smp]
    type = SMP
    full = false
  []
[]

[Executioner]
  type = Transient
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_factor_mat_ordering_type'
  petsc_options_value = 'ilu      rcm'
  nl_rel_tol = 1e-9
  nl_abs_tol = 1e-12
  dt = 0.01
  num_steps = 3
[]

[Outputs]
  csv = true
[]
