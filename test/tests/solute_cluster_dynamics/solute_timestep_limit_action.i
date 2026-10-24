# The continuum solute problem of solute_continuum_action.i with the time step limited by the relative change of the resolvable density
# (timestep_relative_change): the action adds SoluteClusterDynamics_timestep_limit and the time stepper uses it.

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1
[]
[Problem]
  use_hash_table_matrix_assembly = true
  restore_original_nonzero_pattern = false
[]
[SoluteClusterDynamics]
  initial_concentration = 0.0134
  rate_model = interfacial_energy
  diffusivity_model = arrhenius
  temperature = 563.15
  D0 = 6.3e-5
  Q_eV = 2.29
  radiation_enhanced_factor = 2.5e7
  sigma = 0.37
  atomic_volume = 1.1782924e-29
  Omega_kB_K = 6255.0
  DeltaS_kB = 0.866
  resolvable_radius = 5e-10
  num_cluster_sizes = 2000
  tail = fokker_planck
  timestep_relative_change = 0.25
  timestep_cap_fraction = 0.1
  timestep_cap_minimum = 1e-2
  switch_size = 20
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
  line_search = none
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-14
  end_time = 0.5
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-3
    optimal_iterations = 6
    iteration_window = 1
    growth_factor = 1.3
    timestep_limiting_postprocessor = SoluteClusterDynamics_timestep_limit
  []
[]
[Outputs]
  csv = true
[]
