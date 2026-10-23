# ContinuumLoopNodalKernel in time: a loop family that is discrete up to size 20 and a Fokker-Planck
# continuum above it (149 components, largest size 4000), minmod limiter on, starting from an empty
# distribution with a constant monomer source and a linear sink. The time step is capped in proportion
# to the time (0.01 t + 0.5 s): IterationAdaptiveDT alone lets the step grow faster than the transient
# can be followed.

[Mesh]
  type = SingleNodeMesh
  position = '0 0 0'
[]
[Problem]
  use_hash_table_matrix_assembly = true
  restore_original_nonzero_pattern = false
[]
[Variables]
  [u]
    family = LAGRANGE
    order = FIRST
    components = 149
    array = true
  []
[]
[ICs]
  [u_ic]
    type = ArrayConstantIC
    variable = u
    value = '0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0'
  []
[]
[NodalKernels]
  [dot]
    type = ArrayTimeDerivativeNodalKernel
    variable = u
  []
  [dyn]
    type = ContinuumLoopNodalKernel
    variable = u
    generation = 8.699999999999998e-08
    sink = 3.1468192651509726e-06
    temperature = 603.15
    atomic_volume = 1.1664e-29
    lattice_parameter = 0.36e-9
    burgers_vector = 2.078460969082653e-10
    D0 = 6e-5
    Em_eV = 1.35
    Z = 1.0
    Z_loop = 33.0
    alpha = 0.65
    formation_energy_eV = 1.7
    energy_exponent = 0.8
    energy_coefficient_eV = 1.619212219414472
    energy_constant_eV = 0.08078778058552794
    num_cluster_sizes = 4000
    switch_size = 20
    use_limiter = true
  []
[]
[AuxVariables]
  [monomer_aux]
    family = LAGRANGE
    order = FIRST
  []
[]
[AuxKernels]
  [monomer_aux]
    type = ArrayVarReductionAux
    variable = monomer_aux
    array_variable = u
    value_type = sum
    weights = '1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
[Postprocessors]
  [dt_cap]
    type = ParsedPostprocessor
    expression = '0.01 * t + 0.5'
    use_t = true
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [monomer]
    type = NodalVariableValue
    nodeid = 0
    variable = monomer_aux
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
  line_search = none
  scheme = bdf2
  start_time = 0
  end_time = 2000.0
  l_tol = 1.0e-4
  l_max_its = 30
  nl_max_its = 20
  nl_rel_tol = 1.0e-8
  nl_abs_tol = 1.0e-18
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.5
    optimal_iterations = 6
    iteration_window = 1
    growth_factor = 1.3
    timestep_limiting_postprocessor = dt_cap
  []
[]
[Outputs]
  # Only the monomer, at the times where it was compared with an independent CVODE integration of the
  # same rate function (monomer 2.5884e-6, 1.3737e-6, 8.8294e-7, 6.1383e-7, 4.6987e-7, 3.5775e-7).
  [out]
    type = CSV
    show = monomer
    sync_times = '50 100 200 500 1000 2000'
    sync_only = true
  []
[]
