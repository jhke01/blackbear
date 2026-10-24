# Vacancy and interstitial loops with a Fokker-Planck continuum tail (tail = fokker_planck, discrete up to size 20) whose interstitial loops join a cellular dislocation network, [Recovery] with model = cellular_network and loop_transfer = true. The loops of size 3 and up leave the distribution at the rate k = v_l/L and the network density gains their line length.
[Mesh]
  type = SingleNodeMesh
  position = '0 0 0'
[]
[Problem]
  use_hash_table_matrix_assembly = true
  restore_original_nonzero_pattern = false
[]
[DefectDynamics]
  temperature = 603.15
  atomic_volume = 1.1664e-29
  lattice_parameter = 3.6e-10
  dpa_rate = 1.0
  cascade_efficiency = 0.3
  recombination_radius = 7e-10
  tail = fokker_planck
  switch_size = 20
  
  # The interstitial loops join the cellular network of Rest (2006): LoopNetworkTransferNodalKernel
  # removes them from the distribution and the network density gains their line length. The bias
  # is a round number, not a measured one. dpa_rate = 1 and the steps are chosen so that the
  # density changes measurably.
  [Recovery]
    model = cellular_network
    loop_transfer = true
    transfer_minimum_size = 3
    rho0 = 1e10
    bias = 0.1
    ca_crho = 2.0
    interstitial_species = interstitial
    burgers_vector = 2.0784609690826528e-10
  []

  [vacancy]
    geometry = loop
    D0 = 6e-05
    Em = 1.35
    Ef = 1.7
    Eb2 = 0.5
    exponent = 0.8
    Z = 1.0
    Z_loop = 33.0
    alpha = 0.65
    Z_loop_opposite = 63.0
    alpha_opposite = 0.8
    num_cluster_sizes = 2000
    resolvable_radius = 7.5e-10
    cascade_f2 = 0.06
    cascade_f3 = 0.03
    cascade_f4 = 0.02
  []

  [interstitial]
    geometry = loop
    D0 = 1e-07
    Em = 0.45
    Ef = 4.1
    Eb2 = 0.6
    exponent = 0.8
    Z = 1.1
    Z_loop = 63.0
    alpha = 0.8
    Z_loop_opposite = 33.0
    alpha_opposite = 0.65
    num_cluster_sizes = 2000
    resolvable_radius = 7.5e-10
    cascade_f2 = 0.5
    cascade_f3 = 0.2
    cascade_f4 = 0.06
  []
[]
[Preconditioning]
  [smp]
    type = SMP
    full = false
    off_diag_row = 'vacancy_clusters_0 interstitial_clusters_0'
    off_diag_column = 'interstitial_clusters_0 vacancy_clusters_0'
  []
[]
[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_ordering_type'
  petsc_options_value = 'lu       nd'
  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-20
  dt = 0.01
  num_steps = 10
[]
[Postprocessors]
  [dt_cap]
    type = FunctionValuePostprocessor
    function = '0.01'
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
[Outputs]
  csv = true
[]
