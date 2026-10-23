# Twin of cluster_dynamics_vacancy_coupling_supersaturated.i: no vacancy coupling, monomer_diffusivity 10 times larger.
# See that file.

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1
[]

[Problem]
  use_hash_table_matrix_assembly = true
  restore_original_nonzero_pattern = false
[]

[Variables]
  [clusters]
    family = LAGRANGE
    order = FIRST
    components = 8
  []
  [vacancy_clusters]
    family = LAGRANGE
    order = FIRST
    components = 1
    array = true
  []
[]

[Functions]
  [index_fn]
    type = PiecewiseConstant
    x = '0 0.9'
    y = '0.01 0.001'
    direction = left
  []
[]

[ICs]
  [clusters_ic]
    type = ArrayFunctionIndexIC
    variable = clusters
    function = index_fn
  []
  [vacancy_ic]
    type = ArrayConstantIC
    variable = vacancy_clusters
    # C_v_thermal = 4*exp(-Ef/kT) at Ef = 1.6 eV, T = 600 K -- supersaturation ratio = 1
    value = '1.4542916048943378e-12'
  []
[]

[NodalKernels]
  [clusters_dot]
    type = ArrayTimeDerivativeNodalKernel
    variable = clusters
  []
  [clusters_dyn]
    type = ADClusterDynamicsNodalKernel
    variable = clusters
    generation = 1.0e-3
    sink = 0.1
    rate_model = interfacial_energy
    diffusivity_model = constant
    monomer_diffusivity = 1.0e-20
    temperature = 600
    atomic_volume = 1.1782924e-29
    sigma = 0.37
    Omega_kB_K = 6255
    DeltaS_kB = 0.866
  []
  # No kernel besides the time derivative acts on vacancy_clusters, so it stays pinned at its IC.
  [vacancy_dot]
    type = ArrayTimeDerivativeNodalKernel
    variable = vacancy_clusters
  []
[]

[Postprocessors]
  [avg_cluster_radius]
    type = ClusterAverageRadius
    clusters = clusters
    r1 = 1.0
    n_minimum = 2
  []
  [monomer_concentration]
    type = ClusterSizeConcentration
    clusters = clusters
    n_size = 1
  []
  [total_cluster_density]
    type = ClusterTotalDensity
    clusters = clusters
    n_minimum = 1
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
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_factor_mat_ordering_type'
  petsc_options_value = 'ilu      rcm'
  nl_rel_tol = 1.0e-9
  nl_abs_tol = 1.0e-12
  dt = 0.01
  num_steps = 3
[]

[Outputs]
  csv = true
[]
