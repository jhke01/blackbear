# CSVDiff test for ClusterDynamicsNodalKernel's vacancy_coupled_var/vacancy_formation_energy_eV
# radiation-enhanced-diffusion coupling. vacancy_clusters is a free array variable (only
# ArrayTimeDerivativeNodalKernel acts on it, so it stays fixed at its IC for the whole run) whose
# component 0 represents the DefectDynamics vacancy monomer concentration C_v(1,t). Its IC is set
# to the thermal-equilibrium value C_v_thermal = 4*exp(-Ef/kT) at Ef = 1.6 eV, T = 600 K (matching
# clusters_dyn's own temperature), so the dynamic supersaturation ratio C_v(1,t)/C_v_thermal = 1
# and the coupling is a pure no-op -- this run must reproduce cluster_dynamics_energy_models.i's
# own physics bit-for-bit (same rate_model/diffusivity_model/temperature/atomic_volume/sigma/
# Omega_kB_K/DeltaS_kB and the same 8-component index IC).

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
    value = '1.4542916048943378e-13'
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
    monomer_diffusivity = 1.0e-21
    temperature = 600
    atomic_volume = 1.1782924e-29
    sigma = 0.37
    Omega_kB_K = 6255
    DeltaS_kB = 0.866
    vacancy_coupled_var = vacancy_clusters
    vacancy_formation_energy_eV = 1.6
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
