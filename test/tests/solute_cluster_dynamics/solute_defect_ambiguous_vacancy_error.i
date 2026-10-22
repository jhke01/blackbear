# Error-path test: a [DefectDynamics] block with two geometry = void species makes the vacancy
# species to couple to ambiguous (DefectDynamicsAction itself separately rejects two void species
# as physically unsupported, but [SoluteClusterDynamics]'s own ambiguity check runs first here
# since it is declared earlier in the input -- MOOSE actions on the same task run in declaration
# order -- and must raise its own clear error rather than silently picking one).

[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1
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
  explicit_max = 10
[]

[DefectDynamics]
  temperature = 600
  atomic_volume = 1.1782924e-29
  lattice_parameter = 0.36e-9
  dpa_rate = 0
  cascade_efficiency = 0
  dislocation_density = 0

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
  [other_void]
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

[Executioner]
  type = Transient
  num_steps = 1
  dt = 1
[]
