# Dislocation line length and lattice strain of a loop family with a constant concentration per unit
# size, C = 2 on every discrete size and cell, with the radius law r_n = sqrt(n) (r1 = 1, V_at = 1) and
# n_minimum = 2 (the monomer is not a loop). The exact values are
#   line length = 2 pi C (sum_{n=2}^{20} sqrt(n) + int_{20.5}^{2000.5} sqrt(n) dn) = 749579.3449,
#   lattice strain = C (sum_{n=2}^{20} n + int_{20.5}^{2000.5} n dn) / 3 = 1333999.3333.
# The time-derivative kernel alone keeps the initial state.

[Mesh]
  type = SingleNodeMesh
  position = '1 1 1'
[]
[Variables]
  [u]
    family = LAGRANGE
    order = FIRST
    components = 132
    array = true
  []
[]
[ICs]
  [u_ic]
    type = ArrayConstantIC
    variable = u
    value = '2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2'
  []
[]
[NodalKernels]
  [dot]
    type = ArrayTimeDerivativeNodalKernel
    variable = u
  []
[]
[Postprocessors]
  [line_length]
    type = ClusterContinuumMoment
    clusters = u
    quantity = line_length
    num_cluster_sizes = 2000
    switch_size = 20
    n_minimum = 2
    r1 = 1
    radius_exponent = 0.5
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [lattice_strain]
    type = ClusterContinuumMoment
    clusters = u
    quantity = lattice_strain
    num_cluster_sizes = 2000
    switch_size = 20
    n_minimum = 2
    r1 = 1
    radius_exponent = 0.5
    execute_on = 'INITIAL TIMESTEP_END'
  []
[]
[Executioner]
  type = Transient
  solve_type = NEWTON
  dt = 1
  num_steps = 1
[]
[Outputs]
  csv = true
[]
