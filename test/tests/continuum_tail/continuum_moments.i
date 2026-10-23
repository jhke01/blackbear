# Moments of a continuum-layout distribution with a constant concentration per unit size, for which
# the exact values are known: with C = 2 on every discrete size and every cell, the number of clusters
# is C N = 4000, the content C N (N + 1) / 2 = 4002000 and, for the radius law r_n = n, the mean
# radius (N + 1) / 2 = 1000.5. With n_minimum = 11 the number is C (N - 10) = 3980, the content
# is C (N (N + 1) / 2 - 55) = 3.99989e6 and the mean radius is (N + 11) / 2 = 1005.5. The time-derivative
# kernel alone keeps the initial state.

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
  [density]
    type = ClusterContinuumMoment
    clusters = u
    quantity = density
    num_cluster_sizes = 2000
    switch_size = 20
    n_minimum = 1
    r1 = 1
    radius_exponent = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [content]
    type = ClusterContinuumMoment
    clusters = u
    quantity = content
    num_cluster_sizes = 2000
    switch_size = 20
    n_minimum = 1
    r1 = 1
    radius_exponent = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [radius]
    type = ClusterContinuumMoment
    clusters = u
    quantity = average_radius
    num_cluster_sizes = 2000
    switch_size = 20
    n_minimum = 1
    r1 = 1
    radius_exponent = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [density_min]
    type = ClusterContinuumMoment
    clusters = u
    quantity = density
    num_cluster_sizes = 2000
    switch_size = 20
    n_minimum = 11
    r1 = 1
    radius_exponent = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [content_min]
    type = ClusterContinuumMoment
    clusters = u
    quantity = content
    num_cluster_sizes = 2000
    switch_size = 20
    n_minimum = 11
    r1 = 1
    radius_exponent = 1
    execute_on = 'INITIAL TIMESTEP_END'
  []
  [radius_min]
    type = ClusterContinuumMoment
    clusters = u
    quantity = average_radius
    num_cluster_sizes = 2000
    switch_size = 20
    n_minimum = 11
    r1 = 1
    radius_exponent = 1
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
