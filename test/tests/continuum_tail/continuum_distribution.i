# The size distribution of a continuum-layout variable with a constant concentration per unit size:
# a row per discrete size (width 1) and per continuum cell (width and center from the layout). The
# cell widths must sum to N - switch_size and the last cell must end at N + 1/2, and every
# concentration is 2.

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
[VectorPostprocessors]
  [distribution]
    type = ClusterContinuumDistribution
    clusters = u
    num_cluster_sizes = 2000
    switch_size = 20
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
