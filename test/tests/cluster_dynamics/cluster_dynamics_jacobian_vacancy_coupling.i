# Jacobian test for the new off-diagonal block in ClusterDynamicsNodalKernel's
# vacancy_coupled_var coupling (GenericClusterDynamicsNodalKernelTempl<false>::
# computeQpOffDiagJacobian). vacancy_clusters is a second array variable (only
# ArrayTimeDerivativeNodalKernel acts on it, so this is purely a Jacobian structure test, not a
# physical solve) whose component 0 is read as the DefectDynamics vacancy monomer concentration.
# Its IC is set away from the thermal-equilibrium value (S_v = C_v(1)/C_v_thermal != 1) so the
# coupling derivative terms (betaPrime/alphaPrime) are genuinely nonzero and exercised by
# -snes_test_jacobian, verifying the hand-coded cross-variable Jacobian against PETSc's
# finite-difference Jacobian.

[Mesh]
  type = SingleNodeMesh
  position = '1 1 1'
[]

[Problem]
  use_hash_table_matrix_assembly = true
  restore_original_nonzero_pattern = false
[]

[Variables]
  [clusters]
    family = LAGRANGE
    order = FIRST
    components = 11
  []
  [vacancy_clusters]
    family = LAGRANGE
    order = FIRST
    components = 1
    array = true
  []
[]

[ICs]
  [clusters_ic]
    type = ArrayConstantIC
    variable = clusters
    # index 0 = monomer (size 1), indices 1-10 = clusters sizes 2-11
    value = '1.5 0.05 0.04 0.03 0.025 0.02 0.015 0.01 0.008 0.006 0.004'
  []
  [vacancy_ic]
    type = ArrayConstantIC
    variable = vacancy_clusters
    # 3x the thermal-equilibrium concentration at Ef = 1.6 eV, T = 600 K, so S_v != 1
    value = '4.362874814683014e-13'
  []
[]

[NodalKernels]
  [clusters_dot]
    type = ArrayTimeDerivativeNodalKernel
    variable = clusters
  []
  [clusters_dyn]
    type = ClusterDynamicsNodalKernel
    variable = clusters
    generation = 0.05
    sink = 0.15
    rate_model = interfacial_energy
    diffusivity_model = constant
    monomer_diffusivity = 1.0e-21
    temperature = 600
    atomic_volume = 1.1782924e-29
    binding_energy_model = formation_energy_table
    formation_energy_table_eV = '0 0.05 0.12 0.2 0.3 0.41 0.53 0.66 0.8 0.95 1.11'
    vacancy_coupled_var = vacancy_clusters
    vacancy_formation_energy_eV = 1.6
  []
  # No kernel besides the time derivative acts on vacancy_clusters.
  [vacancy_dot]
    type = ArrayTimeDerivativeNodalKernel
    variable = vacancy_clusters
  []
[]

[Preconditioning]
  [smp]
    # SMP's off_diag_row/off_diag_column and coupled_groups do not resolve array-variable logical
    # names (only their per-component names, e.g. vacancy_clusters_0) -- a MOOSE-side limitation
    # unrelated to this kernel. Affordable at this test's small component count.
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_ordering_type'
  petsc_options_value = 'lu       nd'
  dt = 0.01
  num_steps = 1
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-8
[]

[Outputs]
  console = true
[]
