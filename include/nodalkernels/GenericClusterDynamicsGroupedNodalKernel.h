/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/*                       BlackBear                              */
/*                                                              */
/*           (c) 2017 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#pragma once

#include "ClusterGroupingLayout.h"
#include "GenericArrayNodalKernel.h"
#include "MooseEnum.h"

#include <vector>

/**
 * Grouped L0/L1 cluster-dynamics nodal kernel.
 *
 * This object preserves the monomer and a user-selected explicit size range,
 * then stores the large-size tail in grouped bins carrying the first two
 * moments of the in-bin size distribution via the linear coefficients L0 and
 * L1. The grouped-tail dynamics follow the same practical moment formulation
 * used in the standalone Cu grouped-moments solver. The non-AD version assembles the
 * full intra-variable Jacobian, so Problem/use_hash_table_matrix_assembly is required.
 */
template <bool is_ad>
class GenericClusterDynamicsGroupedNodalKernelTempl : public GenericArrayNodalKernel<is_ad>
{
public:
  enum class RateModel
  {
    SIMPLE = 0,
    INTERFACIAL_ENERGY = 1
  };

  enum class DiffusivityModel
  {
    CONSTANT = 0,
    ARRHENIUS = 1
  };

  enum class BindingEnergyModel
  {
    INTERFACIAL_ENERGY = 0,
    CAPILLARY = 1,
    BINDING_ENERGY_TABLE = 2,
    FORMATION_ENERGY_TABLE = 3
  };

  static InputParameters validParams();
  GenericClusterDynamicsGroupedNodalKernelTempl(const InputParameters & parameters);

  virtual void jacobianSetup() override;

protected:
  virtual void computeQpResidual(GenericRealEigenVector<is_ad> & residual) override;
  virtual void computeQpJacobian() override;
  // Not declared "override": ArrayNodalKernel (is_ad = false) declares a virtual
  // computeQpOffDiagJacobian(unsigned int) to override, but ADArrayNodalKernel (is_ad = true)
  // has no such member at all -- only computeOffDiagJacobian(), which is AD-computed and never
  // calls a per-qp hook. Specialized below exactly like computeQpJacobian() already is: a real
  // implementation for is_ad = false, mooseError() for is_ad = true.
  void computeQpOffDiagJacobian(unsigned int jvar);

  Real beta(unsigned int n) const;
  Real alpha(unsigned int n) const;

  /// d(beta(n))/d(vacancy monomer concentration), nonzero only when vacancy_coupled_var is
  /// given. See GenericClusterDynamicsNodalKernelTempl's identical member for the derivation.
  Real betaPrime(unsigned int n) const;

  /// d(alpha(n))/d(vacancy monomer concentration) = betaPrime(n-1)*exp(-E_b(n)/kT).
  Real alphaPrime(unsigned int n) const;

  /// d(monomerDiffusivity())/d(vacancy monomer concentration): a true constant when
  /// vacancy_coupled_var is given.
  Real monomerDiffusivityCouplingDerivative() const;

  Real atomicVolume() const;
  Real monomerDiffusivity() const;

  /// The constant- or Arrhenius-model diffusivity BEFORE radiation_enhanced_factor or the
  /// vacancy-coupling multiplier are applied.
  Real baseDiffusivity() const;

  Real radius(unsigned int n) const;
  Real bindingEnergy(unsigned int n) const;
  void ensureCoefficientCache() const;

  template <typename Scalar>
  void fillConcentrationCache(const std::vector<Scalar> & state,
                              std::vector<Scalar> & concentration) const;

  template <typename Scalar>
  void fillFluxCache(const std::vector<Scalar> & concentration, std::vector<Scalar> & flux) const;

  template <typename Scalar>
  void computeGroupedResidualFromState(const std::vector<Scalar> & state,
                                       std::vector<Scalar> & residual_values) const;

  const Real _generation;
  const Real _sink;
  const RateModel _rate_model;
  const Real _beta0;
  const Real _alpha0;
  const Real _temperature;
  const Real _monomer_diffusivity;
  const DiffusivityModel _diffusivity_model;
  const BindingEnergyModel _binding_energy_model;
  const Real _D0;
  const Real _Q;
  const Real _radiation_enhanced_factor;

  /// radiation_enhanced_factor times the base diffusivity: every input is a constant parameter, so
  /// this is evaluated once in the constructor instead of on every beta()/alpha() call.
  const Real _static_diffusivity;

  /// Optional DefectDynamics vacancy grouped array variable: when coupled, the monomer
  /// diffusivity is additionally scaled by the dynamic vacancy supersaturation ratio
  /// C_v(1,t)/_vacancy_c_thermal. Component 0 of this array is always the vacancy monomer,
  /// regardless of that species' own grouping scheme. Non-AD only, see
  /// GenericClusterDynamicsNodalKernelTempl's identical member for the full rationale.
  const ArrayVariableValue & _vacancy_array;

  /// Whether vacancy_coupled_var was actually provided.
  const bool _vacancy_coupling_active;

  /// DOF-map variable number of the coupled vacancy array, for computeQpOffDiagJacobian's own
  /// jvar guard.
  const unsigned int _vacancy_var_num;

  /// Vacancy thermal-equilibrium monomer concentration, 4*exp(-E_f/kT) at this kernel's own
  /// temperature -- computed once at construction from vacancy_formation_energy_eV.
  const Real _vacancy_c_thermal;

  const Real _sigma;
  const Real _atomic_volume;
  const Real _Omega;
  const Real _DeltaS;
  const Real _binding_energy_constant;
  const Real _binding_energy_enthalpy;
  const Real _binding_energy_entropy;
  const Real _binding_energy_coefficient;
  const Real _binding_energy_exponent;
  const std::vector<Real> _binding_energy_table;
  const std::vector<Real> _formation_energy_table;

  const unsigned int _num_cluster_sizes;
  const bool _use_grouping;
  const unsigned int _explicit_max;
  const ClusterGroupingSpacing _group_spacing;
  const unsigned int _group_linear_width;
  const Real _group_geometric_ratio;
  const bool _enforce_group_nonnegative;
  const Real _group_nonnegative_tolerance_factor;

  const ClusterGroupingLayout _layout;

  mutable std::vector<Real> _beta_cache;
  mutable std::vector<Real> _alpha_cache;
  mutable unsigned int _cache_size = 0;

  /// Whether the preconditioner factorization ordering has been checked
  bool _checked_factor_ordering = false;

  usingGenericArrayNodalKernelMembers;
};

typedef GenericClusterDynamicsGroupedNodalKernelTempl<false> ClusterDynamicsGroupedNodalKernel;
typedef GenericClusterDynamicsGroupedNodalKernelTempl<true> ADClusterDynamicsGroupedNodalKernel;
