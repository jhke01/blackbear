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

#include <vector>

/**
 * Grouped L0/L1 vacancy-void nodal kernel.
 *
 * Vacancy clusters (voids) are three-dimensional: the void radius grows as
 * n^(1/3), the same functional form as the existing Cu-precipitation grouped
 * kernel's spherical-cluster radius. The self-absorption rate coefficient
 * includes a finite-sink-strength correction factor
 * 1/(1 + lattice_parameter/radius), and emission uses a capillary-law
 * formation energy A_v*n^p_v + C_v rather than an interfacial-energy model.
 *
 * In addition to self-reactions, this kernel adds the opposite-species
 * absorption terms from DefectDynamics_Sundials.cpp's
 * addGroupedFamilyRHS/oppositeShrinkFluxAt: an interstitial monomer absorbed
 * onto a void shrinks it by one vacancy (using this family's own radius and
 * finite-sink-strength coefficient, with the opposite species' diffusivity
 * D_i), and the corresponding loss from the interstitial monomer pool, which
 * requires summing the coupled interstitial-loop array using THAT family's
 * own radius/bias-factor coefficients. Point-defect recombination and
 * dislocation sink loss are added by a separate coupling kernel
 * (PointDefectRecombinationNodalKernel) and the sink parameter below,
 * respectively.
 *
 * The opposite-absorption terms are non-AD only, for the same reason as
 * PointDefectRecombinationNodalKernel: Coupleable::coupledArrayValue()
 * returns a plain, non-AD ArrayVariableValue, and
 * ADArrayNodalKernel::setJacobian() unconditionally errors, so an AD variant
 * could never assemble the resulting cross-variable Jacobian block.
 */
template <bool is_ad>
class GenericVacancyVoidGroupedNodalKernelTempl : public GenericArrayNodalKernel<is_ad>
{
public:
  static InputParameters validParams();
  GenericVacancyVoidGroupedNodalKernelTempl(const InputParameters & parameters);

  virtual void jacobianSetup() override;

protected:
  virtual void computeQpResidual(GenericRealEigenVector<is_ad> & residual) override;
  virtual void computeQpJacobian() override;
  // Not declared "override": ArrayNodalKernel (is_ad = false) declares a virtual
  // computeQpOffDiagJacobian(unsigned int) to override, but ADArrayNodalKernel (is_ad = true) has
  // no such member at all -- only the public, non-Qp computeOffDiagJacobian(). The non-AD
  // specialization below still dispatches correctly through ArrayNodalKernel's own virtual call;
  // the AD specialization is a plain no-op member, never invoked by the framework, matching the
  // non-AD-only limitation documented in the class comment above.
  virtual void computeQpOffDiagJacobian(unsigned int jvar);

  Real beta(unsigned int n) const;
  Real alpha(unsigned int n) const;
  Real radius(unsigned int n) const;
  Real formationEnergy(unsigned int n) const;
  Real monomerDiffusivity() const;
  void ensureCoefficientCache() const;

  /// Rate coefficient at which a void of size n absorbs an interstitial monomer, shrinking to
  /// size n-1 (Piece A): the same finite-sink-strength form as beta(n), but with the opposite
  /// species' diffusivity D_i in place of D_v. Matches vacancy.beta_opposite[n] in the reference.
  Real betaOpposite(unsigned int n) const;
  Real oppositeMonomerDiffusivity() const;

  /// Opposite (interstitial loop) radius, matching InterstitialLoopGroupedNodalKernel::radius().
  Real oppositeRadius(unsigned int n) const;
  /// Opposite (interstitial loop) dislocation bias factor for vacancy absorption, matching
  /// InterstitialLoopGroupedNodalKernel::loopBias() but with the vacancy-side parameters
  /// (DD_AMSS_Z_V/Z_LOOP_V/ALPHA_LOOP_V), i.e. loopBiasVacancy() in the reference code.
  Real oppositeBias(unsigned int n) const;
  /// beta_opposite for the interstitial-loop family: the rate coefficient at which a loop of size
  /// n absorbs a vacancy monomer, shrinking to size n-1. Needed to sum the interstitial monomer
  /// loss (Piece B) over the coupled loop array.
  Real oppositeBetaOpposite(unsigned int n) const;
  void ensureOppositeCoefficientCache() const;

  template <typename Scalar>
  void fillConcentrationCache(const std::vector<Scalar> & state,
                              std::vector<Scalar> & concentration) const;

  template <typename Scalar>
  void fillFluxCache(const std::vector<Scalar> & concentration,
                     const std::vector<Scalar> & opposite_concentration,
                     Scalar opposite_monomer,
                     std::vector<Scalar> & flux) const;

  template <typename Scalar>
  void computeGroupedResidualFromState(const std::vector<Scalar> & state,
                                       const std::vector<Scalar> & opposite_state,
                                       std::vector<Scalar> & residual_values) const;

  const Real _generation;
  const Real _sink;
  const Real _temperature;
  const Real _D0;
  const Real _Q;
  const Real _atomic_volume;
  const Real _lattice_parameter;
  const Real _formation_energy_eV;
  const Real _energy_exponent;
  const Real _energy_coefficient_eV;
  const Real _energy_constant_eV;

  const unsigned int _num_cluster_sizes;
  const bool _use_grouping;
  const unsigned int _explicit_max;
  const ClusterGroupingSpacing _group_spacing;
  const unsigned int _group_linear_width;
  const Real _group_geometric_ratio;
  const bool _enforce_group_nonnegative;
  const Real _group_nonnegative_tolerance_factor;

  const ClusterGroupingLayout _layout;

  /// Opposite-species (interstitial) diffusivity prefactor/migration energy, for the void's own
  /// absorption of an interstitial monomer (Piece A).
  const Real _opposite_D0;
  const Real _opposite_Q;

  /// Interstitial-loop geometry/bias parameters needed only to evaluate the loop's own
  /// beta_opposite coefficient when summing the coupled loop array (Piece B). Not used for this
  /// kernel's own radius/absorption physics.
  const Real _opposite_burgers_vector;
  const Real _opposite_z_v;
  const Real _opposite_z_loop_v;
  const Real _opposite_alpha_loop_v;

  const unsigned int _opposite_num_cluster_sizes;
  const bool _opposite_use_grouping;
  const unsigned int _opposite_explicit_max;
  const ClusterGroupingSpacing _opposite_group_spacing;
  const unsigned int _opposite_group_linear_width;
  const Real _opposite_group_geometric_ratio;
  const bool _opposite_enforce_group_nonnegative;
  const Real _opposite_group_nonnegative_tolerance_factor;

  const ClusterGroupingLayout _opposite_layout;

  const ArrayVariableValue & _opposite_array;
  const unsigned int _opposite_var_num;
  /// Whether opposite_coupled_var was actually provided. When false, D0_i/opposite_* default to
  /// values (D0_i = 0, opposite_num_cluster_sizes = 1) that make every opposite-absorption term
  /// evaluate to exactly zero, so the kernel reduces to its pre-existing self-reaction-only
  /// behavior without any extra branching in the residual/Jacobian -- this flag only gates the
  /// constructor's physical-validity checks on the opposite-side parameters.
  const bool _opposite_coupling_active;

  mutable std::vector<Real> _beta_cache;
  mutable std::vector<Real> _alpha_cache;
  mutable std::vector<Real> _beta_opposite_cache;
  mutable unsigned int _cache_size = 0;

  mutable std::vector<Real> _opposite_beta_opposite_cache;
  mutable unsigned int _opposite_cache_size = 0;

  /// Whether the preconditioner factorization ordering has been checked
  bool _checked_factor_ordering = false;

  usingGenericArrayNodalKernelMembers;
};

typedef GenericVacancyVoidGroupedNodalKernelTempl<false> VacancyVoidGroupedNodalKernel;
typedef GenericVacancyVoidGroupedNodalKernelTempl<true> ADVacancyVoidGroupedNodalKernel;
