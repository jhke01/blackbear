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
 * Grouped L0/L1 interstitial-loop self-reaction nodal kernel.
 *
 * Interstitial loops are two-dimensional: the loop radius grows as
 * sqrt(n), not n^(1/3). Unlike the vacancy-void kernel, the absorption rate
 * coefficient here uses a size-dependent dislocation-loop bias factor
 * Z_i(n) in place of a finite-sink-strength correction, and emission uses a
 * capillary-law formation energy A_i*n^p_i + C_i, mirroring the void
 * kernel's capillary emission model but with the loop's own exponent and
 * geometry.
 *
 * In addition to self-reactions, this kernel adds the opposite-species
 * absorption terms from DefectDynamics_Sundials.cpp's
 * addGroupedFamilyRHS/oppositeShrinkFluxAt, mirroring
 * VacancyVoidGroupedNodalKernel: a vacancy monomer absorbed onto a loop
 * shrinks it by one interstitial (using this family's own radius and a
 * vacancy-absorption bias factor loopBiasVacancy, with the opposite
 * species' diffusivity D_v), and the corresponding loss from the vacancy
 * monomer pool, which requires summing the coupled void array using THAT
 * family's own radius/finite-sink-strength coefficient (with this
 * kernel's own D_i). Point-defect recombination and dislocation sink loss
 * are added by a separate coupling kernel (PointDefectRecombinationNodalKernel)
 * and the sink parameter below, respectively.
 *
 * The opposite-absorption terms are non-AD only, for the same reason as
 * PointDefectRecombinationNodalKernel: Coupleable::coupledArrayValue()
 * returns a plain, non-AD ArrayVariableValue, and
 * ADArrayNodalKernel::setJacobian() unconditionally errors, so an AD variant
 * could never assemble the resulting cross-variable Jacobian block.
 */
template <bool is_ad>
class GenericInterstitialLoopGroupedNodalKernelTempl : public GenericArrayNodalKernel<is_ad>
{
public:
  static InputParameters validParams();
  GenericInterstitialLoopGroupedNodalKernelTempl(const InputParameters & parameters);

  virtual void jacobianSetup() override;

protected:
  virtual void computeQpResidual(GenericRealEigenVector<is_ad> & residual) override;
  virtual void computeQpJacobian() override;
  // Not declared "override": see GenericVacancyVoidGroupedNodalKernel.h for why
  // ADArrayNodalKernel has no such base virtual to override, only ArrayNodalKernel does.
  virtual void computeQpOffDiagJacobian(unsigned int jvar);

  Real beta(unsigned int n) const;
  Real alpha(unsigned int n) const;
  Real radius(unsigned int n) const;
  Real formationEnergy(unsigned int n) const;
  Real loopBias(unsigned int n) const;
  Real monomerDiffusivity() const;
  void ensureCoefficientCache() const;

  /// Rate coefficient at which a loop of size n absorbs a vacancy monomer, shrinking to size n-1
  /// (Piece A): same functional form as beta(n) (2*pi*r_i(n)*bias*D/V_at), but with the
  /// vacancy-absorption bias factor loopBiasVacancy() and the opposite species' diffusivity D_v in
  /// place of D_i. Matches interstitial.beta_opposite[n] in the reference.
  Real betaOpposite(unsigned int n) const;
  Real oppositeMonomerDiffusivity() const;
  /// loopBiasVacancy(): same functional form as loopBias() (self-absorption bias), but with the
  /// vacancy-side bias parameters (DD_AMSS_Z_V/Z_LOOP_V/ALPHA_LOOP_V in the reference).
  Real biasVacancy(unsigned int n) const;

  /// Opposite (vacancy void) radius, matching VacancyVoidGroupedNodalKernel::radius().
  Real oppositeRadius(unsigned int n) const;
  /// beta_opposite for the vacancy-void family: the rate coefficient at which a void of size n
  /// absorbs an interstitial monomer, shrinking to size n-1. Uses this kernel's own D_i (no new
  /// diffusivity parameter needed: the void's finite-sink-strength coefficient needs D_i, which
  /// this kernel already computes for its own self-absorption). Needed to sum the vacancy monomer
  /// loss (Piece B) over the coupled void array.
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
  const Real _burgers_vector;
  const Real _z_i;
  const Real _z_loop_i;
  const Real _alpha_loop_i;
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

  /// Opposite-species (vacancy) diffusivity prefactor/migration energy, for the loop's own
  /// absorption of a vacancy monomer (Piece A).
  const Real _opposite_D0;
  const Real _opposite_Q;

  /// Vacancy-absorption bias parameters for this loop family's own Piece A term (loop absorbs
  /// vacancy monomer): loopBiasVacancy() uses these in place of z_i/z_loop_i/alpha_loop_i.
  const Real _z_v;
  const Real _z_loop_v;
  const Real _alpha_loop_v;

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
  /// Whether opposite_coupled_var was actually provided. See
  /// GenericVacancyVoidGroupedNodalKernel.h for why this only gates the constructor's
  /// physical-validity checks on the opposite-side parameters.
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

typedef GenericInterstitialLoopGroupedNodalKernelTempl<false> InterstitialLoopGroupedNodalKernel;
typedef GenericInterstitialLoopGroupedNodalKernelTempl<true> ADInterstitialLoopGroupedNodalKernel;
