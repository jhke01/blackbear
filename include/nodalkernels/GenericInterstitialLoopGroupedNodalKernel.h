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
 * geometry. Only the monomer absorption/emission terms for this family are
 * computed here; point-defect recombination and dislocation sink loss are
 * added by a separate coupling kernel (PointDefectRecombinationNodalKernel)
 * and the sink parameter below, respectively.
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

  Real beta(unsigned int n) const;
  Real alpha(unsigned int n) const;
  Real radius(unsigned int n) const;
  Real formationEnergy(unsigned int n) const;
  Real loopBias(unsigned int n) const;
  Real monomerDiffusivity() const;
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

  mutable std::vector<Real> _beta_cache;
  mutable std::vector<Real> _alpha_cache;
  mutable unsigned int _cache_size = 0;

  /// Whether the preconditioner factorization ordering has been checked
  bool _checked_factor_ordering = false;

  usingGenericArrayNodalKernelMembers;
};

typedef GenericInterstitialLoopGroupedNodalKernelTempl<false> InterstitialLoopGroupedNodalKernel;
typedef GenericInterstitialLoopGroupedNodalKernelTempl<true> ADInterstitialLoopGroupedNodalKernel;
