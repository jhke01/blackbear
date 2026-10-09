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
 * Grouped L0/L1 vacancy-void self-reaction nodal kernel.
 *
 * Vacancy clusters (voids) are three-dimensional: the void radius grows as
 * n^(1/3), the same functional form as the existing Cu-precipitation grouped
 * kernel's spherical-cluster radius. Unlike that kernel, the absorption rate
 * coefficient here includes a finite-sink-strength correction factor
 * 1/(1 + lattice_parameter/radius), and emission uses a capillary-law
 * formation energy A_v*n^p_v + C_v rather than an interfacial-energy model.
 * Only the monomer absorption/emission terms for this family are computed
 * here; point-defect recombination and dislocation sink loss are added by a
 * separate coupling kernel (PointDefectRecombinationNodalKernel) and the
 * sink parameter below, respectively.
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

  Real beta(unsigned int n) const;
  Real alpha(unsigned int n) const;
  Real radius(unsigned int n) const;
  Real formationEnergy(unsigned int n) const;
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

typedef GenericVacancyVoidGroupedNodalKernelTempl<false> VacancyVoidGroupedNodalKernel;
typedef GenericVacancyVoidGroupedNodalKernelTempl<true> ADVacancyVoidGroupedNodalKernel;
