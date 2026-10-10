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

#include "GenericDefectFamilyGroupedNodalKernel.h"
#include "MooseEnum.h"

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
 * respectively. All grouping/opposite-coupling logic that this kernel shares
 * structurally with GenericInterstitialLoopGroupedNodalKernel lives in the
 * shared CRTP base, GenericDefectFamilyGroupedNodalKernel; this class
 * supplies only the void-specific leaf physics below.
 *
 * The opposite-absorption terms are non-AD only, for the same reason as
 * PointDefectRecombinationNodalKernel: Coupleable::coupledArrayValue()
 * returns a plain, non-AD ArrayVariableValue, and
 * ADArrayNodalKernel::setJacobian() unconditionally errors, so an AD variant
 * could never assemble the resulting cross-variable Jacobian block.
 */
template <bool is_ad>
class GenericVacancyVoidGroupedNodalKernelTempl
  : public GenericDefectFamilyGroupedNodalKernel<GenericVacancyVoidGroupedNodalKernelTempl<is_ad>,
                                                 is_ad>
{
public:
  enum class BindingEnergyModel
  {
    CAPILLARY = 0,
    BINDING_ENERGY_TABLE = 1,
    FORMATION_ENERGY_TABLE = 2
  };

  static InputParameters validParams();
  GenericVacancyVoidGroupedNodalKernelTempl(const InputParameters & parameters);

  Real beta(unsigned int n) const;
  Real alpha(unsigned int n) const;
  Real radius(unsigned int n) const;
  Real formationEnergy(unsigned int n) const;
  /// Binding energy E_b(n) for a void of size n >= 2: the table-driven value when
  /// binding_energy_model selects a table and n is within its range, otherwise the capillary-law
  /// value (the default mode, and the fallback once n exceeds either table's range).
  Real bindingEnergy(unsigned int n) const;
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

protected:
  const Real _temperature;
  const Real _D0;
  const Real _Q;
  const Real _atomic_volume;
  const Real _lattice_parameter;
  const Real _formation_energy_eV;
  const Real _energy_exponent;
  const Real _energy_coefficient_eV;
  const Real _energy_constant_eV;

  /// Selected binding-energy model: capillary (default/fallback), or a direct table.
  const BindingEnergyModel _binding_energy_model;
  /// Direct binding energies [eV] indexed as n = index + 2.
  const std::vector<Real> _binding_energy_table_eV;
  /// Formation energies [eV] indexed as n = index + 1.
  const std::vector<Real> _formation_energy_table_eV;

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

  usingGenericDefectFamilyGroupedNodalKernelMembers(GenericVacancyVoidGroupedNodalKernelTempl<is_ad>,
                                                    is_ad);
};

typedef GenericVacancyVoidGroupedNodalKernelTempl<false> VacancyVoidGroupedNodalKernel;
typedef GenericVacancyVoidGroupedNodalKernelTempl<true> ADVacancyVoidGroupedNodalKernel;
