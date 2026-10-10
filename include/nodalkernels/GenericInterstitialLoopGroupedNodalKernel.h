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
 * and the sink parameter below, respectively. All grouping/opposite-coupling
 * logic that this kernel shares structurally with
 * GenericVacancyVoidGroupedNodalKernel lives in the shared CRTP base,
 * GenericDefectFamilyGroupedNodalKernel; this class supplies only the
 * loop-specific leaf physics below.
 *
 * The opposite-absorption terms are non-AD only, for the same reason as
 * PointDefectRecombinationNodalKernel: Coupleable::coupledArrayValue()
 * returns a plain, non-AD ArrayVariableValue, and
 * ADArrayNodalKernel::setJacobian() unconditionally errors, so an AD variant
 * could never assemble the resulting cross-variable Jacobian block.
 */
template <bool is_ad>
class GenericInterstitialLoopGroupedNodalKernelTempl
  : public GenericDefectFamilyGroupedNodalKernel<
        GenericInterstitialLoopGroupedNodalKernelTempl<is_ad>, is_ad>
{
public:
  enum class BindingEnergyModel
  {
    CAPILLARY = 0,
    BINDING_ENERGY_TABLE = 1,
    FORMATION_ENERGY_TABLE = 2
  };

  static InputParameters validParams();
  GenericInterstitialLoopGroupedNodalKernelTempl(const InputParameters & parameters);

  Real beta(unsigned int n) const;
  Real alpha(unsigned int n) const;
  Real radius(unsigned int n) const;
  Real formationEnergy(unsigned int n) const;
  /// Binding energy E_b(n) for a loop of size n >= 2: the table-driven value when
  /// binding_energy_model selects a table and n is within its range, otherwise the capillary-law
  /// value (the default mode, and the fallback once n exceeds either table's range).
  Real bindingEnergy(unsigned int n) const;
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

protected:
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

  /// Selected binding-energy model: capillary (default/fallback), or a direct table.
  const BindingEnergyModel _binding_energy_model;
  /// Direct binding energies [eV] indexed as n = index + 2.
  const std::vector<Real> _binding_energy_table_eV;
  /// Formation energies [eV] indexed as n = index + 1.
  const std::vector<Real> _formation_energy_table_eV;

  /// Opposite-species (vacancy) diffusivity prefactor/migration energy, for the loop's own
  /// absorption of a vacancy monomer (Piece A).
  const Real _opposite_D0;
  const Real _opposite_Q;

  /// Vacancy-absorption bias parameters for this loop family's own Piece A term (loop absorbs
  /// vacancy monomer): loopBiasVacancy() uses these in place of z_i/z_loop_i/alpha_loop_i.
  const Real _z_v;
  const Real _z_loop_v;
  const Real _alpha_loop_v;

  usingGenericDefectFamilyGroupedNodalKernelMembers(
      GenericInterstitialLoopGroupedNodalKernelTempl<is_ad>, is_ad);
};

typedef GenericInterstitialLoopGroupedNodalKernelTempl<false> InterstitialLoopGroupedNodalKernel;
typedef GenericInterstitialLoopGroupedNodalKernelTempl<true> ADInterstitialLoopGroupedNodalKernel;
