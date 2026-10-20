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

#include "GenericArrayNodalKernel.h"
#include "MooseEnum.h"

#include <vector>

/**
 * Combined cluster dynamics array nodal kernel for all cluster sizes (1 through N).
 *
 * Handles the monomer (size-1) and all larger clusters in a single array
 * variable. Array component index i corresponds to cluster size n = i+1:
 *   i=0  ->  n=1  (monomer)
 *   i=1  ->  n=2  (dimer)
 *   ...
 *   i=N-1 -> n=N
 *
 * Unlike the paired ClusterSizeOneNodalKernel + ClusterDynamicsArrayNodalKernel
 * approach, no separate scalar variable is needed. The non-AD version assembles the
 * full intra-variable Jacobian, so Problem/use_hash_table_matrix_assembly is required.
 */
template <bool is_ad>
class GenericClusterDynamicsNodalKernelTempl : public GenericArrayNodalKernel<is_ad>
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
  GenericClusterDynamicsNodalKernelTempl(const InputParameters & parameters);

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

  /// Absorption rate coefficient for the selected rate model
  Real beta(unsigned int n) const;

  /// Emission rate coefficient for the selected rate model
  Real alpha(unsigned int n) const;

  /// d(beta(n))/d(vacancy monomer concentration), nonzero only when vacancy_coupled_var is
  /// given. monomerDiffusivity() is linear in the coupled vacancy monomer value when the
  /// coupling is active, so this derivative is itself independent of n_comp/solution state
  /// beyond that one linear factor -- see monomerDiffusivityCouplingDerivative().
  Real betaPrime(unsigned int n) const;

  /// d(alpha(n))/d(vacancy monomer concentration) = betaPrime(n-1)*exp(-E_b(n)/kT), the same
  /// chain rule alpha(n) itself already uses (alpha(n) = beta(n-1)*exp(-E_b(n)/kT), and E_b(n)
  /// does not depend on diffusivity).
  Real alphaPrime(unsigned int n) const;

  /// d(monomerDiffusivity())/d(vacancy monomer concentration): a true constant when
  /// vacancy_coupled_var is given, since monomerDiffusivity() = radiation_enhanced_factor *
  /// diffusivity_base * (vacancy monomer value / _vacancy_c_thermal), and diffusivity_base
  /// (constant or Arrhenius) does not itself depend on the coupled variable.
  Real monomerDiffusivityCouplingDerivative() const;

  /// Atomic volume used by the interfacial-energy model
  Real atomicVolume() const;

  /// Monomer diffusivity used by the interfacial-energy model
  Real monomerDiffusivity() const;

  /// The constant- or Arrhenius-model diffusivity BEFORE radiation_enhanced_factor or the
  /// vacancy-coupling multiplier are applied. Factored out of monomerDiffusivity() so
  /// monomerDiffusivityCouplingDerivative() can reuse it without re-deriving the same
  /// diffusivity_model branch.
  Real baseDiffusivity() const;

  /// Radius of a cluster of size n used by the interfacial-energy model
  Real radius(unsigned int n) const;

  /// Binding energy for a cluster of size n used by the interfacial-energy model
  Real bindingEnergy(unsigned int n) const;

  /// Build or refresh cached size-dependent coefficients for the current array size
  void ensureCoefficientCache(unsigned int n_comp) const;

  /// Monomer generation rate G_1
  const Real _generation;

  /// Linear sink coefficient k_s
  const Real _sink;

  /// Selected rate model
  const RateModel _rate_model;

  /// Base absorption coefficient beta_0 for the simple model
  const Real _beta0;

  /// Base emission coefficient alpha_0 for the simple model
  const Real _alpha0;

  /// Temperature T for the interfacial-energy model
  const Real _temperature;

  /// Monomer diffusion coefficient for the interfacial-energy model
  const Real _monomer_diffusivity;

  /// Selected diffusivity model for the interfacial-energy rate model
  const DiffusivityModel _diffusivity_model;

  /// Selected binding-energy model for detailed-balance emission
  const BindingEnergyModel _binding_energy_model;

  /// Diffusion prefactor D0 for the Arrhenius diffusivity model
  const Real _D0;

  /// Activation energy Q for the Arrhenius diffusivity model [J]
  const Real _Q;

  /// Radiation-enhanced multiplier applied to the monomer diffusivity
  const Real _radiation_enhanced_factor;

  /// Optional DefectDynamics vacancy grouped array variable: when coupled, the monomer
  /// diffusivity is additionally scaled by the dynamic vacancy supersaturation ratio
  /// C_v(1,t)/_vacancy_c_thermal. Component 0 of this array is always the vacancy monomer,
  /// regardless of that species' own grouping scheme (ClusterGroupingLayout::monomerComponent()
  /// is always 0), so no grouping-layout parameters are needed here. Non-AD only, the same
  /// pattern PointDefectRecombinationNodalKernel already uses: Coupleable::coupledArrayValue()
  /// returns a plain, non-AD ArrayVariableValue, and ADArrayNodalKernel::setJacobian()
  /// unconditionally errors, so an AD variant could never assemble the cross-variable Jacobian
  /// block this coupling needs.
  const ArrayVariableValue & _vacancy_array;

  /// Whether vacancy_coupled_var was actually provided.
  const bool _vacancy_coupling_active;

  /// DOF-map variable number of the coupled vacancy array, for computeQpOffDiagJacobian's own
  /// jvar guard.
  const unsigned int _vacancy_var_num;

  /// Vacancy thermal-equilibrium monomer concentration, 4*exp(-E_f/kT) at this kernel's own
  /// temperature -- computed once at construction from vacancy_formation_energy_eV, using the
  /// SAME kB/eV constants this kernel's own anonymous namespace already defines, so this is
  /// bit-identical to DefectDynamicsAction::monomerConcentration()'s own formula and constants,
  /// not merely numerically close.
  const Real _vacancy_c_thermal;

  /// Interfacial energy sigma for interfacial-energy and thermal-surface binding models
  const Real _sigma;

  /// Atomic volume V_at for the interfacial-energy model
  const Real _atomic_volume;

  /// Enthalpy term Omega for the interfacial-energy model
  const Real _Omega;

  /// Non-configurational entropy term DeltaS for the interfacial-energy model
  const Real _DeltaS;

  /// Generic capillary binding-energy constant [J]
  const Real _binding_energy_constant;

  /// Generic capillary enthalpy term [J]
  const Real _binding_energy_enthalpy;

  /// Generic capillary entropy term [J/K]
  const Real _binding_energy_entropy;

  /// Generic capillary binding-energy coefficient [J]
  const Real _binding_energy_coefficient;

  /// Generic capillary size exponent
  const Real _binding_energy_exponent;

  /// Direct binding energies indexed as n = index + 2 [J]
  const std::vector<Real> _binding_energy_table;

  /// Formation energies indexed as n = index + 1 [J]
  const std::vector<Real> _formation_energy_table;

  /// Cached beta coefficients indexed by cluster size n
  mutable std::vector<Real> _beta_cache;

  /// Cached alpha coefficients indexed by cluster size n
  mutable std::vector<Real> _alpha_cache;

  /// Number of array components represented by the current coefficient cache
  mutable unsigned int _cache_size = 0;

  /// Whether the preconditioner factorization ordering has been checked
  bool _checked_factor_ordering = false;

  usingGenericArrayNodalKernelMembers;
};

typedef GenericClusterDynamicsNodalKernelTempl<false> ClusterDynamicsNodalKernel;
typedef GenericClusterDynamicsNodalKernelTempl<true> ADClusterDynamicsNodalKernel;
