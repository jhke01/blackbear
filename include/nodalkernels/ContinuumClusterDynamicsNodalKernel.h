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

#include "ArrayNodalKernel.h"
#include "ClusterContinuumFlux.h"
#include "ClusterContinuumLayout.h"

/**
 * Solute clustering (for example Cu precipitation) with a size distribution that is discrete up to
 * switch_size and a Fokker-Planck continuum above it, the alternative to the L0/L1 grouping of
 * ClusterDynamicsGroupedNodalKernel. The rate coefficients are those of that kernel for
 * rate_model = interfacial_energy: beta(n) = 4 pi (r(1) + r(n)) D / V_at for spherical clusters
 * r(n) = (3 n V_at / (4 pi))^(1/3), and detailed-balance emission alpha(n) = beta(n - 1)
 * exp(-E_b(n) / kT) from the interfacial-energy (or capillary) binding energy. The array variable
 * has one component for the monomer, one for each discrete size 2 .. switch_size, and one for each
 * continuum cell (ClusterContinuumLayout).
 *
 * The monomer diffusivity can be scaled by the vacancy supersaturation C_v(1, t) / C_v_thermal of a
 * coupled array variable (vacancy_coupled_var), with the off-diagonal Jacobian column of that
 * variable's monomer.
 */
class ContinuumClusterDynamicsNodalKernel : public ArrayNodalKernel
{
public:
  static InputParameters validParams();
  ContinuumClusterDynamicsNodalKernel(const InputParameters & parameters);

protected:
  virtual void computeQpResidual(RealEigenVector & residual) override;
  virtual void computeQpJacobian() override;
  virtual void computeQpOffDiagJacobian(unsigned int jvar) override;

  /// Monomer diffusivity, with the radiation-enhanced factor and the vacancy supersaturation
  Real monomerDiffusivity() const;
  /// Monomer diffusivity before the vacancy supersaturation
  Real staticDiffusivity() const { return _static_diffusivity; }

  /// Radius of a spherical cluster of n atoms
  Real radius(unsigned int n) const;
  /// Binding energy [J] of a cluster of n atoms, n >= 2
  Real bindingEnergy(unsigned int n) const;
  /// beta(n) / D and alpha(n) / D, which do not depend on the (possibly changing) diffusivity
  Real betaShape(unsigned int n) const;
  Real alphaShape(unsigned int n) const;

  /// Build the coefficient shapes on first use
  void ensureShapes();
  /// Scale the shapes by the current diffusivity into the coefficient tables
  void applyDiffusivity();

  const Real _generation;
  const Real _sink;
  const Real _temperature;
  const Real _atomic_volume;
  const Real _static_diffusivity;

  const bool _simple;
  const Real _beta0;
  const Real _alpha0;

  const Real _sigma;
  const Real _omega;
  const Real _delta_s;
  const Real _binding_constant;
  const Real _binding_enthalpy;
  const Real _binding_entropy;
  const Real _binding_coefficient;
  const Real _binding_exponent;
  const bool _capillary;

  const ArrayVariableValue & _vacancy_array;
  const bool _vacancy_active;
  const unsigned int _vacancy_var_num;
  const Real _vacancy_c_thermal;

  const ClusterContinuumLayout _layout;
  const ClusterContinuumFlux _flux;

  /// Coefficient tables at unit diffusivity, and scaled by the current diffusivity
  ClusterContinuumCoefficients _shapes, _coefficients;
  bool _have_shapes = false;
  Real _applied_diffusivity = -1.0;
  std::vector<Real> _rate;
  std::vector<ClusterContinuumFlux::Triplet> _triplets;
};
