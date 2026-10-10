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

#include "ClusterDynamicsSolverCheck.h"
#include "ClusterGroupingLayout.h"
#include "FEProblemBase.h"
#include "GenericArrayNodalKernel.h"
#include "NonlinearSystemBase.h"

#include "libmesh/libmesh_common.h"

#include <vector>

/**
 * Shared L0/L1 grouped-tail logic for a two-species defect-dynamics family (vacancy voids or
 * interstitial loops), common to GenericVacancyVoidGroupedNodalKernel and
 * GenericInterstitialLoopGroupedNodalKernel.
 *
 * This base owns the grouping layout, the opposite-family coupling parameters, and every piece
 * of residual/Jacobian logic that is structurally identical between the two families: filling
 * concentration/flux caches from a grouped array state, assembling the grouped residual
 * (self-absorption/emission plus the opposite-species shrink terms, see "Piece A"/"Piece B" in
 * [theory/defect_dynamics.md]), and the exact intra- and cross-variable Jacobians obtained by
 * chaining through those fluxes. It dispatches to the derived class (CRTP, via derived()) for
 * exactly the leaf physics that differs between the two families: the cluster radius law (n^1/3
 * for voids, sqrt(n) for loops), the self- and opposite-absorption rate coefficients, and the
 * corresponding opposite-family coefficient used to sum the Piece B drain. Member templates are
 * only instantiated where actually called, so Derived need only be complete at that point (the
 * usual CRTP header-only pattern); the body therefore lives entirely in this header.
 *
 * The opposite-absorption terms are non-AD only, for the same reason as
 * PointDefectRecombinationNodalKernel: Coupleable::coupledArrayValue() returns a plain, non-AD
 * ArrayVariableValue, and ADArrayNodalKernel::setJacobian() unconditionally errors, so an AD
 * variant could never assemble the resulting cross-variable Jacobian block.
 */
template <typename Derived, bool is_ad>
class GenericDefectFamilyGroupedNodalKernel : public GenericArrayNodalKernel<is_ad>
{
public:
  static InputParameters validParams();
  GenericDefectFamilyGroupedNodalKernel(const InputParameters & parameters);

  virtual void jacobianSetup() override;

protected:
  virtual void computeQpResidual(GenericRealEigenVector<is_ad> & residual) override;
  virtual void computeQpJacobian() override;
  // Not declared "override": ArrayNodalKernel (is_ad = false) declares a virtual
  // computeQpOffDiagJacobian(unsigned int) to override, but ADArrayNodalKernel (is_ad = true) has
  // no such member at all -- only the public, non-Qp computeOffDiagJacobian(). The non-AD
  // specialization still dispatches correctly through ArrayNodalKernel's own virtual call; the AD
  // specialization is a plain no-op member, never invoked by the framework, matching the
  // non-AD-only limitation documented in the class comment above.
  virtual void computeQpOffDiagJacobian(unsigned int jvar);

  Derived & derived() { return static_cast<Derived &>(*this); }
  const Derived & derived() const { return static_cast<const Derived &>(*this); }

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

  const unsigned int _num_cluster_sizes;
  const bool _use_grouping;
  const unsigned int _explicit_max;
  const ClusterGroupingSpacing _group_spacing;
  const unsigned int _group_linear_width;
  const Real _group_geometric_ratio;
  const bool _enforce_group_nonnegative;
  const Real _group_nonnegative_tolerance_factor;

  const ClusterGroupingLayout _layout;

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
  /// Whether opposite_coupled_var was actually provided. When false, the opposite-side rate
  /// coefficients default to values that make every opposite-absorption term evaluate to exactly
  /// zero, so the kernel reduces to its self-reaction-only behavior without any extra branching
  /// in the residual/Jacobian -- this flag only gates the constructor's physical-validity checks
  /// on the opposite-side parameters, performed by the derived class.
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

#define usingGenericDefectFamilyGroupedNodalKernelMembers(Derived, is_ad)                        \
  using GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::_generation;                      \
  using GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::_sink;                             \
  using GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::_num_cluster_sizes;                \
  using GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::_use_grouping;                     \
  using GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::_explicit_max;                     \
  using GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::_group_spacing;                    \
  using GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::_group_linear_width;               \
  using GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::_group_geometric_ratio;            \
  using GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::_enforce_group_nonnegative;        \
  using GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::_group_nonnegative_tolerance_factor; \
  using GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::_layout;                           \
  using GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::_opposite_num_cluster_sizes;       \
  using GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::_opposite_use_grouping;            \
  using GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::_opposite_explicit_max;            \
  using GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::_opposite_group_spacing;           \
  using GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::_opposite_group_linear_width;      \
  using GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::_opposite_group_geometric_ratio;   \
  using GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::_opposite_enforce_group_nonnegative; \
  using GenericDefectFamilyGroupedNodalKernel<Derived,                                           \
                                              is_ad>::_opposite_group_nonnegative_tolerance_factor; \
  using GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::_opposite_layout;                  \
  using GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::_opposite_array;                   \
  using GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::_opposite_var_num;                 \
  using GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::_opposite_coupling_active;         \
  using GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::_beta_cache;                       \
  using GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::_alpha_cache;                      \
  using GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::_beta_opposite_cache;              \
  using GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::_cache_size;                       \
  using GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::_opposite_beta_opposite_cache;     \
  using GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::_opposite_cache_size;              \
  using GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::_checked_factor_ordering;          \
  usingGenericArrayNodalKernelMembers

template <typename Derived, bool is_ad>
InputParameters
GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::validParams()
{
  InputParameters params = GenericArrayNodalKernel<is_ad>::validParams();
  MooseEnum group_spacing("linear geometric", "linear");
  params.addParam<Real>(
      "sink", 0.0, "Linear dislocation sink coefficient for this family's monomer loss.");
  params.addRequiredParam<unsigned int>(
      "num_cluster_sizes", "Maximum physical cluster size represented by the grouped distribution.");
  params.addParam<bool>(
      "use_grouping", true, "Whether to activate grouped-tail storage above explicit_max.");
  params.addParam<unsigned int>(
      "explicit_max", 100, "Maximum cluster size retained explicitly before grouped bins.");
  params.addParam<MooseEnum>(
      "group_spacing", group_spacing, "Grouped-tail spacing mode: linear or geometric.");
  params.addParam<unsigned int>(
      "group_linear_width", 25, "Fixed grouped-bin width when group_spacing = linear.");
  params.addParam<Real>("group_geometric_ratio",
                        1.35,
                        "Geometric grouped-bin growth factor when group_spacing = geometric.");
  params.addParam<bool>(
      "enforce_group_nonnegative",
      true,
      "Whether to limit grouped L0/L1 reconstruction to keep bin endpoint concentrations "
      "nonnegative (up to the tolerance factor).");
  params.addParam<Real>(
      "group_nonnegative_tolerance_factor",
      0.0,
      "Relaxation factor for grouped nonnegativity enforcement. Endpoint concentrations are "
      "allowed down to -factor * L0 before limiting.");

  params.addCoupledVar(
      "opposite_coupled_var",
      0.0,
      "The opposite-family grouped array variable. Optional: when not provided, the "
      "opposite-species absorption terms (Piece A and Piece B) are inactive and this kernel "
      "reduces exactly to its self-reaction-only behavior.");
  params.addParam<unsigned int>(
      "opposite_num_cluster_sizes",
      1,
      "Maximum physical cluster size represented by the coupled grouped distribution. Only "
      "meaningful when opposite_coupled_var is provided.");
  params.addParam<bool>(
      "opposite_use_grouping", true, "Whether the coupled array uses grouped-tail storage.");
  params.addParam<unsigned int>(
      "opposite_explicit_max", 100, "Maximum cluster size stored explicitly in the coupled array.");
  params.addParam<MooseEnum>(
      "opposite_group_spacing", group_spacing, "Grouped-tail spacing mode for the coupled array.");
  params.addParam<unsigned int>("opposite_group_linear_width",
                                25,
                                "Fixed grouped-bin width for the coupled array when "
                                "opposite_group_spacing = linear.");
  params.addParam<Real>("opposite_group_geometric_ratio",
                        1.35,
                        "Geometric grouped-bin growth factor for the coupled array when "
                        "opposite_group_spacing = geometric.");
  params.addParam<bool>("opposite_enforce_group_nonnegative",
                        true,
                        "Whether the coupled array's grouped reconstruction enforces "
                        "nonnegative bin endpoints.");
  params.addParam<Real>("opposite_group_nonnegative_tolerance_factor",
                        0.0,
                        "Relaxation factor for the coupled array's grouped nonnegativity "
                        "enforcement.");
  return params;
}

template <typename Derived, bool is_ad>
GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::GenericDefectFamilyGroupedNodalKernel(
    const InputParameters & parameters)
  : GenericArrayNodalKernel<is_ad>(parameters),
    _generation(this->template getParam<Real>("generation")),
    _sink(this->template getParam<Real>("sink")),
    _num_cluster_sizes(this->template getParam<unsigned int>("num_cluster_sizes")),
    _use_grouping(this->template getParam<bool>("use_grouping")),
    _explicit_max(this->template getParam<unsigned int>("explicit_max")),
    _group_spacing(this->template getParam<MooseEnum>("group_spacing")
                       .template getEnum<ClusterGroupingSpacing>()),
    _group_linear_width(this->template getParam<unsigned int>("group_linear_width")),
    _group_geometric_ratio(this->template getParam<Real>("group_geometric_ratio")),
    _enforce_group_nonnegative(this->template getParam<bool>("enforce_group_nonnegative")),
    _group_nonnegative_tolerance_factor(
        this->template getParam<Real>("group_nonnegative_tolerance_factor")),
    _layout(_num_cluster_sizes,
            _explicit_max,
            _use_grouping,
            _group_spacing,
            _group_linear_width,
            _group_geometric_ratio),
    _opposite_num_cluster_sizes(this->template getParam<unsigned int>("opposite_num_cluster_sizes")),
    _opposite_use_grouping(this->template getParam<bool>("opposite_use_grouping")),
    _opposite_explicit_max(this->template getParam<unsigned int>("opposite_explicit_max")),
    _opposite_group_spacing(this->template getParam<MooseEnum>("opposite_group_spacing")
                                .template getEnum<ClusterGroupingSpacing>()),
    _opposite_group_linear_width(
        this->template getParam<unsigned int>("opposite_group_linear_width")),
    _opposite_group_geometric_ratio(
        this->template getParam<Real>("opposite_group_geometric_ratio")),
    _opposite_enforce_group_nonnegative(
        this->template getParam<bool>("opposite_enforce_group_nonnegative")),
    _opposite_group_nonnegative_tolerance_factor(
        this->template getParam<Real>("opposite_group_nonnegative_tolerance_factor")),
    _opposite_layout(_opposite_num_cluster_sizes,
                     _opposite_explicit_max,
                     _opposite_use_grouping,
                     _opposite_group_spacing,
                     _opposite_group_linear_width,
                     _opposite_group_geometric_ratio),
    _opposite_array(this->coupledArrayValue("opposite_coupled_var")),
    _opposite_var_num(this->coupled("opposite_coupled_var")),
    _opposite_coupling_active(this->isCoupled("opposite_coupled_var"))
{
  // Grouped and explicit components are coupled through the monomer and the inter-bin fluxes,
  // so the Jacobian contains entries outside the default diagonal-block sparsity pattern.
  if (!this->_fe_problem.useHashTableMatrixAssembly())
    mooseError(this->type(), " requires Problem/use_hash_table_matrix_assembly = true");
  if (this->_tid == 0)
    ClusterDynamics::checkSolverSetup(*this,
                                      this->_fe_problem,
                                      this->_var.sys().number(),
                                      this->_var.number(),
                                      this->_var.count());

  if (_group_nonnegative_tolerance_factor < 0.0)
    mooseError(this->type(), " requires group_nonnegative_tolerance_factor >= 0.");
  // The opposite-absorption parameters only need to be physically valid when the opposite
  // coupling is actually active: when opposite_coupled_var is not provided, this kernel reduces
  // exactly to its pre-existing self-reaction-only behavior (see validParams), and the unused
  // defaults are never evaluated by any residual/Jacobian path -- all of the Piece A/Piece B
  // terms are gated on _opposite_coupling_active. The derived class checks its own opposite-side
  // rate-coefficient parameters (diffusivity, bias factors, etc.) after this base constructor
  // runs; this base only checks the one parameter it owns directly.
  if (_opposite_coupling_active && _opposite_group_nonnegative_tolerance_factor < 0.0)
    mooseError(this->type(),
              " requires opposite_group_nonnegative_tolerance_factor >= 0 when "
              "opposite_coupled_var is provided.");
}

template <typename Derived, bool is_ad>
template <typename Scalar>
void
GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::fillConcentrationCache(
    const std::vector<Scalar> & state, std::vector<Scalar> & concentration) const
{
  concentration.assign(_num_cluster_sizes + 1, Scalar(0.0));
  for (unsigned int n = 1; n <= _num_cluster_sizes; ++n)
    concentration[n] = _layout.concentration(
        state, n, _enforce_group_nonnegative, _group_nonnegative_tolerance_factor);
}

template <typename Derived, bool is_ad>
template <typename Scalar>
void
GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::fillFluxCache(
    const std::vector<Scalar> & concentration,
    const std::vector<Scalar> & opposite_concentration,
    Scalar opposite_monomer,
    std::vector<Scalar> & flux) const
{
  // combined_flux[n] = self_flux[n] - K[n+1], where K[n] = betaOpposite(n)*opposite_monomer*C_n is
  // the rate at which an opposite-species monomer absorbed onto a cluster of size n shrinks it to
  // n-1 (oppositeShrinkFluxAt in the reference). Substituting combined_flux for the plain
  // self-flux reproduces the reference's explicit-size and grouped-bin residual formulas
  // unchanged; only the monomer equation needs an extra explicit term beyond this substitution.
  flux.assign(_num_cluster_sizes + 1, Scalar(0.0));
  const Scalar c1 = concentration[1];
  for (unsigned int n = 1; n < _num_cluster_sizes; ++n)
  {
    const Scalar self_flux = derived().beta(n) * c1 * concentration[n] -
                             derived().alpha(n + 1) * concentration[n + 1];
    const Scalar shrink_flux = n + 1 >= 2 ? derived().betaOpposite(n + 1) * opposite_monomer *
                                               concentration[n + 1]
                                         : Scalar(0.0);
    flux[n] = self_flux - shrink_flux;
  }
}

template <typename Derived, bool is_ad>
template <typename Scalar>
void
GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::computeGroupedResidualFromState(
    const std::vector<Scalar> & state,
    const std::vector<Scalar> & opposite_state,
    std::vector<Scalar> & residual_values) const
{
  std::vector<Scalar> concentration;
  fillConcentrationCache(state, concentration);

  std::vector<Scalar> opposite_concentration(_opposite_num_cluster_sizes + 1, Scalar(0.0));
  for (unsigned int n = 1; n <= _opposite_num_cluster_sizes; ++n)
    opposite_concentration[n] = _opposite_layout.concentration(
        opposite_state, n, _opposite_enforce_group_nonnegative,
        _opposite_group_nonnegative_tolerance_factor);
  const Scalar opposite_monomer = opposite_concentration[1];

  std::vector<Scalar> flux;
  fillFluxCache(concentration, opposite_concentration, opposite_monomer, flux);

  residual_values.assign(_layout.componentCount(), Scalar(0.0));

  // Monomer: unlike the explicit/grouped sizes, the monomer equation keeps the PLAIN self-flux
  // (not combined_flux) and gets one extra +K[2] term (a dimer-or-larger cluster shrinking to a
  // monomer adds back exactly one unit; oppositeShrinkFluxAt is not doubled the way the self-flux
  // dimerization term is). It also gains the Piece B drain: the rate at which the ENTIRE coupled
  // opposite-family array is absorbing this kernel's own monomer.
  Scalar self_monomer_flux = 2.0 * (derived().beta(1) * concentration[1] * concentration[1] -
                                    derived().alpha(2) * concentration[2]);
  for (unsigned int n = 2; n < _num_cluster_sizes; ++n)
    self_monomer_flux += derived().beta(n) * concentration[1] * concentration[n] -
                         derived().alpha(n + 1) * concentration[n + 1];
  const Scalar shrink_to_monomer = _num_cluster_sizes >= 2
                                       ? derived().betaOpposite(2) * opposite_monomer *
                                             concentration[2]
                                       : Scalar(0.0);

  Scalar opposite_monomer_loss = Scalar(0.0);
  for (unsigned int n = 2; n <= _opposite_num_cluster_sizes; ++n)
    opposite_monomer_loss +=
        derived().oppositeBetaOpposite(n) * concentration[1] * opposite_concentration[n];

  residual_values[_layout.monomerComponent()] =
      -(_generation - _sink * concentration[1] - self_monomer_flux + shrink_to_monomer) +
      opposite_monomer_loss;

  for (unsigned int n = 2; n <= _layout.explicitMax(); ++n)
    residual_values[_layout.explicitComponent(n)] =
        -((n < _num_cluster_sizes ? flux[n - 1] - flux[n] : flux[n - 1]));

  for (const auto & bin : _layout.groups())
  {
    const auto j_left = flux[bin.start - 1];
    const auto j_right = bin.end < _num_cluster_sizes ? flux[bin.end] : Scalar(0.0);

    residual_values[bin.l0_component] = -((j_left - j_right) / bin.width);

    if (bin.sigma2 > 0.0 && bin.width > 1.0)
    {
      Scalar j_avg = 0.0;
      for (unsigned int n = bin.start; n < bin.end; ++n)
        j_avg += flux[n];
      j_avg /= (bin.width - 1.0);

      residual_values[bin.l1_component] =
          (j_left - 2.0 * j_avg + j_right) * (bin.width - 1.0) / (2.0 * bin.sigma2 * bin.width);
    }
    else
      residual_values[bin.l1_component] = 0.0;
  }
}

template <typename Derived, bool is_ad>
void
GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::jacobianSetup()
{
  GenericArrayNodalKernel<is_ad>::jacobianSetup();

  // The preconditioner is configured from the PETSc options only once the solve starts, so the
  // factorization ordering is checked before the first Jacobian is factored
  if (_checked_factor_ordering || this->_tid != 0)
    return;
  if (auto * nl = dynamic_cast<NonlinearSystemBase *>(&this->_sys))
    _checked_factor_ordering = ClusterDynamics::checkFactorOrdering(*this, *nl);
}

template <typename Derived, bool is_ad>
void
GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::computeQpResidual(
    GenericRealEigenVector<is_ad> & residual)
{
  const auto n_comp = _u[_qp].size();
  if (n_comp != _layout.componentCount())
    mooseError(this->type(),
              " expected ",
              _layout.componentCount(),
              " array components from the grouping layout, but the variable has ",
              n_comp,
              ".");

  derived().ensureCoefficientCache();
  derived().ensureOppositeCoefficientCache();
  std::vector<GenericReal<is_ad>> state(n_comp, GenericReal<is_ad>(0.0));
  for (unsigned int i = 0; i < n_comp; ++i)
    state[i] = _u[_qp](i);

  const auto opposite_n_comp = _opposite_array[_qp].size();
  if (opposite_n_comp != _opposite_layout.componentCount())
    mooseError(this->type(),
              " expected ",
              _opposite_layout.componentCount(),
              " array components from the opposite grouping layout, but opposite_coupled_var "
              "has ",
              opposite_n_comp,
              ".");
  std::vector<GenericReal<is_ad>> opposite_state(opposite_n_comp, GenericReal<is_ad>(0.0));
  for (unsigned int i = 0; i < opposite_n_comp; ++i)
    opposite_state[i] = _opposite_array[_qp](i);

  std::vector<GenericReal<is_ad>> residual_values;
  computeGroupedResidualFromState(state, opposite_state, residual_values);
  residual.resize(n_comp);
  for (unsigned int i = 0; i < n_comp; ++i)
    residual(i) = residual_values[i];
}

// Non-AD Jacobian: every residual row is a linear combination of the fluxes
// J_n = beta(n)*C_1*C_n - alpha(n+1)*C_{n+1}, and every reconstructed concentration is linear in
// the state, C_n = a0(n)*L0 + a1(n)*L1 for grouped sizes (with the limiter coefficients frozen).
// The chain rule through the fluxes gives the exact intra-variable Jacobian in O(N) work.
template <typename Derived, bool is_ad>
void
GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::computeQpJacobian()
{
  if constexpr (is_ad)
  {
    mooseError("Internal error: computeQpJacobian should never be called for the AD version");
  }
  else
  {
    const auto n_comp = static_cast<unsigned int>(_u[_qp].size());
    if (n_comp != _layout.componentCount())
      mooseError(this->type(),
                " expected ",
                _layout.componentCount(),
                " array components from the grouping layout, but the variable has ",
                n_comp,
                ".");

    derived().ensureCoefficientCache();

    const RealEigenVector & state = _u[_qp];
    std::vector<Real> concentration(_num_cluster_sizes + 1, 0.0);
    std::vector<Real> a0(_num_cluster_sizes + 1, 1.0);
    std::vector<Real> a1(_num_cluster_sizes + 1, 0.0);
    for (unsigned int n = 1; n <= _num_cluster_sizes; ++n)
    {
      if (_layout.isExplicitSize(n))
      {
        concentration[n] = state(_layout.explicitComponent(n));
        continue;
      }
      const auto & bin = _layout.groupForSize(n);
      const Real l0 = state(bin.l0_component);
      const Real l1 = state(bin.l1_component);
      ClusterGroupingLayout::reconstructionCoefficients(bin,
                                                        n,
                                                        l0,
                                                        l1,
                                                        _enforce_group_nonnegative,
                                                        _group_nonnegative_tolerance_factor,
                                                        a0[n],
                                                        a1[n]);
      concentration[n] = a0[n] * l0 + a1[n] * l1;
    }
    const Real c1 = concentration[1];

    derived().ensureOppositeCoefficientCache();
    const Real opposite_monomer = _opposite_array[_qp](0);
    Real opposite_coefficient_sum = 0.0;
    for (unsigned int n = 2; n <= _opposite_num_cluster_sizes; ++n)
    {
      const Real opposite_c_n = _opposite_layout.concentration(
          _opposite_array[_qp], n, _opposite_enforce_group_nonnegative,
          _opposite_group_nonnegative_tolerance_factor);
      opposite_coefficient_sum += derived().oppositeBetaOpposite(n) * opposite_c_n;
    }

    // Accumulate one Jacobian row at a time, visiting only the columns it touches
    std::vector<Real> row_values(n_comp, 0.0);
    std::vector<bool> row_touched(n_comp, false);
    std::vector<unsigned int> row_columns;
    auto addEntry = [&](unsigned int column, Real value)
    {
      if (!row_touched[column])
      {
        row_touched[column] = true;
        row_columns.push_back(column);
      }
      row_values[column] += value;
    };
    auto addConcentrationDerivative = [&](unsigned int n, Real value)
    {
      if (_layout.isExplicitSize(n))
        addEntry(_layout.explicitComponent(n), value);
      else
      {
        const auto & bin = _layout.groupForSize(n);
        addEntry(bin.l0_component, value * a0[n]);
        addEntry(bin.l1_component, value * a1[n]);
      }
    };
    // Plain self-flux derivative (alpha only): used ONLY for the monomer row's self_monomer_flux
    // term, which by construction (see computeGroupedResidualFromState) does not include the
    // opposite-absorption K terms -- those are added to the monomer row as separate explicit
    // terms below, not folded into this helper.
    auto addFluxDerivative = [&](unsigned int n, Real weight)
    {
      // J_0 and J_N vanish at the ends of the truncated size space
      if (n == 0 || n >= _num_cluster_sizes)
        return;
      addConcentrationDerivative(1, weight * derived().beta(n) * concentration[n]);
      addConcentrationDerivative(n, weight * derived().beta(n) * c1);
      addConcentrationDerivative(n + 1, -weight * derived().alpha(n + 1));
    };
    // Combined-flux derivative (alpha + opposite-shrink): used for explicit-size and grouped-bin
    // rows, where the residual's flux[n] means combined_flux[n] = self_flux[n] -
    // betaOpposite(n+1)*opposite_monomer*C_{n+1}. The opposite-shrink term only ever modifies the
    // existing alpha(n+1) coefficient's column (C_{n+1}), so it is folded in as an "effective
    // alpha" -- opposite_monomer is a fixed real value within this Jacobian evaluation, so this
    // remains a purely intra-variable (own-array) derivative.
    auto addCombinedFluxDerivative = [&](unsigned int n, Real weight)
    {
      if (n == 0 || n >= _num_cluster_sizes)
        return;
      addConcentrationDerivative(1, weight * derived().beta(n) * concentration[n]);
      addConcentrationDerivative(n, weight * derived().beta(n) * c1);
      const Real effective_alpha =
          derived().alpha(n + 1) + derived().betaOpposite(n + 1) * opposite_monomer;
      addConcentrationDerivative(n + 1, -weight * effective_alpha);
    };
    auto setRow = [&](unsigned int row)
    {
      for (const auto column : row_columns)
      {
        setJacobian(row, column, row_values[column]);
        row_values[column] = 0.0;
        row_touched[column] = false;
      }
      row_columns.clear();
    };

    // Monomer: F = -(G_1 - k_s*C_1 - self_monomer_flux + shrink_to_monomer) + opposite_monomer_loss
    addEntry(_layout.monomerComponent(), _sink);
    for (unsigned int n = 1; n < _num_cluster_sizes; ++n)
      addFluxDerivative(n, n == 1 ? 2.0 : 1.0);
    // -shrink_to_monomer = -betaOpposite(2)*opposite_monomer*C_2: intra-variable derivative wrt C_2
    if (_num_cluster_sizes >= 2)
      addConcentrationDerivative(2, -derived().betaOpposite(2) * opposite_monomer);
    // +opposite_monomer_loss = C_1 * sum_n oppositeBetaOpposite(n)*C_opposite,n: intra-variable
    // derivative wrt C_1 is the (fixed, within this call) coefficient sum itself.
    addConcentrationDerivative(1, opposite_coefficient_sum);
    setRow(_layout.monomerComponent());

    // Explicit sizes: F = -(J_{n-1} - J_n), J now meaning combined_flux
    for (unsigned int n = 2; n <= _layout.explicitMax(); ++n)
    {
      addCombinedFluxDerivative(n - 1, -1.0);
      addCombinedFluxDerivative(n, 1.0);
      setRow(_layout.explicitComponent(n));
    }

    for (const auto & bin : _layout.groups())
    {
      // L0: F = -(J_left - J_right) / width
      addCombinedFluxDerivative(bin.start - 1, -1.0 / bin.width);
      addCombinedFluxDerivative(bin.end, 1.0 / bin.width);
      setRow(bin.l0_component);

      // L1: F = (J_left - 2*J_avg + J_right) * (width - 1) / (2*sigma2*width)
      if (bin.sigma2 > 0.0 && bin.width > 1.0)
      {
        const Real factor = (bin.width - 1.0) / (2.0 * bin.sigma2 * bin.width);
        addCombinedFluxDerivative(bin.start - 1, factor);
        for (unsigned int n = bin.start; n < bin.end; ++n)
          addCombinedFluxDerivative(n, -2.0 * factor / (bin.width - 1.0));
        addCombinedFluxDerivative(bin.end, factor);
        setRow(bin.l1_component);
      }
    }
  }
}

// Cross-variable (opposite-array) Jacobian for the opposite-absorption terms. Non-AD only, same
// reasoning as PointDefectRecombinationNodalKernel: the AD variant's computeOffDiagJacobian only
// re-runs the AD residual graph, which never touched the opposite array's own AD derivatives
// (coupledArrayValue strips them), so it would silently assemble a zero block here.
template <typename Derived, bool is_ad>
void
GenericDefectFamilyGroupedNodalKernel<Derived, is_ad>::computeQpOffDiagJacobian(
    unsigned int jvar)
{
  if constexpr (is_ad)
  {
    // Intentionally a no-op: see the class comment and PointDefectRecombinationNodalKernel for
    // why the AD variant cannot assemble this cross-variable block.
    return;
  }
  else
  {
    if (jvar != _opposite_var_num)
      return;

    derived().ensureCoefficientCache();
    derived().ensureOppositeCoefficientCache();

    const RealEigenVector & state = _u[_qp];
    std::vector<Real> concentration(_num_cluster_sizes + 1, 0.0);
    for (unsigned int n = 1; n <= _num_cluster_sizes; ++n)
      concentration[n] = _layout.concentration(
          state, n, _enforce_group_nonnegative, _group_nonnegative_tolerance_factor);

    const auto opposite_n_comp = static_cast<unsigned int>(_opposite_array[_qp].size());
    std::vector<Real> opposite_a0(_opposite_num_cluster_sizes + 1, 1.0);
    std::vector<Real> opposite_a1(_opposite_num_cluster_sizes + 1, 0.0);
    std::vector<Real> opposite_concentration(_opposite_num_cluster_sizes + 1, 0.0);
    for (unsigned int n = 1; n <= _opposite_num_cluster_sizes; ++n)
    {
      if (_opposite_layout.isExplicitSize(n))
      {
        opposite_concentration[n] = _opposite_array[_qp](_opposite_layout.explicitComponent(n));
        continue;
      }
      const auto & bin = _opposite_layout.groupForSize(n);
      const Real l0 = _opposite_array[_qp](bin.l0_component);
      const Real l1 = _opposite_array[_qp](bin.l1_component);
      ClusterGroupingLayout::reconstructionCoefficients(bin,
                                                        n,
                                                        l0,
                                                        l1,
                                                        _opposite_enforce_group_nonnegative,
                                                        _opposite_group_nonnegative_tolerance_factor,
                                                        opposite_a0[n],
                                                        opposite_a1[n]);
      opposite_concentration[n] = opposite_a0[n] * l0 + opposite_a1[n] * l1;
    }

    std::vector<Real> row_values(opposite_n_comp, 0.0);
    std::vector<bool> row_touched(opposite_n_comp, false);
    std::vector<unsigned int> row_columns;
    auto addEntry = [&](unsigned int column, Real value)
    {
      if (!row_touched[column])
      {
        row_touched[column] = true;
        row_columns.push_back(column);
      }
      row_values[column] += value;
    };
    // d(combined_flux[n])/d(opposite monomer) = -betaOpposite(n+1)*C_{n+1}(own): a single entry
    // at the opposite array's monomer component (0).
    auto addFluxOffDiagDerivative = [&](unsigned int n, Real weight)
    {
      if (n == 0 || n >= _num_cluster_sizes)
        return;
      addEntry(0, -weight * derived().betaOpposite(n + 1) * concentration[n + 1]);
    };
    auto addOppositeConcentrationDerivative = [&](unsigned int n, Real value)
    {
      if (_opposite_layout.isExplicitSize(n))
        addEntry(_opposite_layout.explicitComponent(n), value);
      else
      {
        const auto & bin = _opposite_layout.groupForSize(n);
        addEntry(bin.l0_component, value * opposite_a0[n]);
        addEntry(bin.l1_component, value * opposite_a1[n]);
      }
    };
    auto setRow = [&](unsigned int row)
    {
      for (const auto column : row_columns)
      {
        setJacobian(row, column, row_values[column]);
        row_values[column] = 0.0;
        row_touched[column] = false;
      }
      row_columns.clear();
    };

    // Monomer row: d(-shrink_to_monomer)/d(opposite monomer) = -betaOpposite(2)*C_2(own); plus
    // d(opposite_monomer_loss)/d(opposite C_n) = C_1(own) * oppositeBetaOpposite(n) for every
    // opposite size n >= 2 (dense row against the opposite array, mirroring the existing
    // intra-variable Jacobian's own dense monomer row against its own array).
    if (_num_cluster_sizes >= 2)
      addEntry(0, -derived().betaOpposite(2) * concentration[2]);
    for (unsigned int n = 2; n <= _opposite_num_cluster_sizes; ++n)
      addOppositeConcentrationDerivative(n, concentration[1] * derived().oppositeBetaOpposite(n));
    setRow(_layout.monomerComponent());

    // Explicit sizes and grouped bins: only the single opposite-monomer entry from
    // addFluxOffDiagDerivative; no other dependence on the opposite array.
    for (unsigned int n = 2; n <= _layout.explicitMax(); ++n)
    {
      addFluxOffDiagDerivative(n - 1, -1.0);
      addFluxOffDiagDerivative(n, 1.0);
      setRow(_layout.explicitComponent(n));
    }
    for (const auto & bin : _layout.groups())
    {
      addFluxOffDiagDerivative(bin.start - 1, -1.0 / bin.width);
      addFluxOffDiagDerivative(bin.end, 1.0 / bin.width);
      setRow(bin.l0_component);

      if (bin.sigma2 > 0.0 && bin.width > 1.0)
      {
        const Real factor = (bin.width - 1.0) / (2.0 * bin.sigma2 * bin.width);
        addFluxOffDiagDerivative(bin.start - 1, factor);
        for (unsigned int n = bin.start; n < bin.end; ++n)
          addFluxOffDiagDerivative(n, -2.0 * factor / (bin.width - 1.0));
        addFluxOffDiagDerivative(bin.end, factor);
        setRow(bin.l1_component);
      }
    }
  }
}
