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

#include "PointDefectRecombinationNodalKernel.h"

registerMooseObject("BlackBearApp", PointDefectRecombinationNodalKernel);

InputParameters
PointDefectRecombinationNodalKernel::validParams()
{
  InputParameters params = ArrayNodalKernel::validParams();
  params.addClassDescription(
      "Point-defect recombination loss on this variable's monomer component (component 0), "
      "coupled to the other defect family's own array variable. Apply one instance per family "
      "(each coupled to the other) to represent recombination symmetrically.");
  params.addRequiredParam<Real>(
      "recombination_rate",
      "Recombination rate coefficient [m^3/s], 4*pi*r_rec*(D_i+D_v)/V_at.");
  params.addRequiredCoupledVar(
      "coupled_var", "The other defect family's grouped array variable.");
  params.addParam<bool>(
      "assemble_opposite_jacobian_in_own_block",
      false,
      "EXPERIMENTAL. If true, the cross-variable Jacobian entry against coupled_var is written "
      "from this kernel's own-variable Jacobian call, using the other variable's dof indices, "
      "so that the preconditioner needs no off_diag_row/off_diag_column pair.");
  return params;
}

PointDefectRecombinationNodalKernel::PointDefectRecombinationNodalKernel(
    const InputParameters & parameters)
  : ArrayNodalKernel(parameters),
    _recombination_rate(getParam<Real>("recombination_rate")),
    _coupled_monomer(coupledArrayValue("coupled_var")),
    _coupled_var_num(coupled("coupled_var")),
    _assemble_in_own_block(getParam<bool>("assemble_opposite_jacobian_in_own_block"))
{
  if (_recombination_rate < 0.0)
    mooseError("PointDefectRecombinationNodalKernel requires recombination_rate >= 0.");
  // The off-diagonal Jacobian is only called for a variable other than this one, so coupling the
  // kernel to its own variable would silently drop the second term of d(C_i*C_v)/dC.
  if (_coupled_var_num == _var.number())
    paramError("coupled_var", "must be the other defect family's variable, not this one.");
}

void
PointDefectRecombinationNodalKernel::computeQpResidual(RealEigenVector & residual)
{
  residual.resize(_count);
  residual.setZero();
  residual(0) = _recombination_rate * _u[_qp](0) * _coupled_monomer[_qp](0);
}

void
PointDefectRecombinationNodalKernel::computeQpJacobian()
{
  setJacobian(0, 0, _recombination_rate * _coupled_monomer[_qp](0));
  if (_assemble_in_own_block)
    addJacobianElement(_assembly,
                       _recombination_rate * _u[_qp](0),
                       _var.dofIndices()[0],
                       getVariable(_coupled_var_num).dofIndices()[0],
                       _var.arrayScalingFactor()[0]);
}

void
PointDefectRecombinationNodalKernel::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar != _coupled_var_num || _assemble_in_own_block)
    return;

  setJacobian(0, 0, _recombination_rate * _u[_qp](0));
}
