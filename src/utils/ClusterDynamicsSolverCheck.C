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

#include "ClusterDynamicsSolverCheck.h"

#include "FEProblemBase.h"
#include "MooseObject.h"
#include "NonlinearSystemBase.h"

#include "libmesh/coupling_matrix.h"
#include "libmesh/dense_matrix.h"

#include <petscsnes.h>

#include <algorithm>
#include <mutex>
#include <set>

namespace
{
/// Value of the PETSc option name under prefix in the options database, or fallback if unset
std::string
petscOption(::PetscOptions options,
            const std::string & prefix,
            const std::string & name,
            const std::string & fallback)
{
  char value[256];
  PetscBool set = PETSC_FALSE;
  LibmeshPetscCallA(
      PETSC_COMM_SELF,
      PetscOptionsGetString(options, prefix.c_str(), name.c_str(), value, sizeof(value), &set));
  return set ? std::string(value) : fallback;
}

/// Estimated memory [bytes] above which the cross-variable storage warning is issued
constexpr double cross_variable_warning_bytes = 1.0e9;

/// Families of local Jacobian/residual block storage that MOOSE sizes by the square of the number
/// of components (_sub_Kee, _sub_Keg, _sub_Ken, _sub_Kne, _sub_Knn, _sub_Kll, _sub_Kle, _sub_Kln,
/// _sub_Kel and _sub_Knl in Assembly)
constexpr unsigned int block_families = 10;

/// True if a component of one variable is coupled to a component of a different variable
bool
hasCrossVariableCoupling(const libMesh::CouplingMatrix & cm, NonlinearSystemBase & nl)
{
  for (const auto * var : nl.getVariables(0))
  {
    const unsigned int first = var->number();
    const unsigned int last = first + var->count();
    for (unsigned int component = first; component < last; ++component)
      for (const auto other : libMesh::ConstCouplingRow(component, cm))
        if (other < first || other >= last)
          return true;
  }
  return false;
}
} // namespace

namespace ClusterDynamics
{
void
checkSolverSetup(const MooseObject & object,
                 const FEProblemBase & problem,
                 unsigned int nl_sys_num,
                 unsigned int first_component,
                 unsigned int n_components)
{
  if (problem.useHashTableMatrixAssembly() && problem.restoreOriginalNonzeroPattern())
    object.mooseWarning(
        "Set Problem/restore_original_nonzero_pattern = false. With hash table matrix assembly, "
        "MOOSE otherwise rebuilds the matrix from the hash table for every Jacobian, which "
        "costs time proportional to the square of the number of cluster sizes because of the "
        "dense monomer row. The cluster-dynamics sparsity pattern does not change, so the "
        "pattern from the first assembly can be kept.");

  const auto * cm = problem.couplingMatrix(nl_sys_num);
  if (cm && n_components > 1 && (*cm)(first_component, first_component + 1))
    object.mooseWarning(
        "The preconditioner couples the array components of the cluster variable with each "
        "other. This is the default when there is no Preconditioning block (an SMP "
        "preconditioner with full = true), and building and traversing that coupling costs "
        "time proportional to the square of the number of cluster sizes. The kernel adds its "
        "Jacobian entries directly, so use an SMP preconditioner with full = false:\n"
        "  [Preconditioning]\n    [smp]\n      type = SMP\n      full = false\n    []\n  []");

  // MOOSE allocates local block storage for every pair of components in the nonlinear system as
  // soon as two different variables are coupled (Assembly::init), so the memory grows with the
  // square of the total number of components of ALL variables, not only of the coupled ones.
  if (cm)
  {
    auto & nl = const_cast<FEProblemBase &>(problem).getNonlinearSystemBase(nl_sys_num);
    const double n_components_total = nl.nVariables();
    // The time tag is not yet registered when kernels are built, so assume at least two tags
    const double n_tags = std::max<double>(problem.numMatrixTags(), 2);
    const double bytes = block_families * n_tags * n_components_total * n_components_total *
                         sizeof(libMesh::DenseMatrix<double>);
    static std::mutex warned_mutex;
    static std::set<const FEProblemBase *> warned;
    if (bytes > cross_variable_warning_bytes && hasCrossVariableCoupling(*cm, nl))
    {
      bool first_time;
      {
        std::lock_guard<std::mutex> lock(warned_mutex);
        first_time = warned.insert(&problem).second;
      }
      if (first_time)
        object.mooseWarning(
            "The preconditioner couples different array variables, and the nonlinear system has ",
            static_cast<unsigned int>(n_components_total),
            " array components in total. MOOSE then allocates local block storage for every pair "
            "of components of the whole system, which needs about ",
            bytes / 1.0e9,
            " GB however few pairs are coupled and is paid before the first solve. Reduce the "
            "number of components (coarser grouping, a smaller maximum cluster size), or leave the "
            "off_diag_row/off_diag_column pair out of the preconditioner (SMP with full = false) "
            "and set assemble_opposite_jacobian_in_own_block = true, which keeps the exact "
            "cross-variable Jacobian. Leaving the pair out WITHOUT that setting drops the "
            "cross-variable Jacobian.");
    }
  }
}

void
checkOppositeCoupling(const MooseObject & object,
                      const FEProblemBase & problem,
                      unsigned int nl_sys_num,
                      unsigned int first_component,
                      unsigned int opposite_first_component,
                      bool assembled_in_own_block)
{
  if (assembled_in_own_block)
    return;
  const auto * cm = problem.couplingMatrix(nl_sys_num);
  // No coupling matrix at all means every variable is coupled to every other one
  if (!cm || (*cm)(first_component, opposite_first_component))
    return;
  object.mooseWarning(
      "This kernel has opposite_coupled_var set, but the preconditioner does not couple the two "
      "variables, so MOOSE never asks it for the cross-variable Jacobian and that block is "
      "missing. The Newton iteration is then approximate: in the cases tested it needed several "
      "times as many Newton iterations or time steps, and with a Fokker-Planck tail the results "
      "differed from those with the exact Jacobian. Either list the pair in the preconditioner "
      "(SMP, full = false, off_diag_row and off_diag_column naming the monomer components of the "
      "two arrays), or set assemble_opposite_jacobian_in_own_block = true and leave the pair "
      "out, which keeps the exact Jacobian without the memory cost of the pair.");
}

bool
checkFactorOrdering(const MooseObject & object, NonlinearSystemBase & nl)
{
  SNES snes = nl.getSNES();
  if (!snes)
    return false;

  KSP ksp;
  PC pc;
  PCType pc_type;
  LibmeshPetscCallA(PETSC_COMM_SELF, SNESGetKSP(snes, &ksp));
  LibmeshPetscCallA(PETSC_COMM_SELF, KSPGetPC(ksp, &pc));
  LibmeshPetscCallA(PETSC_COMM_SELF, PCGetType(pc, &pc_type));
  if (!pc_type)
    return false;

  ::PetscOptions options;
  const char * pc_prefix;
  LibmeshPetscCallA(PETSC_COMM_SELF, PetscObjectGetOptions((PetscObject)pc, &options));
  LibmeshPetscCallA(PETSC_COMM_SELF, PCGetOptionsPrefix(pc, &pc_prefix));
  std::string prefix = pc_prefix ? pc_prefix : "";
  std::string type = pc_type;

  // Block preconditioners factor each block with a sub preconditioner, ILU by default
  if (type == PCBJACOBI || type == PCASM || type == PCGASM)
  {
    type = petscOption(options, prefix, "-sub_pc_type", PCILU);
    prefix += "sub_";
  }

  const bool complete = type == PCLU || type == PCCHOLESKY;
  const bool incomplete = type == PCILU || type == PCICC;
  if (!complete && !incomplete)
    return true;

  // External factorization packages choose their own ordering
  if (petscOption(options, prefix, "-pc_factor_mat_solver_type", MATSOLVERPETSC) != MATSOLVERPETSC)
    return true;

  // PETSc orders complete factorizations with nested dissection and incomplete ones naturally
  const auto ordering = petscOption(options,
                                    prefix,
                                    "-pc_factor_mat_ordering_type",
                                    complete ? MATORDERINGND : MATORDERINGNATURAL);
  if (ordering != MATORDERINGNATURAL)
    return true;

  const std::string option = "-" + prefix + "pc_factor_mat_ordering_type";
  const std::string advice =
      "Set the PETSc option " + option +
      " to 'nd' or 'rcm', for example with petsc_options_iname = '-pc_factor_mat_ordering_type' "
      "and petsc_options_value = 'rcm' in the Executioner block.";
  if (complete)
    object.mooseError("The cluster-dynamics Jacobian has a dense monomer row and column. Its '",
                      type,
                      "' factorization with the natural ordering fills in to a dense matrix, "
                      "which needs memory proportional to the square of the number of cluster "
                      "sizes. ",
                      advice);
  else
  {
    const bool has_fill = petscOption(options, prefix, "-pc_factor_levels", "0") != "0";
    object.mooseWarning("The cluster-dynamics Jacobian has a dense monomer row and column. Its '",
                        type,
                        "' factorization with the natural ordering, which is the PETSc default "
                        "for incomplete factorizations, costs time",
                        has_fill ? " and memory" : "",
                        " proportional to the square of the number of cluster sizes. ",
                        advice);
  }
  return true;
}
} // namespace ClusterDynamics
