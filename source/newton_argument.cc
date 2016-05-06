#include "../include/newton_argument.h"

int NewtonArgument::setup_jacobian_prec(Vector<double> const &)
{
  Assert(false, ExcPureFunctionCalled());
  return 0;
}

int NewtonArgument::jacobian_prec(dealii::Vector<double> &,
                                  dealii::Vector<double> const &,
                                  dealii::Vector<double> const &)
{
  Assert(false, ExcPureFunctionCalled());
  return 0;
}

int NewtonArgument::jacobian_prec_prod(dealii::Vector<double> &,
                                       dealii::Vector<double> const &,
                                       dealii::Vector<double> const &)
{
  Assert(false, ExcPureFunctionCalled());
  return 0;
}

int NewtonArgument::jacobian(dealii::Vector<double> &,
                             dealii::Vector<double> const &,
                             dealii::Vector<double> const &)
{
  Assert(false, ExcPureFunctionCalled());
  return 0;
}





