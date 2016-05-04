#include "../include/ode_argument.h"

int OdeArgument::setup_jacobian_prec(double,
				     Vector<double> const&,
				     Vector<double> const&,
				     double) 
{
  Assert(false, ExcPureFunctionCalled());
  return 0;
}

int OdeArgument::jacobian_prec(double, dealii::Vector<double>&,
			       dealii::Vector<double> const&,
			       dealii::Vector<double> const&,
			       dealii::Vector<double> const&,
			       double) 
{
  Assert(false, ExcPureFunctionCalled());
  return 0;
}

int OdeArgument::jacobian(double,
			  dealii::Vector<double>&,
			  dealii::Vector<double> const&,
			  dealii::Vector<double> const&,
			  dealii::Vector<double> const&,
			  double)
{
  Assert(false, ExcPureFunctionCalled());
  return 0;
}


Vector<double> & OdeArgument::differential_components() 
{
  Assert(false, ExcPureFunctionCalled());
  static Vector<double> tmp;
  return tmp;
}


Vector<double> & OdeArgument::get_diameters() 
{
  static Vector<double> tmp;
  if(tmp.size() == 0) {
     tmp.reinit(n_dofs());
     tmp = 1;
  }
  return tmp;
}
