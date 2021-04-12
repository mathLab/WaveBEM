#ifndef newton_solver_h
#define newton_solver_h

#include <deal.II/base/logstream.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/parameter_handler.h>


#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#  include "mpi.h"
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif



#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"

#include "NOX.H"
#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "NOX_Epetra_Interface_Preconditioner.H"
#include "NOX_Epetra_LinearSystem_AztecOO.H"
#include "NOX_Epetra_Group.H"
#include "newton_argument.h"
#include <deal.II/lac/solver_control.h>







// ==========================================================================
// JacobianOperator is the class that implements the user provided jacobian
// as an EpetraOperator
// ==========================================================================
class JacobianOperator :

  public Epetra_Operator
{
public:

  // The constructor accepts an initial guess and the ode_argument class
  // containing the specific function to be zeroed.  We make
  // deep copies of each.
  JacobianOperator(NewtonArgument &Solver,
                   const Vector<double> &current_solution,
                   const Epetra_Map &Map,
                   const Epetra_Comm &Comm) :
    solver(Solver),
    y(current_solution),
    comm(Comm),
    map(Map)
  {
  }

  virtual ~JacobianOperator()
  {
  }

  int Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
  {

    Assert(X.Map().NumGlobalElements() == int(solver.n_dofs()),
           ExcDimensionMismatch(X.Map().NumGlobalElements(), solver.n_dofs()));
    Assert(Y.Map().NumGlobalElements() == int(solver.n_dofs()),
           ExcDimensionMismatch(Y.Map().NumGlobalElements(), solver.n_dofs()));
    //const VectorView<double> x(X.Map().NumGlobalElements(), &X[0][0]);
    //VectorView<double> dst(Y.Map().NumGlobalElements(), &Y[0][0]);
    
    Vector<double> x(solver.n_dofs());
    Vector<double> dst(solver.n_dofs());
    for (unsigned int i=0; i<solver.n_dofs(); ++i)
        x(i) = X[0][i];
    

    solver.jacobian(dst,y,x);
    
    for (unsigned int i=0; i<solver.n_dofs(); ++i)
        Y[0][i] = dst(i);
    
    
    //cout<<"X: "<<X[0][0]<<" "<<X[0][1]<<" (t)"<<endl;
    //cout<<"x: "<<x(0)<<" "<<x(1)<<" (d)"<<endl;
    //cout<<"Y: "<<Y[0][0]<<" "<<Y[0][1]<<" (t)"<<endl;
    //cout<<"dst: "<<dst(0)<<" "<<dst(1)<<" (d)"<<endl;
    //cout<<"y: "<<y(0)<<" "<<y(1)<<" (d)"<<endl;
    return 0;

  }

  void vmult(Vector<double> &dst, const Vector<double> &src) const
  {

    Assert(dst.size() == solver.n_dofs(),
           ExcDimensionMismatch(dst.size(), solver.n_dofs()));
    Assert(src.size() == solver.n_dofs(),
           ExcDimensionMismatch(src.size(), solver.n_dofs()));

    solver.jacobian(dst,y,src);

  }

  int ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
  {
    Assert(false, ExcPureFunctionCalled());
    return 0;
  }


  double NormInf() const
  {
    Assert(false, ExcPureFunctionCalled());
    return 0;
  }

  int SetUseTranspose(bool UseTranspose)
  {
    UseTranspose = false;
    return 0;
  }

  const char *Label() const
  {
    return ("jacobian_operator");
  }

  bool UseTranspose() const
  {
    return false;
  }

  bool HasNormInf() const
  {
    return false;
  }

  const Epetra_Comm &Comm() const
  {
    return comm;
  }

  const Epetra_Map &OperatorDomainMap() const
  {
    return map;
  }

  const Epetra_Map &OperatorRangeMap() const
  {
    return map;
  }



private:
  NewtonArgument &solver;
  const Vector<double> &y;
  const Epetra_Comm &comm;
  const Epetra_Map &map;

};


// ==========================================================================
// PreconditionerOperator is the class that implements the user provided jacobian
// as an EpetraOperator
// ==========================================================================
class PreconditionerOperator :

  public Epetra_Operator
{
public:

  // The constructor accepts an initial guess and the ode_argument class
  // containing the specific function to be zeroed.  We make
  // deep copies of each.
  PreconditionerOperator(NewtonArgument &Solver,
                         const Vector<double> &current_solution,
                         const Epetra_Map &Map,
                         const Epetra_Comm &Comm) :
    solver(Solver),
    y(current_solution),
    comm(Comm),
    map(Map)
  {
  }

  virtual ~PreconditionerOperator()
  {
  }

  int Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
  {
    Assert(X.Map().NumGlobalElements() == int(solver.n_dofs()),
           ExcDimensionMismatch(X.Map().NumGlobalElements(), solver.n_dofs()));
    Assert(Y.Map().NumGlobalElements() == int(solver.n_dofs()),
           ExcDimensionMismatch(Y.Map().NumGlobalElements(), solver.n_dofs()));
    //const VectorView<double> x(X.Map().NumGlobalElements(), &X[0][0]);
    //VectorView<double> dst(Y.Map().NumGlobalElements(), &Y[0][0]);


    Vector<double> x(solver.n_dofs());
    Vector<double> dst(solver.n_dofs());
    for (unsigned int i=0; i<solver.n_dofs(); ++i)
        x(i) = X[0][i];
    

    solver.jacobian_prec_prod(dst,y,x);
    
    for (unsigned int i=0; i<solver.n_dofs(); ++i)
        Y[0][i] = dst(i);


    

    return 0;

  }


  int ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
  {
    //cout<<"X: "<<X[0][0]<<" "<<X[0][1]<<" (t)"<<endl;
    //cout<<"Y: "<<Y[0][0]<<" "<<Y[0][1]<<" (t)"<<endl;

    //Y[0][0] = 5.0;

    //cout<<"*X: "<<X[0][0]<<" "<<X[0][1]<<" (t)"<<endl;
    //cout<<"*Y: "<<Y[0][0]<<" "<<Y[0][1]<<" (t)"<<endl;

    Assert(X.Map().NumGlobalElements() == int(solver.n_dofs()),
           ExcDimensionMismatch(X.Map().NumGlobalElements(), solver.n_dofs()));
    Assert(Y.Map().NumGlobalElements() == int(solver.n_dofs()),
           ExcDimensionMismatch(Y.Map().NumGlobalElements(), solver.n_dofs()));
    //const VectorView<double> x(X.Map().NumGlobalElements(), &X[0][0]);
    //VectorView<double> dst(Y.Map().NumGlobalElements(), &Y[0][0]);

    Vector<double> x(solver.n_dofs());
    Vector<double> dst(solver.n_dofs());
    for (unsigned int i=0; i<solver.n_dofs(); ++i)
        x(i) = X[0][i];
    

    solver.jacobian_prec(dst,y,x);
    
    for (unsigned int i=0; i<solver.n_dofs(); ++i)
        Y[0][i] = dst(i);




    //cout<<"X: "<<X[0][0]<<" "<<X[0][1]<<" (t)"<<endl;
    //cout<<"x: "<<x(0)<<" "<<x(1)<<" (d)"<<endl;
    //cout<<"Y: "<<Y[0][0]<<" "<<Y[0][1]<<" (t)"<<endl;
    //cout<<"dst: "<<dst(0)<<" "<<dst(1)<<" (d)"<<endl;
    //cout<<"y: "<<y(0)<<" "<<y(1)<<" (d)"<<endl;
    return 0;
  }

  void vmult(Vector<double> &dst, const Vector<double> &src) const
  {

    Assert(dst.size() == solver.n_dofs(),
           ExcDimensionMismatch(dst.size(), solver.n_dofs()));
    Assert(src.size() == solver.n_dofs(),
           ExcDimensionMismatch(src.size(), solver.n_dofs()));

    solver.jacobian_prec(dst,y,src);

  }


  double NormInf() const
  {
    Assert(false, ExcPureFunctionCalled());
    return 0;
  }

  int SetUseTranspose(bool UseTranspose)
  {
    UseTranspose = false;
    return 0;
  }

  const char *Label() const
  {
    return ("preconditioner_operator");
  }

  bool UseTranspose() const
  {
    return false;
  }

  bool HasNormInf() const
  {
    return false;
  }

  const Epetra_Comm &Comm() const
  {
    return comm;
  }

  const Epetra_Map &OperatorDomainMap() const
  {
    return map;
  }

  const Epetra_Map &OperatorRangeMap() const
  {
    return map;
  }



private:
  NewtonArgument &solver;
  const Vector<double> &y;
  const Epetra_Comm &comm;
  const Epetra_Map &map;


};



// ==========================================================================
// ProblemInterface, the problem interface in this example,
// defines the interface between NOX and our nonlinear problem to
// solve.
// ==========================================================================
class ProblemInterface :
  public NOX::Epetra::Interface::Required,
  public NOX::Epetra::Interface::Jacobian,
  public NOX::Epetra::Interface::Preconditioner
{
public:

  // The constructor accepts an initial guess and the ode_argument class
  // containing the specific function to be zeroed.  We make
  // deep copies of each.
  ProblemInterface (NewtonArgument &Solver,
                    Vector<double> &current_sol,
                    PreconditionerOperator &prec_oper,
                    JacobianOperator &jac_oper) :
    solver(Solver),
    y(current_sol),
    preconditioner_operator(prec_oper),
    jacobian_operator(jac_oper)
  {

  }

  virtual ~ProblemInterface()
  {
  }

  // Compute f := F(x), where x is the input vector and f the output
  // vector.
  bool
  computeF (const Epetra_Vector &x,
            Epetra_Vector &f,
            NOX::Epetra::Interface::Required::FillType F)
  {

    Assert(x.Map().NumGlobalElements() == int(solver.n_dofs()),
           ExcDimensionMismatch(x.Map().NumGlobalElements(), solver.n_dofs()));
    Assert(f.Map().NumGlobalElements() == int(solver.n_dofs()),
           ExcDimensionMismatch(f.Map().NumGlobalElements(), solver.n_dofs()));


    //const VectorView<double> yy(x.Map().NumGlobalElements(), &x[0]);
    //VectorView<double> residual(f.Map().NumGlobalElements(), &f[0]);
    
    
    Vector<double> yy(solver.n_dofs());
    Vector<double> residual(solver.n_dofs());
    for (unsigned int i=0; i<solver.n_dofs(); ++i)
        yy(i) = x[i];
    

    solver.residual(residual, yy);
    
    for (unsigned int i=0; i<solver.n_dofs(); ++i)
        f[i] = residual(i);
    
    return true;
  };

  bool
  computeJacobian(const Epetra_Vector &x, Epetra_Operator &Jac)
  {

    Assert(x.Map().NumGlobalElements() == int(solver.n_dofs()),
           ExcDimensionMismatch(x.Map().NumGlobalElements(), solver.n_dofs()));

    //const VectorView<double> yy(x.Map().NumGlobalElements(), &x[0]);
    
    for (unsigned int i=0; i<solver.n_dofs(); ++i)
        y(i) = x[i];
    
    Jac = jacobian_operator;
    //Epetra_Vector y(x);
    //cout<<"Test "<<endl,
    //Jac.Apply(x,y);

    return true;
  }



  bool computePrecMatrix (const Epetra_Vector &x,
                          Epetra_RowMatrix &M)
  {
    throw std::runtime_error ("*** ProblemInterface does not implement "
                              "computing an explicit preconditioner from an "
                              "Epetra_RowMatrix ***");
  }

  bool computePreconditioner (const Epetra_Vector &x,
                              Epetra_Operator &Prec,
                              Teuchos::ParameterList * =0)
  {
    Assert(x.Map().NumGlobalElements() == int(solver.n_dofs()),
           ExcDimensionMismatch(x.Map().NumGlobalElements(), solver.n_dofs()));

    //const VectorView<double> yy(x.Map().NumGlobalElements(), &x[0]);
    Vector<double> yy(solver.n_dofs());
    for (unsigned int i=0; i<solver.n_dofs(); ++i)
        yy(i) = x[i];
    
    solver.setup_jacobian_prec(yy);
    y = yy;
    Prec = preconditioner_operator;

    return true;
  }



private:

  NewtonArgument &solver;
  Vector<double> &y;
  PreconditionerOperator &preconditioner_operator;
  JacobianOperator &jacobian_operator;

};

class NewtonSolver
{
public:
  /** Constructor for the NewtonSolver class. The Solver class is
   * required to have a Solver.solve(Vector<double> &dst, const
   * Vector<double> &src) method that will be called by the time
   * integrator to find out about the solution to a given src. */
  NewtonSolver(NewtonArgument &solver);


  /** House cleaning. */
  ~NewtonSolver();

  /** Declare parameters for this class to function properly. */
  static void declare_parameters(ParameterHandler &prm);

  /** Parse a parameter handler. */
  void parse_parameters(ParameterHandler &prm);

  /** Solve. This function returns the final number of steps. */
  unsigned int solve(Vector<double> &solution,
                     const unsigned int max_steps);




private:
  /** The bubble membrane poperties. */
  NewtonArgument &solver;
  Vector<double> y;
#ifdef EPETRA_MPI
  Epetra_MpiComm comm;
#else
  Epetra_SerialComm comm;
#endif
  Epetra_Map map;
  PreconditionerOperator *preconditioner_operator;
  JacobianOperator *jacobian_operator;
  std::string linear_solver_name;
  bool provide_jac;
  bool provide_jac_prec;
  double rel_tol;
  double linear_rel_tol;
  SolverControl solver_control;


};

#endif
