#include "../include/newton_solver.h"
#include "../include/newton_argument.h"
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>


#include <deal.II/base/utilities.h>
#include <iostream>

using namespace dealii;
using namespace std;



NewtonSolver::NewtonSolver(NewtonArgument &bubble) :
    solver(bubble),
    y(solver.n_dofs()),
    comm(MPI_COMM_WORLD),
    map((int)solver.n_dofs(), 0, comm)
    
{
    preconditioner_operator = new PreconditionerOperator(solver,y,map,comm);
    jacobian_operator = new JacobianOperator(solver,y,map,comm);
    // assigning default values to solver options
    linear_solver_name = "GMRES";
    provide_jac = true;
    provide_jac_prec = true;
    linear_rel_tol = 1e-8;
    rel_tol = 1e-7; 
    

}


NewtonSolver::~NewtonSolver()
{
}


void NewtonSolver::declare_parameters(ParameterHandler &prm) {
    prm.enter_subsection("Newton Solver Params");
    {
    prm.declare_entry("Linear solver kind", "GMRES", Patterns::Selection("GMRES|CG|CGS|TFQMR|BiCGStab|LU"));
    prm.declare_entry("Provide jacobian product", "false", Patterns::Bool());
    prm.declare_entry("Provide jacobian preconditioner", "false", Patterns::Bool());
    prm.declare_entry("Relative nonlinear error tolerance", "1e-5", Patterns::Double());
    prm.declare_entry("Relative linear error tolerance", "1e-8", Patterns::Double());
    }
    prm.leave_subsection();
    prm.enter_subsection("Solver");
    SolverControl::declare_parameters(prm);
    prm.leave_subsection();
}


void NewtonSolver::parse_parameters(ParameterHandler &prm) {
    prm.enter_subsection("Newton Solver Params");
    {
    linear_solver_name = prm.get("Linear solver kind");
    provide_jac   = prm.get_bool("Provide jacobian product");
    provide_jac_prec = prm.get_bool("Provide jacobian preconditioner");
    rel_tol = prm.get_double("Relative nonlinear error tolerance");
    linear_rel_tol = prm.get_double("Relative linear error tolerance");
    }
    prm.leave_subsection();
    prm.enter_subsection("Solver");
    solver_control.parse_parameters(prm);
    prm.leave_subsection();
}

unsigned int NewtonSolver::solve(Vector<double> &solution, 
                                 const unsigned int max_steps)
{     


    AssertThrow(solution.size() == solver.n_dofs(),
		ExcDimensionMismatch(solution.size(), solver.n_dofs()));




    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;


    //Let's resize the Epetra Map (if the solution)
    map = Epetra_Map((int)solver.n_dofs(), 0, comm);
    Epetra_Vector sol(map);
    Epetra_Vector init_guess(map);

    y = solution;
    for (unsigned int i=0; i<solver.n_dofs(); ++i)
        {
        init_guess[i] = solution(i);
        sol[i] = solution(i);
        }



    // Set up the problem interface.  Your application will define
    // its own problem interface.  ProblemInterface is our
    // example interface, which you can use as a model.
    // 
    // Our ProblemInterface makes a deep copy of the initial
    // guess.
    Teuchos::RCP<ProblemInterface> interface = 
    Teuchos::rcp(new ProblemInterface(solver, y, *preconditioner_operator, *jacobian_operator) );

    RCP<Epetra_Operator> preconditioner(preconditioner_operator);
    RCP<Epetra_Operator> jacobian(jacobian_operator);
    /*
    cout<<"========================================================="<<endl;
    cout<<"Deal test"<<endl;
    SolverGMRES<Vector<double> > gmres(solver_control,
         SolverGMRES<Vector<double> >::AdditionalData(100));

    Vector<double> res(solver.n_dofs());
    solver.residual(res,solution);
    res *= -1.0;
    
    gmres.solve (*(interface->jacobian_operator), solution, res, PreconditionIdentity());// *(interface->preconditioner_operator));
    cout<<"Passed"<<endl;
    cout<<"========================================================="<<endl;
    */





    // Create the top-level parameter list to control NOX.
    //
    // "parameterList" (lowercase initial "p") is a "nonmember
    // constructor" that returns an RCP<ParameterList> with the
    // given name.
    RCP<ParameterList> params = parameterList ("NOX");

    // Tell the nonlinear solver to use line search.
    params->set ("Nonlinear Solver", "Line Search Based");

    //
    // Set the printing parameters in the "Printing" sublist.
    //
    ParameterList& printParams = params->sublist ("Printing");
    printParams.set ("MyPID", comm.MyPID ()); 
    printParams.set ("Output Precision", 3);
    printParams.set ("Output Processor", 0);

    // Set verbose=true to see a whole lot of intermediate status
    // output, during both linear and nonlinear iterations.
    const bool verbose = false;
    if (verbose) {
        printParams.set ("Output Information", 
                         NOX::Utils::OuterIteration + 
                         NOX::Utils::OuterIterationStatusTest + 
                         NOX::Utils::InnerIteration +
                         NOX::Utils::Parameters + 
                         NOX::Utils::Details + 
                         NOX::Utils::Warning);
      } else {
        printParams.set ("Output Information", NOX::Utils::Warning);
      }

      //
      // Set the nonlinear solver parameters.
      //

      // Line search parameters.
      ParameterList& searchParams = params->sublist ("Line Search");
      searchParams.set ("Method", "Full Step");

      // Parameters for picking the search direction.
      ParameterList& dirParams = params->sublist ("Direction");
      // Use Newton's method to pick the search direction.
      dirParams.set ("Method", "Newton");

      // Parameters for Newton's method.
      ParameterList& newtonParams = dirParams.sublist ("Newton");
      newtonParams.set ("Forcing Term Method", "Constant");

      //
      // Newton's method invokes a linear solver repeatedly.
      // Set the parameters for the linear solver.
      //
      ParameterList& lsParams = newtonParams.sublist ("Linear Solver");

      // Use Aztec's implementation of GMRES, with at most 800
      // iterations, a residual tolerance of 1.0e-4, with output every
      // 50 iterations, and Aztec's native ILU preconditioner.
      lsParams.set ("Aztec Solver", linear_solver_name);  
      lsParams.set ("Max Iterations", 800);  
      lsParams.set ("Tolerance", linear_rel_tol);
      lsParams.set ("Output Frequency", 1);
      lsParams.set ("Preconditioner Operator", "Finite Difference");
      lsParams.set ("Aztec Preconditioner", "ilu"); 
      //lsParams.set("Preconditioner", "User Defined");
      lsParams.set("Preconditioner", "User Defined");
      //lsParams.set("Preconditioner Reuse Policy", "Reuse");


      // Our ProblemInterface implements both Required and
      // Jacobian, so we can use the same object for each.
      RCP<NOX::Epetra::Interface::Required> iReq = interface;
      RCP<NOX::Epetra::Interface::Jacobian> iJac = interface;
      RCP<NOX::Epetra::Interface::Preconditioner> iPrec = interface;


      RCP<NOX::Epetra::LinearSystemAztecOO> linSys;

      if ((provide_jac==false) && (provide_jac_prec==false) )
         linSys = rcp (new NOX::Epetra::LinearSystemAztecOO (printParams, lsParams,
                                                             iReq, init_guess));
      else if ((provide_jac==true) && (provide_jac_prec==false) )
         linSys = rcp (new NOX::Epetra::LinearSystemAztecOO (printParams, lsParams,
                                                             iReq, iJac, jacobian, init_guess));
      else if ((provide_jac==false) && (provide_jac_prec==true) )
         linSys = rcp (new NOX::Epetra::LinearSystemAztecOO (printParams, lsParams,
                                                             iReq, iPrec, preconditioner, init_guess));
      else if ((provide_jac==true) && (provide_jac_prec==true) )
         {
         linSys = rcp (new NOX::Epetra::LinearSystemAztecOO (printParams, lsParams,
                                                             iJac, jacobian, iPrec, preconditioner, init_guess)); 
         }

      // Need a NOX::Epetra::Vector for constructor.
      NOX::Epetra::Vector noxInitGuess (init_guess, NOX::DeepCopy);
      

      RCP<NOX::Epetra::Group> group = 
        rcp (new NOX::Epetra::Group (printParams, iReq, noxInitGuess, linSys));

      //
      // Set up NOX's iteration stopping criteria ("status tests").
      //

      // ||F(X)||_2 / N < 1.0e-4, where N is the length of F(X).
      //
      // NormF has many options for setting up absolute vs. relative
      // (scaled by the norm of the initial guess) tolerances, scaling
      // or not scaling by the length of F(X), and choosing a
      // different norm (we use the 2-norm here).
      RCP<NOX::StatusTest::NormF> testNormF = 
        rcp (new NOX::StatusTest::NormF (rel_tol));

      // At most 20 (nonlinear) iterations.
      RCP<NOX::StatusTest::MaxIters> testMaxIters = 
        rcp (new NOX::StatusTest::MaxIters (max_steps));

      // Combine the above two stopping criteria (normwise
      // convergence, and maximum number of nonlinear iterations).
      // The result tells NOX to stop if at least one of them is
      // satisfied.
      RCP<NOX::StatusTest::Combo> combo = 
        rcp (new NOX::StatusTest::Combo (NOX::StatusTest::Combo::OR, 
                                         testNormF, testMaxIters));

      // Create the NOX nonlinear solver.
      RCP<NOX::Solver::Generic> solver = 
        NOX::Solver::buildSolver (group, combo, params);

      // Solve the nonlinear system.
      NOX::StatusTest::StatusType status = solver->solve();

      // Print the result.
      //
      // For this particular example, Comm contains only one MPI
      // process.  However, we check for Comm.MyPID() == 0 here just
      // so that the example is fully general.  (If you're solving a
      // larger nonlinear problem, you could safely use the code
      // below.)
      if (comm.MyPID() == 0) {
        cout << endl << "-- Parameter List From Solver --" << endl;
        solver->getList ().print (cout);
      }

      // Get the Epetra_Vector with the final solution from the solver.
      const NOX::Epetra::Group& finalGroup = 
        dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());

      const Epetra_Vector& finalSolution = 
        dynamic_cast<const NOX::Epetra::Vector&> (finalGroup.getX ()).getEpetraVector ();

      for (unsigned int i=0; i<solution.size(); ++i)
          {
          solution(i) = finalSolution[i];
          }      

      //if (Comm.MyPID() == 0) {
        //cout << "Computed solution: " << endl;
      //}
      // Epetra objects know how to print themselves politely when
      // their operator<<(std::ostream&) is invoked on all MPI
      // process(es) in the communicator to which they are associated.
      //cout << finalSolution;

    return 0;
}



