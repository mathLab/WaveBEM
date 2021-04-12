#include "../include/dae_time_integrator.h"
#include "../include/ode_argument.h"
#include <deal.II/sundials/copy.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/array_view.h>
#include <iostream>

using namespace dealii;
using namespace std;

int dae_residual(realtype tt, N_Vector yy, N_Vector yp,
                 N_Vector rr, void *user_data)
{
  OdeArgument &solver = *static_cast<OdeArgument *>(user_data);
  Assert(NV_LENGTH_S(yy) == solver.n_dofs(),
         ExcDimensionMismatch(NV_LENGTH_S(yy), solver.n_dofs()));
  Assert(NV_LENGTH_S(yp) == solver.n_dofs(),
         ExcDimensionMismatch(NV_LENGTH_S(yp), solver.n_dofs()));
//    if ((NV_LENGTH_S(yy) == NV_LENGTH_S(yp)) & (NV_LENGTH_S(yy) =! NV_LENGTH_S(rr)))
  NV_LENGTH_S(rr) = NV_LENGTH_S(yy);

  Assert(NV_LENGTH_S(rr) == solver.n_dofs(),
         ExcDimensionMismatch(NV_LENGTH_S(rr), solver.n_dofs()));
  

  GrowingVectorMemory<Vector<double> > mem;

  VectorMemory<Vector<double> >::Pointer src_yy(mem);
  src_yy->reinit(solver.n_dofs());

  VectorMemory<Vector<double> >::Pointer src_yp(mem);
  src_yp->reinit(solver.n_dofs());

  VectorMemory<Vector<double> >::Pointer residual(mem);
  residual->reinit(solver.n_dofs());

  for (std::size_t i = 0; i < solver.n_dofs(); ++i)
      {
      (*src_yy)(i) = NV_Ith_S(yy, i);
      (*src_yp)(i) = NV_Ith_S(yp, i);
      }
  
  int err = solver.residual(tt, *residual, *src_yy, *src_yp);

  for (std::size_t i = 0; i < solver.n_dofs(); ++i)
      {
      NV_Ith_S(rr, i) = (*residual)(i);
      }
  
  return err;


//  const ArrayView<double> src_yy(solver.n_dofs(), NV_DATA_S(yy));
//  const ArrayView<double> src_yp(solver.n_dofs(), NV_DATA_S(yp));
//  ArrayView<double> residual(solver.n_dofs(), NV_DATA_S(rr));
//  return solver.residual(tt, residual, src_yy, src_yp);
}

int dae_setup_jac(realtype tt,
                  N_Vector yy,
                  N_Vector yp,
                  N_Vector rr,
                  realtype cj,
                  void *user_data)
{  
return 0;
}


int dae_setup_prec(realtype tt, // time
                   N_Vector yy,
                   N_Vector yp,
                   N_Vector rr, // Current residual
                   realtype alpha, // J = dG/dyy + alpha dG/dyp
                   void *user_data)//, // the pointer to the correct class
                   //N_Vector /*tmp1*/, // temporary storage
                   //N_Vector /*tmp2*/,
                   //N_Vector /*tmp3*/)
{
  OdeArgument &solver = *static_cast<OdeArgument *>(user_data);

  Assert(NV_LENGTH_S(yy) == solver.n_dofs(),
         ExcDimensionMismatch(NV_LENGTH_S(yy), solver.n_dofs()));
  Assert(NV_LENGTH_S(yp) == solver.n_dofs(),
         ExcDimensionMismatch(NV_LENGTH_S(yp), solver.n_dofs()));
  Assert(NV_LENGTH_S(rr) == solver.n_dofs(),
         ExcDimensionMismatch(NV_LENGTH_S(rr), solver.n_dofs()));

  GrowingVectorMemory<Vector<double> > mem;

  VectorMemory<Vector<double> >::Pointer src_yy(mem);
  src_yy->reinit(solver.n_dofs());

  VectorMemory<Vector<double> >::Pointer src_yp(mem);
  src_yp->reinit(solver.n_dofs());

  VectorMemory<Vector<double> >::Pointer residual(mem);
  residual->reinit(solver.n_dofs());

  for (std::size_t i = 0; i < solver.n_dofs(); ++i)
      {
      (*src_yy)(i) = NV_Ith_S(yy, i);
      (*src_yp)(i) = NV_Ith_S(yp, i);
      }

  return solver.setup_jacobian_prec(tt, *src_yy, *src_yp, alpha);

//  // A previous call to residual has already been done.
//  const ArrayView<double> src_yy(solver.n_dofs(), NV_DATA_S(yy));
//  const ArrayView<double> src_yp(solver.n_dofs(), NV_DATA_S(yp));
//  const ArrayView<double> residual(solver.n_dofs(), NV_DATA_S(rr));
//  return solver.setup_jacobian_prec(tt, src_yy, src_yp, alpha);
}

int dae_jtimes(realtype tt, N_Vector yy, N_Vector yp,
               N_Vector rr, // Current residual
               N_Vector src, // right hand side to solve for
               N_Vector dst, // computed output
               realtype alpha, // J = dG/dyy + alpha dG/dyp
               void *user_data, // the pointer to the correct class
               N_Vector /*tmp*/,
               N_Vector /*tmp2*/) // Storage
{
  OdeArgument &solver = *static_cast<OdeArgument *>(user_data);

  Assert(NV_LENGTH_S(yy) == solver.n_dofs(),
         ExcDimensionMismatch(NV_LENGTH_S(yy), solver.n_dofs()));
  Assert(NV_LENGTH_S(yp) == solver.n_dofs(),
         ExcDimensionMismatch(NV_LENGTH_S(yp), solver.n_dofs()));
  Assert(NV_LENGTH_S(src) == solver.n_dofs(),
         ExcDimensionMismatch(NV_LENGTH_S(src), solver.n_dofs()));
  Assert(NV_LENGTH_S(dst) == solver.n_dofs(),
         ExcDimensionMismatch(NV_LENGTH_S(dst), solver.n_dofs()));

  GrowingVectorMemory<Vector<double> > mem;

  VectorMemory<Vector<double> >::Pointer src_yy(mem);
  src_yy->reinit(solver.n_dofs());

  VectorMemory<Vector<double> >::Pointer src_yp(mem);
  src_yp->reinit(solver.n_dofs());

  VectorMemory<Vector<double> >::Pointer residual(mem);
  residual->reinit(solver.n_dofs());
  
  VectorMemory<Vector<double> >::Pointer src_v(mem);
  src_v->reinit(solver.n_dofs());

  VectorMemory<Vector<double> >::Pointer dst_v(mem);
  dst_v->reinit(solver.n_dofs());

  for (std::size_t i = 0; i < solver.n_dofs(); ++i)
      {
      (*src_yy)(i) = NV_Ith_S(yy, i);
      (*src_yp)(i) = NV_Ith_S(yp, i);
      (*src_v)(i) = NV_Ith_S(src, i);
      }
   int err = solver.jacobian(tt, *dst_v, *src_yy, *src_yp, *src_v, alpha);
   
   for (std::size_t i = 0; i < solver.n_dofs(); ++i)
      {
      NV_Ith_S(rr, i) = (*residual)(i);
      NV_Ith_S(dst, i) = (*dst_v)(i);
      }
      
   return err;
//  //double* a = NV_DATA_S(src);
//  //for (unsigned int i=0; i<NV_LENGTH_S(src); ++i)
//  //    if (!(numbers::is_finite(a[i])))
//  //       cout<<i<<endl;
//  OdeArgument &solver = *static_cast<OdeArgument *>(user_data);


//  // A previous call to residual has already been done.
//  const ArrayView<double> src_yy(solver.n_dofs(), NV_DATA_S(yy));
//  const ArrayView<double> src_yp(solver.n_dofs(), NV_DATA_S(yp));
//  const ArrayView<double> residual(solver.n_dofs(), NV_DATA_S(rr));
//  const ArrayView<double> src_v(solver.n_dofs(), NV_DATA_S(src));
//  VectorView<double> dst_v(solver.n_dofs(), NV_DATA_S(dst));
//  return solver.jacobian(tt, dst_v, src_yy, src_yp, src_v, alpha);
}


int dae_prec(realtype tt, N_Vector yy, N_Vector yp,
             N_Vector rr, // Current residual
             N_Vector rvec, // right hand side to solve for
             N_Vector zvec, // computed output
             realtype alpha, // J = dG/dyy + alpha dG/dyp
             realtype /*delta*/, // input tolerance. The residual rr - Pz has to be smaller than delta
             void *user_data) // the pointer to the correct class
             //N_Vector /*tmp*/) // Storage
{
  OdeArgument &solver = *static_cast<OdeArgument *>(user_data);
  Assert(NV_LENGTH_S(yy) == solver.n_dofs(),
         ExcDimensionMismatch(NV_LENGTH_S(yy), solver.n_dofs()));
  Assert(NV_LENGTH_S(yp) == solver.n_dofs(),
         ExcDimensionMismatch(NV_LENGTH_S(yp), solver.n_dofs()));
  Assert(NV_LENGTH_S(rr) == solver.n_dofs(),
         ExcDimensionMismatch(NV_LENGTH_S(rr), solver.n_dofs()));
  // A previous call to residual has already been done.
  
    GrowingVectorMemory<Vector<double> > mem;

  VectorMemory<Vector<double> >::Pointer src_yy(mem);
  src_yy->reinit(solver.n_dofs());

  VectorMemory<Vector<double> >::Pointer src_yp(mem);
  src_yp->reinit(solver.n_dofs());

  VectorMemory<Vector<double> >::Pointer residual(mem);
  residual->reinit(solver.n_dofs());
  
  VectorMemory<Vector<double> >::Pointer rhs(mem);
  rhs->reinit(solver.n_dofs());

  VectorMemory<Vector<double> >::Pointer output(mem);
  output->reinit(solver.n_dofs());

  for (std::size_t i = 0; i < solver.n_dofs(); ++i)
      {
      (*src_yy)(i) = NV_Ith_S(yy, i);
      (*src_yp)(i) = NV_Ith_S(yp, i);
      (*rhs)(i) = NV_Ith_S(rvec, i);
      }
  
  int err = solver.jacobian_prec(tt, *output, *src_yy, *src_yp, *rhs, alpha);
  
   for (std::size_t i = 0; i < solver.n_dofs(); ++i)
      {
      NV_Ith_S(zvec, i) = (*output)(i);
      }
   return err;
//  const VectorView<double> src_yy(solver.n_dofs(), NV_DATA_S(yy));
//  const VectorView<double> src_yp(solver.n_dofs(), NV_DATA_S(yp));
//  const VectorView<double> residual(solver.n_dofs(), NV_DATA_S(rr));
//  const VectorView<double> rhs(solver.n_dofs(), NV_DATA_S(rvec));
//  VectorView<double> output(solver.n_dofs(), NV_DATA_S(zvec));
//  return solver.jacobian_prec(tt, output, src_yy, src_yp, rhs, alpha);
}



DAETimeIntegrator::DAETimeIntegrator(OdeArgument &bubble) :
  solver(bubble),
  is_initialized(false)
{
  initial_step_size = 1e-4;
  min_step_size = 1e-6;

  abs_tol = 1e-6;
  rel_tol = 1e-8;

  ida_mem = IDACreate();
  is_initialized = true;

}


DAETimeIntegrator::DAETimeIntegrator(OdeArgument &bubble,
                                     double initial_step_size,
                                     double min_step_size,
                                     double initial_time,
                                     double final_time,
                                     double abs_tol,
                                     double rel_tol,
                                     unsigned int max_n_steps,
                                     double outputs_period,
                                     unsigned int ic_type,
                                     bool use_iterative,
                                     bool provide_jac,
                                     bool provide_jac_prec) :
  final_time(final_time),
  initial_time(initial_time),
  solver(bubble),
  initial_step_size(initial_step_size),
  min_step_size(min_step_size),
  abs_tol(abs_tol),
  rel_tol(rel_tol),
  max_n_steps(max_n_steps),
  outputs_period(outputs_period),
  ic_type(ic_type),
  is_initialized(false),
  use_iterative(use_iterative),
  provide_jac(provide_jac),
  provide_jac_prec(provide_jac_prec)
{

  ida_mem = IDACreate();
  is_initialized = true;

}



DAETimeIntegrator::~DAETimeIntegrator()
{
  if (ida_mem)
    IDAFree(&ida_mem);
}


void DAETimeIntegrator::declare_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("IDA Solver Params");
  {
    prm.declare_entry("Use iterative algorithm", "false", Patterns::Bool());
    prm.declare_entry("Provide jacobian product", "false", Patterns::Bool());
    prm.declare_entry("Provide jacobian preconditioner", "false", Patterns::Bool());
    prm.declare_entry("Initial step size", "1e-4", Patterns::Double());
    prm.declare_entry("Min step size", "5e-5", Patterns::Double());
    prm.declare_entry("Absolute error tolerance", "1e-4", Patterns::Double());
    prm.declare_entry("Relative error tolerance", "1e-3", Patterns::Double());
    prm.declare_entry("Initial time", "0", Patterns::Double());
    prm.declare_entry("Final time", "100000", Patterns::Double());
    prm.declare_entry("Seconds between each output", "1e-1", Patterns::Double());
    prm.declare_entry("Initial condition type", "0", Patterns::Integer());
  }
  prm.leave_subsection();
}


void DAETimeIntegrator::parse_parameters(ParameterHandler &prm)
{
  prm.enter_subsection("IDA Solver Params");
  {
    use_iterative = prm.get_bool("Use iterative algorithm");
    provide_jac   = prm.get_bool("Provide jacobian product");
    provide_jac_prec = prm.get_bool("Provide jacobian preconditioner");
    initial_step_size = prm.get_double("Initial step size");
    min_step_size     = prm.get_double("Min step size");
    abs_tol = prm.get_double("Absolute error tolerance");
    rel_tol = prm.get_double("Relative error tolerance");
    initial_time = prm.get_double("Initial time");
    final_time = prm.get_double("Final time");
    outputs_period = prm.get_double("Seconds between each output");
    ic_type = prm.get_integer("Initial condition type");
  }
  prm.leave_subsection();
}

unsigned int DAETimeIntegrator::start_ode(Vector<double> &solution,
                                          Vector<double> &solution_dot,
                                          const unsigned int max_steps)
{


  AssertThrow(solution.size() == solver.n_dofs(),
              ExcDimensionMismatch(solution.size(), solver.n_dofs()));

  AssertThrow(is_initialized, ExcMessage("Not Initialized!"));

  double t = initial_time;
  double h = initial_step_size;
  unsigned int step_number = 0;

  int status;

  // The solution is stored in
  // solution. Here we take only a
  // view of it.

  yy = N_VNewEmpty_Serial(solver.n_dofs());
  NV_DATA_S(yy) = solution.begin();

  yp = N_VNewEmpty_Serial(solver.n_dofs());
  NV_DATA_S(yp) = solution_dot.begin();

  diff_id = N_VNewEmpty_Serial(solver.n_dofs());
  NV_DATA_S(diff_id) = solver.differential_components().begin();

  Vector<double> tolerances = solver.get_diameters();
  tolerances*=(1/tolerances.linfty_norm()*abs_tol);
  abs_tolls = N_VNewEmpty_Serial(solver.n_dofs());
  NV_DATA_S(abs_tolls) = tolerances.begin();

  status = IDAInit(ida_mem, dae_residual, initial_time, yy, yp);
  //status += IDASStolerances(ida_mem, rel_tol, abs_tol);
  status += IDASVtolerances(ida_mem, rel_tol, abs_tolls);
  status += IDASetInitStep(ida_mem, step_number);
  status += IDASetUserData(ida_mem, (void *) &solver);
  //status += IDASetMaxNonlinIters(ida_mem, 60);
//AssertThrow(status == 0, ExcMessage("Error in IDA Solver"));
//    status += IDASetNonlinConvCoef(ida_mem, 10.0);
  //status += IDASetMaxOrd(ida_mem, 2);

  reset_ode(solution, solution_dot, initial_time, initial_step_size, max_steps);

  double next_time = 0;
  while ((t<final_time) && (step_number < max_steps))
    {

      next_time += outputs_period;
      cout << t <<"---->"<<next_time<<endl;
      status = IDASolve(ida_mem, next_time, &t, yy, yp, IDA_NORMAL);

      status = IDAGetLastStep(ida_mem, &h);
      AssertThrow(status == 0, ExcMessage("Error in IDA Solver"));
      cout << "Step " << step_number
           << ", t = " << t
           << ", h = " << h << endl;

      // Check the solution
      bool reset = solver.solution_check(solution, solution_dot, t, step_number, h);


      solver.output_step(solution, solution_dot, t, step_number, h);

      if ( reset == true )
        {
          NV_LENGTH_S(yy) = solution.size();
          NV_DATA_S(yy) = solution.begin();
          NV_LENGTH_S(yp) = solution_dot.size();
          NV_DATA_S(yp) = solution_dot.begin();

          double frac = 0;
          int k = 0;
          IDAGetLastOrder(ida_mem, &k);
          frac = std::pow((double)k,2.);
          reset_ode(solution, solution_dot, t,
                    h/frac, max_steps);
        }


      step_number++;
    }

  // Free the vectors which are no longer used.
  N_VDestroy_Serial(yy);
  N_VDestroy_Serial(yp);
  N_VDestroy_Serial(abs_tolls);
  N_VDestroy_Serial(diff_id);

  return step_number;
}


void DAETimeIntegrator::reset_ode(Vector<double> &solution,
                                  Vector<double> &solution_dot,
                                  double current_time,
                                  double current_time_step,
                                  unsigned int max_steps)
{
  if (ida_mem)
    IDAFree(&ida_mem);

  ida_mem = IDACreate();

  int status;
  Assert(solution.size() == solver.n_dofs(),
         ExcDimensionMismatch(solution.size(), solver.n_dofs()));

  Assert(solution_dot.size() == solver.n_dofs(),
         ExcDimensionMismatch(solution_dot.size(), solver.n_dofs()));

  // Free the vectors which are no longer used.
  if (yy)
    N_VDestroy_Serial(yy);
  if (yp)
    N_VDestroy_Serial(yp);
  if (abs_tolls)
    N_VDestroy_Serial(abs_tolls);
  if (diff_id)
    N_VDestroy_Serial(diff_id);


  // The solution is stored in
  // solution. Here we take only a
  // view of it.
  yy = N_VNewEmpty_Serial(solver.n_dofs());
  NV_DATA_S(yy) = solution.begin();

  //N_VPrint_Serial(yy);
  //solution_dot.print();
  yp = N_VNewEmpty_Serial(solver.n_dofs());
  NV_DATA_S(yp) = solution_dot.begin();
  //N_VPrint_Serial(yp);

  diff_id = N_VNewEmpty_Serial(solver.n_dofs());
  NV_DATA_S(diff_id) = solver.differential_components().begin();

  Vector<double> tolerances = solver.get_diameters();
  tolerances*=(1/tolerances.linfty_norm()*abs_tol);
  abs_tolls = N_VNewEmpty_Serial(solver.n_dofs());
  NV_DATA_S(abs_tolls) = tolerances.begin();
  //N_VPrint_Serial(abs_tolls);

  status = IDAInit(ida_mem, dae_residual, current_time, yy, yp);
  //status += IDASStolerances(ida_mem, rel_tol, abs_tol);
  status += IDASVtolerances(ida_mem, rel_tol, abs_tolls);
  status += IDASetInitStep(ida_mem, current_time_step);
  status += IDASetUserData(ida_mem, (void *) &solver);

  status += IDASetId(ida_mem, diff_id);
  status += IDASetSuppressAlg(ida_mem, true);

  status += IDASetMaxNumSteps(ida_mem, max_steps);
  status += IDASetStopTime(ida_mem, final_time);

  status += IDASetMaxNonlinIters(ida_mem, 10);


  if (use_iterative == true)
    {
      unsigned int maxl = 16;
      SUNLinearSolver LS = SUNSPGMR(yy, PREC_LEFT, maxl);

      /* IDA recommends allowing up to 5 restarts (default is 0) */
      status += SUNSPGMRSetMaxRestarts(LS, 5);

      /* Attach the linear sovler */
      status += IDASpilsSetLinearSolver(ida_mem, LS);


      
      
      //status += IDASpgmr(ida_mem, solver.n_dofs());
      if (provide_jac)
        //status += IDASpilsJacTimesVecFn(ida_mem, dae_jtimes);
        status += IDASpilsSetJacTimes(ida_mem, dae_setup_jac, dae_jtimes);
      if (provide_jac_prec)
        status += IDASpilsSetPreconditioner(ida_mem, dae_setup_prec, dae_prec);
    }
  else
    {
      //status += IDALapackDense(ida_mem, solver.n_dofs());
      AssertThrow(false , ExcMessage("Not Currently Implemented."));
    }

  status += IDASetMaxOrd(ida_mem, 5);
  //std::cout<<"???1"<<std::endl;

  AssertThrow(status == 0, ExcMessage("Error initializing IDA."));
  //std::cout<<"???1"<<std::endl;
  if (ic_type == 2)
    {
      // (re)initialization of the vectors
      //solution_dot = 0;
      if (current_time !=0)
        IDACalcIC(ida_mem, IDA_Y_INIT, current_time+current_time_step);
      IDAGetConsistentIC(ida_mem, yy, yp);
    }
  else if (ic_type == 1)
    {
      IDACalcIC(ida_mem, IDA_YA_YDP_INIT, current_time+current_time_step);
      IDAGetConsistentIC(ida_mem, yy, yp);
    }
  else if (ic_type == 3)
    {
      IDAGetConsistentIC(ida_mem, yy, yp);
      std::cout << "Using consistent conditions type 3" << std::endl;
    }
  Vector<double> resid(solver.n_dofs());
  solver.residual(current_time,resid,solution,solution_dot);

  //solution_dot.reinit(solver.n_dofs());
  //Vector<double> res(solver.n_dofs());
  //solver.residual(0,res,solution_dot,solution);
  //solution_dot -= res;
  //solver.output_step(solution, solution_dot, 0, 0, current_time_step);
}
