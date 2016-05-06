#ifndef dae_time_integrator_h
#define dae_time_integrator_h

#include <deal.II/base/logstream.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/parameter_handler.h>

// For time integration.
#include <ida/ida.h>
#include <ida/ida_dense.h>
#include <ida/ida_lapack.h>
#include <ida/ida_spgmr.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>


#include "ode_argument.h"


class DAETimeIntegrator
{
public:
  /** Constructor for the DAETimeIntegrator class. The Solver class is
   * required to have a Solver.solve(Vector<double> &dst, const
   * Vector<double> &src) method that will be called by the time
   * integrator to find out about the solution to a given src. */
  DAETimeIntegrator(OdeArgument &solver);

  /** In cases in which we don't get  the integrator parameters
   * from a file, we use a constructor in which they are directly
   * specified. */
  DAETimeIntegrator(OdeArgument &solver,
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
                    bool provide_jac_prec);

  /** House cleaning. */
  ~DAETimeIntegrator();

  /** Declare parameters for this class to function properly. */
  static void declare_parameters(ParameterHandler &prm);

  /** Parse a parameter handler. */
  void parse_parameters(ParameterHandler &prm);

  /** Evolve. This function returns the final number of steps. */
  unsigned int start_ode(Vector<double> &solution,
                         Vector<double> &solution_dot,
                         const unsigned int max_steps);

  /** Clear internal memory, and
  start with clean
  objects. This is useful if
  you need to refine your
  mesh between stesp. */
  void reset_ode(Vector<double> &y, Vector<double> &yp,
                 double t, double h, unsigned int max_steps);


  /** Final time. */
  double final_time;

  /** Initial time for the ode.*/
  double initial_time;

private:
  /** The bubble membrane poperties. */
  OdeArgument &solver;

  /** Initial step size. */
  double initial_step_size;

  /** Minimum step size. */
  double min_step_size;

  /** Absolute error tolerance for adaptive time stepping. */
  double abs_tol;

  /** Relative error tolerance for adaptive time stepping. */
  double rel_tol;

  /** Maximum number of time steps. */
  unsigned int max_n_steps;

  double outputs_period;

  unsigned int ic_type;

  /** Initialization flag.*/
  bool is_initialized;

  /** If true, we use
  preconditioned gmres. */
  bool use_iterative;

  /** Provides the Jacobian vector
  product. If this is false,
  then finite difference is
  used. */
  bool provide_jac;
  /** Provide preconditioning for
  Jacobian. If not set, then no
  preconditioning is used. */
  bool provide_jac_prec;


  /** Ida memory object. */
  void *ida_mem;

  /** Ida solution vector. */
  N_Vector yy;
  /** Ida solution derivative vector. */
  N_Vector yp;
  /** Ida absolute tolerances vector. */
  N_Vector abs_tolls;
  /** Ida differential components vector. */
  N_Vector diff_id;


};

#endif
