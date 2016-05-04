#ifndef __thiwi__ode_argument_h
#define __thiwi__ode_argument_h

#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_view.h>

using namespace dealii;

/** Base class that needs to be inherited by any function that wants
 * to use the time integrator class. */

class OdeArgument {
public :
    /** Returns the number of degrees of freedom. Pure virtual function. */
    virtual unsigned int n_dofs() const = 0;
       
    /** This function is called at the end of each iteration step for
     * the ode solver. Once again, the conversion between pointers and
     * other forms of vectors need to be done inside the inheriting
     * class. */
    virtual void output_step(Vector<double> & solution,
	                     Vector<double> &solution_dot,
                             const double t, 
			     const unsigned int step_number,
			     const double h) = 0;

    /** This function will check the behaviour of the solution. If it
     * is converged or if it is becoming unstable the time integrator
     * will be stopped. If the convergence is not achived the
     * calculation will be continued. If necessary, it can also reset
     * the time stepper. */
    virtual bool solution_check( Vector<double> &solution,
				 Vector<double> &solution_dot,
				 const double t,
				 const unsigned int step_number,
				 const double h) = 0;

				     /** For dae problems, we need a
					 residual function. */
    virtual int residual(const double t, 
			 Vector<double> &dst,  
			 const Vector<double> &src_yy,
			 const Vector<double> &src_yp) = 0;

				     /** Jacobian vector product. */
    virtual int jacobian(const double t,
			 Vector<double> &dst,  
			 const Vector<double> &src_yy,
			 const Vector<double> &src_yp,
			 const Vector<double> &src,
			 const double alpha);
    
				     /** Setup Jacobian preconditioner. */
    virtual int setup_jacobian_prec(const double t,
				    const Vector<double> &src_yy,
				    const Vector<double> &src_yp,
				    const double alpha);    
    
				     /** Jacobian preconditioner
					 vector product. */
    virtual int jacobian_prec(const double t,
			      Vector<double> &dst,  
			      const Vector<double> &src_yy,
			      const Vector<double> &src_yp,
			      const Vector<double> &src,
			      const double alpha);    
    
				     /** And an identification of the
					 differential components. This
					 has to be 1 if the
					 corresponding variable is a
					 differential component, zero
					 otherwise.  */
    virtual Vector<double> & differential_components();

    virtual Vector<double> & get_diameters();

    virtual ~OdeArgument() {};

    bool reset_time_integrator;

    bool stop_time_integrator;

};

#endif
