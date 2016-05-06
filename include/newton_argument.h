#ifndef newton_argument_h
#define newton_argument_h

#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_view.h>

using namespace dealii;

/** Base class that needs to be inherited by any function that wants
 * to use the newton solver class. */

class NewtonArgument
{
public :
  /** Returns the number of degrees of freedom. Pure virtual function. */
  virtual unsigned int n_dofs() const = 0;

  /** This function is called at the end of each iteration step for
   * the newton solver. Once again, the conversion between pointers and
   * other forms of vectors need to be done inside the inheriting
   * class. */
  virtual void output_step(Vector<double> &solution,
                           const unsigned int step_number) = 0;

  /** For dae problems, we need a
  residual function. */
  virtual int residual(Vector<double> &dst,
                       const Vector<double> &src_yy) = 0;

  /** Jacobian vector product. */
  virtual int jacobian(Vector<double> &dst,
                       const Vector<double> &src_yy,
                       const Vector<double> &src);

  /** Setup Jacobian preconditioner (builds the
                              Jacobian preconditioner matrix). */
  virtual int setup_jacobian_prec(const Vector<double> &src_yy);

  /** Inverse Jacobian preconditioner
  vector product. */
  virtual int jacobian_prec(Vector<double> &dst,
                            const Vector<double> &src_yy,
                            const Vector<double> &src);

  /** Jacobian preconditioner
  vector product. */
  virtual int jacobian_prec_prod(Vector<double> &dst,
                                 const Vector<double> &src_yy,
                                 const Vector<double> &src);

  virtual ~NewtonArgument() {};

};

#endif
