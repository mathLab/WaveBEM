//----------------------------  step-34.cc  ---------------------------
//    $Id: step-34.cc 18734 2009-04-25 13:36:48Z heltai $
//    Version: $Name$
//
//    Copyright (C) 2009, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//    Authors: Luca Heltai, Cataldo Manigrasso
//
//----------------------------  step-34.cc  ---------------------------

#ifndef bem_problem_h
#define bem_problem_h
// @sect3{Include files}

// The program starts with including a bunch
// of include files that we will use in the
// various parts of the program.


#include <deal.II/base/smartpointer.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/quadrature_selector.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_direct.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>


#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1_eulerian.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/solution_transfer.h>

// And here are a few C++ standard header
// files that we will need:
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <map>

#include "../include/octree_block.h"
#include "../include/local_expansion.h"
#include "../include/multipole_expansion.h"
#include "../include/ass_leg_function.h"
#include "../include/computational_domain.h"
#include "../include/bem_fma.h"
#include "../include/constrained_matrix.h"



using namespace dealii;


template <int dim>
class BEMProblem
{
public:

  typedef typename DoFHandler<dim-1,dim>::active_cell_iterator cell_it;

  BEMProblem(ComputationalDomain<dim> &comp_dom,BEMFMA<dim> &fma);

  void solve(Vector<double> &phi, Vector<double> &dphi_dn,
             const Vector<double> &tmp_rhs);

  void reinit();

  void compute_constraints(AffineConstraints<double> &constraints, const Vector<double> &tmp_rhs);

  //  private:

  void declare_parameters(ParameterHandler &prm);

  void parse_parameters(ParameterHandler &prm);

  // To be commented

  void compute_alpha();

  void assemble_system();

  // The next three methods are
  // needed by the GMRES solver:
  // the first provides result of
  // the product of the system
  // matrix (a combination of Neumann
  // and Dirichlet matrices) by the
  // vector src. The result is stored
  // in the vector dst.

  void vmult(Vector<double> &dst, const Vector<double> &src) const;

  // The second method computes the
  // right hand side vector of the
  // system.

  void compute_rhs(Vector<double> &dst, const Vector<double> &src) const;

  // The third method computes the
  // product between the solution vector
  // and the (fully populated) sytstem
  // matrix.

  void assemble_preconditioner();

  void compute_surface_gradients(const Vector<double> &tmp_rhs);

  void solve_system(Vector<double> &phi, Vector<double> &dphi_dn,
                    const Vector<double> &tmp_rhs);

  void residual(Vector<double> &res, const Vector<double> &phi,
                const Vector<double> &dphi_dn);

  void output_results(const std::string);

  ComputationalDomain<dim> &comp_dom;

  BEMFMA<dim> &fma;

  FullMatrix<double>    neumann_matrix;
  FullMatrix<double>    dirichlet_matrix;
  Vector<double>        system_rhs;

  Vector<double>              sol;
  Vector<double>              alpha;

  Vector<double>              serv_phi;
  Vector<double>              serv_dphi_dn;
  Vector<double>              serv_tmp_rhs;

  AffineConstraints<double>     constraints;

  std::string solution_method;

  SolverControl solver_control;

  SparseDirectUMFPACK preconditioner;

  SparsityPattern preconditioner_sparsity_pattern;

  int preconditioner_band;

  bool is_preconditioner_initialized;

  std::vector<Point<dim> > node_surface_gradients;

};

#endif
