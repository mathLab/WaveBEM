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


				 // @sect3{Include files}

				 // The program starts with including a bunch
				 // of include files that we will use in the
				 // various parts of the program. Most of them
				 // have been discussed in previous tutorials
				 // already:
#ifndef free_surface_h
#define free_surface_h

			 
#include <base/smartpointer.h>
#include <base/convergence_table.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/quadrature_selector.h>
#include <base/parsed_function.h>
#include <base/utilities.h>

#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
//#include <lac/matrix_lib.h>
#include <lac/solver_control.h>
//#include <lac/solver_gmres.h>
//#include <lac/precondition.h>
#include <lac/constraint_matrix.h>


#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/grid_generator.h>
#include <grid/grid_in.h>
#include <grid/grid_out.h>
#include <grid/tria_boundary_lib.h>
#include <grid/grid_refinement.h>

#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <dofs/dof_renumbering.h>

#include <fe/fe_q.h>
#include <fe/fe_values.h>
#include <fe/fe_system.h>
#include <fe/mapping_q1_eulerian.h>
#include <fe/mapping_q1.h>

#include <numerics/data_out.h>
#include <numerics/vector_tools.h>
#include <numerics/solution_transfer.h>
#include <numerics/fe_field_function.h>

				 // And here are a few C++ standard header
				 // files that we will need:
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <set>


#include <dae_time_integrator.h>
//#include <newton_solver.h>
#include <newton_argument.h>

#include "../include/bem_problem.h"
#include "../include/numerical_towing_tank.h"


template <int dim>
class FreeSurface:
public OdeArgument, 
public NewtonArgument 
{
public:
  FreeSurface(NumericalTowingTank &comp_dom, BEMProblem<dim> &bem) :
		  wind(dim),  comp_dom(comp_dom), bem(bem)
      {dofs_number = 0,
       output_frequency = 1;
       stop_time_integrator = false;
       reset_time_integrator = false;}
    
  virtual unsigned int n_dofs() const;

  virtual void output_step(Vector<double> & solution,
                           const unsigned int step_number);
  
  virtual void output_step(Vector<double> & solution,
	                   Vector<double> &solution_dot,
			   const double t,
	                   const unsigned int step_number,
		           const double  h);

  virtual bool solution_check(Vector<double> &solution,
			      Vector<double> &solution_dot,
			      const double t,
			      const unsigned int step_number,
			      const double h);


				     /** For dae problems, we need a
					 residual function. This one is computed with
                                         sacado to allow computations of derivatives for jacobian*/
  virtual int residual(const double t, 
     		       Vector<double> &dst,  
		       const Vector<double> &src_yy,
		       const Vector<double> &src_yp);

				     /** For newton solver, we need a
					 residual function. This one is computed with
                                         sacado to allow computations of derivatives for jacobian*/
  virtual int residual(Vector<double> &dst,  
	               const Vector<double> &src_yy);

				     /** Jacobian vector product for dae. */
  virtual int jacobian(const double t,
		       Vector<double> &dst,  
		       const Vector<double> &src_yy,
		       const Vector<double> &src_yp,
		       const Vector<double> &src,
		       const double alpha);


				     /** Jacobian vector product for newton solver. */
  virtual int jacobian(Vector<double> &dst,  
	               const Vector<double> &src_yy,
		       const Vector<double> &src);

				     /** This function computes either DAE residual
                                         or corrensponding Jacobian matrix vector product
                                         with vector src. Boolean is_jacobian determines if
                                         this function is used for residual (false) or for
                                         Jacobian (true). In residual case, src and alpha
                                         values assigned are disregarded.*/
  int residual_and_jacobian(const double t,
		            Vector<double> &dst,  
		            const Vector<double> &src_yy,
		            const Vector<double> &src_yp,
		            const Vector<double> &src,
		            const double alpha,
                            const bool is_jacobian);
    

				     /** Setup Jacobian preconditioner for Newton. */
  virtual int setup_jacobian_prec(const Vector<double> &src_yy);    

				     /** Setup Jacobian preconditioner for DAE. */
  virtual int setup_jacobian_prec(const double t,
				  const Vector<double> &src_yy,
				  const Vector<double> &src_yp,
				  const double alpha);    

  void compute_constraints(ConstraintMatrix &constraints, 
			   ConstraintMatrix &vector_constraints);
    
				     /** Jacobian inverse preconditioner
					 vector product for dae. */
  virtual int jacobian_prec(const double t,
			    Vector<double> &dst,  
			    const Vector<double> &src_yy,
			    const Vector<double> &src_yp,
			    const Vector<double> &src,
			    const double alpha);    
    
				     /** Jacobian inverse preconditioner
					 vector product for newton solver. */
  virtual int jacobian_prec(Vector<double> &dst,  
			    const Vector<double> &src_yy,
			    const Vector<double> &src);

				     /** Jacobian preconditioner
					 vector product for newton solver. */
  virtual int jacobian_prec_prod(Vector<double> &dst,  
			    const Vector<double> &src_yy,
			    const Vector<double> &src);


				     /** And an identification of the
					 differential components. This
					 has to be 1 if the
					 corresponding variable is a
					 differential component, zero
					 otherwise.  */
  virtual Vector<double> & differential_components();

				     /** Method to enforce mesh conformity
                                         at edges at each restart */
  void make_edges_conformal(Vector<double> & solution,
		            Vector<double> &solution_dot,
			    const double t,
	                    const unsigned int step_number,
		            const double  h);

                                      // in the first layer of water cells past
                                      // the transom there can't be hanging nodes:
                                      // this method removes them
  void remove_transom_hanging_nodes(Vector<double> & solution,
                                    Vector<double> &solution_dot,
                                    const double t,
                                    const unsigned int step_number,
                                    const double  h);

				     /** Method to make sure residual is
					 null at each (re)start of the
					 computation */
  void prepare_restart(const double t, Vector<double> &y, Vector<double> &yp, bool restart_flag=true);

  typedef typename DoFHandler<dim-1,dim>::active_cell_iterator cell_it;
  typedef typename Triangulation<dim-1,dim>::active_cell_iterator tria_it;  

  void declare_parameters(ParameterHandler &prm);
  
  void parse_parameters(ParameterHandler &prm);
  
  void initial_conditions(Vector<double> &dst);

  void reinit();

  void prepare_bem_vectors(double time,
			   Vector<double> &bem_bc,
			   Vector<double> &dphi_dn) const;

                                     // In this routine we use the
				     // solution obtained to compute the 
				     // DXDt and DphiDt on the domain
				     // boundary. 


  void compute_DXDt_and_DphiDt(double time,
                               const Vector<double> & phi,
                               const Vector<double> & dphi_dn,
                               const Vector<double> & nodes_velocities);
    
                                     // In this routine we compute the
                                     // complete potential gradient,
                                     // surface potential gradient,
                                     // surface normal

  void vmult(Vector<double> &dst, const Vector<double> &src) const; 

  void compute_potential_gradients(Vector<double> &complete_potential_gradients,
                                   const Vector<double> & phi,
                                   const Vector<double> & dphi_dn);


  void compute_keel_smoothing(Vector<double> & smoothing);

    
  void compute_pressure(Vector<double> & pressure, 
                        Vector<double> & comp_1, Vector<double> & comp_2, Vector<double> & comp_3, Vector<double> & comp_4, 
			  const double t,
			  const Vector<double> & solution,
			  const Vector<double> & solution_dot);

  void output_results(const std::string,
                      const double t,
                      const Vector<double> & solution,
                      const Vector<double> & pressure);

  Vector<double>& get_diameters();

  inline unsigned int Rhs_evaluations_counter()
        {return rhs_evaluations_counter;}

    void dump_solution(const Vector<double> &y,
		       const Vector<double> &yp,
		       const std::string fname) const;

    void restore_solution(Vector<double> &y,
			  Vector<double> &yp,
			  const std::string fname);
    

    void enforce_partial_geometry_constraints(const double blend_factor);

    void enforce_full_geometry_constraints(); 
    
    bool initial_condition_from_dump;
    std::string restore_filename;
  Vector<double> sys_comp;
private:  

  std::string output_file_name;
  double remeshing_period;
  double dumping_period;
    
  Functions::ParsedFunction<dim> wind;
  Functions::ParsedFunction<dim> initial_wave_shape;
  Functions::ParsedFunction<dim> initial_wave_potential;
  Functions::ParsedFunction<dim> initial_norm_potential_grad;
  Functions::ParsedFunction<dim> inflow_norm_potential_grad;

  
  std::string node_displacement_type;
  
  SolverControl solver_control;

  NumericalTowingTank &comp_dom;

  BEMProblem<dim> &bem;
  
  unsigned int dofs_number;
  
  unsigned int output_frequency;

  unsigned int refinement_level_on_boat;

                              // vectors for the DAE solver  
  Vector<double> DXDt_and_DphiDt_vector;
  Vector<double> diff_comp;
  Vector<double> alg_comp;
  Vector<double> temp_src;
  Vector<double> bem_residual;
  Vector<double> bem_phi;
  Vector<double> bem_dphi_dn;


  
  double initial_time;

    

                              // set of vectors and (mass) matrices
                              // needed in compute_DXDt_and_DphiDt 
  ConstraintMatrix     constraints;
  ConstraintMatrix     vector_constraints;

  SparsityPattern      DphiDt_sparsity_pattern;
  SparseMatrix<double> DphiDt_sys_matrix;
  SparseMatrix<double> DphiDt_sys_matrix_2;

  Vector<double>       DphiDt_sys_solution;
  Vector<double>       DphiDt_sys_solution_2;
  Vector<double>       DphiDt_sys_solution_3;

  Vector<double>       DphiDt_sys_rhs;
  Vector<double>       DphiDt_sys_rhs_2;
  Vector<double>       DphiDt_sys_rhs_3;
  Vector<double>       DphiDt_sys_rhs_4;
   
  SparsityPattern      vector_sparsity_pattern;
  SparseMatrix<double> vector_sys_matrix;
  SparseMatrix<double> vector_beltrami_matrix;
 
  Vector<double>       vector_sys_solution;
  Vector<double>       vector_sys_solution_2;

  Vector<double>       vector_sys_rhs;
  Vector<double>       vector_sys_rhs_2;


  SparsityPattern      jacobian_sparsity_pattern;
  SparseMatrix<double> jacobian_matrix;
  SparseMatrix<double> jacobian_dot_matrix;
  SparseMatrix<double> jacobian_preconditioner_matrix;
  SparseMatrix<double> preconditioner_preconditioner_matrix;
                              // vector storing the diameters of cells around each node
                              // needed to estimate local tolerances for 
                              // the ode solver

  Vector<double>       diameters;

                              // max y coodiante value (needed for
  		              // semi-lagrangian free surface
			      // deformation)
     
  double max_y_coor_value;

  unsigned int rhs_evaluations_counter;

    double min_diameter;

    unsigned int max_number_of_dofs;

    double coarsening_fraction;

    double refinement_fraction;

    bool sync_bem_with_geometry;

    bool restart_flag;

    double last_remesh_time;
                              // working copy of the map points vector 
  Vector<double>  working_map_points;
                              // working copy of the map points vector 
  Vector<double>  working_nodes_velocities;
                              // vector for position residuals
  Vector<double>  nodes_pos_res;
                              // vector for geometric treatment residuals
                              // (to be passed to num_ tow_tank)
  Vector<double>  nodes_ref_surf_dist;

                              // vector for jacobian X delta vector
                              // on differential nodes
  Vector<double>  nodes_diff_jac_x_delta;
                              // vector for jacobian X delta vector
                              // on algebraic nodes
  Vector<double>  nodes_alg_jac_x_delta;
                              // residual of dae equation
  Vector<double>  dae_nonlin_residual;
                              // residual of dae linear step
  Vector<double>  dae_linear_step_residual;

  double alpha;

  Vector<double> current_sol;
  
  Vector<double> current_sol_dot;

  double current_time;

  double ref_transom_wet_surface;

  Vector <double> ref_cell_areas;

  Vector <double> break_wave_press;

  Vector <double> transom_pressure_patch; // temporary

                              // this number sets the bottom limit of the
                              // adaptive refinement. when set to 2.0, it
                              // limits the refinements so that all the cells
                              // have a higher diameter w/r to the min_diameter
                              // cell of the initial mesh. if set to 1.0, it allows
                              // cells to be refined to half of the initial min_diameter,
                              // and so on. On the other side, a value of 4.0 means
                              // that cells being 2 times bigger than initial min_diameter
                              // won't be refined.
  double adaptive_ref_limit;

  //to be removed soooooon!!
  Vector<double> cell_flag_vector_before;


};

#endif

