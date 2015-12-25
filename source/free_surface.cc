//----------------------------  step-34.cc  ---------------------------
//    $Id: step-34.cc 18734 2009-04-25 13:36:48Z heltai $
//    Version: $Name$ 
//
//    Copyright (C) 2009, 2010, 2011 by the deal.II authors 
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//    Authors: Luca Heltai, Andrea Mola
//
//----------------------------  step-34.cc  ---------------------------


				 // @sect3{Include files}

				 // The program starts with including a bunch
				 // of include files that we will use in the
				 // various parts of the program. Most of them
				 // have been discussed in previous tutorials
				 // already:



//#include <deal.II/lac/vector.h>		 
//#include <deal.II/lac/sparse_direct.h>
//#include <deal.II/lac/constraint_matrix.h>
//#include <deal.II/lac/solver_cg.h>
//#include <deal.II/lac/vector_view.h>
//#include <deal.II/grid/grid_refinement.h>
//#include <deal.II/numerics/matrices.h>

#include <Sacado.hpp>
#include "../include/free_surface.h"
#include "../include/restart_nonlinear_problem_diff.h"
#include "../include/restart_nonlinear_problem_alg.h"
#include "newton_solver.h"

#include <deal.II/numerics/fe_field_function.h>

#include <GeomPlate_BuildPlateSurface.hxx>
#include <GeomPlate_PointConstraint.hxx>
#include <GeomPlate_MakeApprox.hxx>
#include <Geom_Surface.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <GeomPlate_Surface.hxx>
#include <BRepAdaptor_Curve.hxx>

typedef Sacado::Fad::DFad<double> fad_double;

using namespace dealii;
using namespace OpenCascade;


template <int dim>
void FreeSurface<dim>::declare_parameters(ParameterHandler &prm) {

  prm.declare_entry("Output file name", "waveResult", Patterns::Anything());
  prm.declare_entry("Remeshing period", "1", Patterns::Double());
  prm.declare_entry("Dumping period", "1", Patterns::Double());
  prm.declare_entry("Restore from file", "", Patterns::Anything());
  prm.declare_entry("Keep BEM in sync with geometry", "true", Patterns::Bool());

  prm.declare_entry("Is x-translation imposed", "true", Patterns::Bool());
  prm.enter_subsection("Hull x-axis translation");
  {
    Functions::ParsedFunction<1>::declare_parameters(prm, 1);
    prm.set("Function expression", "0");
  }
  prm.leave_subsection();

  prm.declare_entry("Is y-translation imposed", "true", Patterns::Bool());
  prm.enter_subsection("Hull y-axis translation");
  {
    Functions::ParsedFunction<1>::declare_parameters(prm, 1);
    prm.set("Function expression", "0");
  }
  prm.leave_subsection();

  prm.declare_entry("Is z-translation imposed", "true", Patterns::Bool());
  prm.enter_subsection("Hull z-axis translation");
  {
    Functions::ParsedFunction<1>::declare_parameters(prm, 1);
    prm.set("Function expression", "0");
  }
  prm.leave_subsection();

  prm.enter_subsection("Wind function 2d");
  {
    Functions::ParsedFunction<2>::declare_parameters(prm, 2);
    prm.set("Function expression", "1; 1");
  }
  prm.leave_subsection();

  prm.enter_subsection("Wind function 3d");
  {
    Functions::ParsedFunction<3>::declare_parameters(prm, 3);
    prm.set("Function expression", "1; 1; 1");
  }
  prm.leave_subsection();

  prm.enter_subsection("Initial Wave Shape 2d");
  {
    Functions::ParsedFunction<2>::declare_parameters(prm);
    prm.set("Function expression", "0");
  }
  prm.leave_subsection();

  prm.enter_subsection("Initial Wave Shape 3d");
  {
    Functions::ParsedFunction<3>::declare_parameters(prm);
    prm.set("Function expression", "0");
  }
  prm.leave_subsection();
  
  prm.enter_subsection("Initial Wave Potential 2d");
  {
    Functions::ParsedFunction<2>::declare_parameters(prm);
    prm.set("Function expression", "0");
  }
  prm.leave_subsection();

  prm.enter_subsection("Initial Wave Potential 3d");
  {
    Functions::ParsedFunction<3>::declare_parameters(prm);
    prm.set("Function expression", "0");
  }
  prm.leave_subsection();

  prm.enter_subsection("Initial Wave Normal Potential Gradient 2d");
  {
    Functions::ParsedFunction<2>::declare_parameters(prm);
    prm.set("Function expression", "0");
  }
  prm.leave_subsection();

  prm.enter_subsection("Initial Wave Normal Potential Gradient 3d");
  {
    Functions::ParsedFunction<3>::declare_parameters(prm);
    prm.set("Function expression", "0");
  }
  prm.leave_subsection();

  prm.enter_subsection("Inflow Normal Potential Gradient 2d");
  {
    Functions::ParsedFunction<2>::declare_parameters(prm);
    prm.set("Function expression", "0");
  }
  prm.leave_subsection();

  prm.enter_subsection("Inflow Normal Potential Gradient 3d");
  {
    Functions::ParsedFunction<3>::declare_parameters(prm);
    prm.set("Function expression", "0");
  }
  prm.leave_subsection(); 

  prm.declare_entry("Node displacement type", "lagrangian", 
                    Patterns::Selection("lagrangian|semilagrangian"));
		    

  prm.enter_subsection("Local Refinement");
  
  prm.declare_entry("Maximum number of dofs", "4500", 
                    Patterns::Integer());

  prm.declare_entry("Coarsening fraction", ".1", 
                    Patterns::Double());

  prm.declare_entry("Refinement fraction", ".3", 
                    Patterns::Double());

  prm.declare_entry("Adaptive refinement limit", "2.0", 
                    Patterns::Double());

  prm.leave_subsection();
  
  
  prm.enter_subsection("Solver");
  SolverControl::declare_parameters(prm);
  prm.leave_subsection();

}

template <int dim>
void FreeSurface<dim>::parse_parameters(ParameterHandler &prm) {

  output_file_name = prm.get("Output file name");
  remeshing_period = prm.get_double("Remeshing period");
  dumping_period = prm.get_double("Dumping period");

  restore_filename = prm.get("Restore from file");
  initial_condition_from_dump = (restore_filename != "");
  
  sync_bem_with_geometry = prm.get_bool("Keep BEM in sync with geometry");
  
  prm.enter_subsection(std::string("Wind function ")+
		       Utilities::int_to_string(dim)+std::string("d"));
  {
    wind.parse_parameters(prm);
  }
  prm.leave_subsection();

  is_hull_x_translation_imposed = prm.get_bool("Is x-translation imposed");
  prm.enter_subsection("Hull x-axis translation");
  {
    hull_x_axis_translation.parse_parameters(prm);
  }
  prm.leave_subsection();

  is_hull_y_translation_imposed = prm.get_bool("Is y-translation imposed");
  prm.enter_subsection("Hull y-axis translation");
  {
    hull_y_axis_translation.parse_parameters(prm);
  }
  prm.leave_subsection();

  is_hull_z_translation_imposed = prm.get_bool("Is z-translation imposed");
  prm.enter_subsection("Hull z-axis translation");
  {
    hull_z_axis_translation.parse_parameters(prm);
  }
  prm.leave_subsection();

  prm.enter_subsection(std::string("Initial Wave Shape ")+
		       Utilities::int_to_string(dim)+std::string("d"));
  {
    initial_wave_shape.parse_parameters(prm);
  }
  prm.leave_subsection();
 
  prm.enter_subsection(std::string("Initial Wave Potential ")+
		       Utilities::int_to_string(dim)+std::string("d"));
  {
    initial_wave_potential.parse_parameters(prm);
  }
  prm.leave_subsection();

  prm.enter_subsection(std::string("Initial Wave Normal Potential Gradient ")+
		       Utilities::int_to_string(dim)+std::string("d"));
  {
    initial_norm_potential_grad.parse_parameters(prm);
  }
  prm.leave_subsection();
 
  prm.enter_subsection(std::string("Inflow Normal Potential Gradient ")+
		       Utilities::int_to_string(dim)+std::string("d"));
  {
    inflow_norm_potential_grad.parse_parameters(prm);
  }
  prm.leave_subsection();

  node_displacement_type = prm.get("Node displacement type");
  
  prm.enter_subsection("Solver");
  solver_control.parse_parameters(prm);
  prm.leave_subsection();


  prm.enter_subsection("Local Refinement");
  
  max_number_of_dofs = prm.get_integer("Maximum number of dofs");
  coarsening_fraction = prm.get_double("Coarsening fraction");
  refinement_fraction = prm.get_double("Refinement fraction");
  adaptive_ref_limit = prm.get_double("Adaptive refinement limit");

  prm.leave_subsection();
  
}


template <int dim>
void FreeSurface<dim>::reinit() {
    
  dofs_number = comp_dom.vector_dh.n_dofs() + // nodes positions
                comp_dom.dh.n_dofs() +        // nodes potential
                comp_dom.dh.n_dofs() +        // nodes normal potential grandient
                6 +                           // these are the dofs needed for hull rigid velocities and displacements
                7;                            // these are the dofs needed for hull rigid angular velocities and quaternions
  
  DXDt_and_DphiDt_vector.reinit(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()); 


                                      // we initialize here all the
                                      // sparsity patterns

  compute_constraints(constraints, vector_constraints);

  DphiDt_sparsity_pattern.reinit (comp_dom.dh.n_dofs(),
                                   comp_dom.dh.n_dofs(),
                                   comp_dom.dh.max_couplings_between_dofs());
  vector_sparsity_pattern.reinit (comp_dom.vector_dh.n_dofs(),
                                     comp_dom.vector_dh.n_dofs(),
                                     comp_dom.vector_dh.max_couplings_between_dofs());


  DoFTools::make_sparsity_pattern (comp_dom.dh, DphiDt_sparsity_pattern, constraints);
  DphiDt_sparsity_pattern.compress();

  DoFTools::make_sparsity_pattern (comp_dom.vector_dh, vector_sparsity_pattern, vector_constraints);
  vector_sparsity_pattern.compress();


  working_map_points.reinit (comp_dom.vector_dh.n_dofs());
  working_nodes_velocities.reinit(comp_dom.vector_dh.n_dofs());
  nodes_pos_res.reinit(comp_dom.vector_dh.n_dofs());
  nodes_ref_surf_dist.reinit(comp_dom.vector_dh.n_dofs());
  nodes_diff_jac_x_delta.reinit(dofs_number);
  nodes_alg_jac_x_delta.reinit(dofs_number);
  dae_nonlin_residual.reinit(dofs_number);
  dae_linear_step_residual.reinit(dofs_number);
  current_sol.reinit(dofs_number);
  current_sol_dot.reinit(dofs_number);
  bem_residual.reinit(comp_dom.dh.n_dofs());
  bem_phi.reinit(comp_dom.dh.n_dofs());
  bem_dphi_dn.reinit(comp_dom.dh.n_dofs());
  temp_src.reinit(comp_dom.vector_dh.n_dofs());
  break_wave_press.reinit(comp_dom.dh.n_dofs());
                                      // we initialize the rhs evaluations counter
  rhs_evaluations_counter = 0;
} 


template <int dim>
void FreeSurface<dim>::initial_conditions(Vector<double> &dst) {

  // da trattare il caso di dumped solution (continuazione conto salvato)
  initial_time = 0.0;
  initial_wave_shape.set_time(initial_time);
  initial_wave_potential.set_time(initial_time);
  initial_norm_potential_grad.set_time(initial_time);
  inflow_norm_potential_grad.set_time(initial_time);
  wind.set_time(initial_time);

  dst.reinit(dofs_number);

  Point<1> time(initial_time);
  Point<1> delta_t(1e-5);

  if (is_hull_x_translation_imposed)
     {
     restart_hull_displacement(0) = hull_x_axis_translation.value(time);
     dst(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+3) = restart_hull_displacement(0);
     //double vel = ( hull_x_axis_translation.value(time+delta_t) -
     //               hull_x_axis_translation.value(time+(-1.0*delta_t)) ) / 2.0 / delta_t(0);
     double vel = 0.0;
     cout<<"VELX: "<<vel<<endl;
     dst(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()) = vel;
     }
  else
   restart_hull_displacement(0) = 0.0;

  if (is_hull_y_translation_imposed)
     {
     restart_hull_displacement(1) = hull_y_axis_translation.value(time);
     dst(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+4) = restart_hull_displacement(1);
     //double vel = ( hull_y_axis_translation.value(time+delta_t) -
     //               hull_y_axis_translation.value(time+(-1.0*delta_t)) ) / 2.0 / delta_t(0);
     double vel = 0.0;
     cout<<"VELY: "<<vel<<endl;
     dst(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+1) = vel;
     }
  else
   restart_hull_displacement(1) = 0.0;

  if (is_hull_z_translation_imposed)
     {
     restart_hull_displacement(2) = hull_z_axis_translation.value(time);
     dst(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+5) = restart_hull_displacement(2);
     double vel = ( hull_z_axis_translation.value(time+delta_t) -
                    hull_z_axis_translation.value(time+(-1.0*delta_t)) ) / 2.0 / delta_t(0);
     //double vel = 0.0; //0.01*2.0*dealii::numbers::PI*cos(2.0*dealii::numbers::PI*initial_time);
     cout<<"VELZ: "<<vel<<endl;
     dst(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+2) = vel;
     }
  else
   restart_hull_displacement(2) = 0.0;

  Point<3> init_omega(0.0,0.0,0.0);
  Point<3> quaternion_vect(0.0,0,0.0);
  double quaternion_scalar = 1.0;
  for (unsigned int d=0; d<3; ++d)
      {
      dst(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+6+d) = init_omega(d);
      dst(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+9+d) = quaternion_vect(d);
      }
  dst(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+12) = quaternion_scalar;

  gp_Trsf curr_transf = comp_dom.boat_model.set_current_position(restart_hull_displacement,
                                                                 quaternion_scalar,
                                                                 quaternion_vect);

  restart_hull_location = curr_transf;
  restart_transom_center_point = comp_dom.boat_model.CurrentPointCenterTransom; 
  restart_transom_left_point = comp_dom.boat_model.CurrentPointLeftTransom;
  restart_transom_right_point = comp_dom.boat_model.CurrentPointRightTransom;
  restart_transom_left_tangent = comp_dom.boat_model.current_left_transom_tangent;
  restart_transom_right_tangent = comp_dom.boat_model.current_right_transom_tangent;

  
  Vector<double> instantWindValue(dim);
  Point<dim> zero(0,0,0);
  wind.vector_value(zero,instantWindValue);
  Point<dim> Vinf;
  for (unsigned int i = 0; i < dim; i++)
      Vinf(i) = instantWindValue(i);
  std::cout<<std::endl<<"Initial conditions: simulation time= "<<initial_time<<"   Vinf= ";
  instantWindValue.print(std::cout,4,false,true);

  std::vector<Point<dim> > support_points(comp_dom.dh.n_dofs());
  DoFTools::map_dofs_to_support_points<dim-1, dim>( *comp_dom.mapping, comp_dom.dh, support_points);
  std::vector<Point<dim> > vector_support_points(comp_dom.vector_dh.n_dofs());
  DoFTools::map_dofs_to_support_points<dim-1, dim>( *comp_dom.mapping, comp_dom.vector_dh, vector_support_points);


  for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
      {
      if ((comp_dom.flags[i] & water) ||  
          (comp_dom.flags[i] & near_water) )
         comp_dom.map_points(3*i+2) =  initial_wave_shape.value(vector_support_points[3*i]);
      }

  unsigned int j = dim-1;  
  Vector<double> geom_res(comp_dom.vector_dh.n_dofs());
  if (!comp_dom.no_boat)
     comp_dom.evaluate_ref_surf_distances(geom_res,false);

    //enforce_full_geometry_constraints(); //*********************
    //comp_dom.update_support_points(); //*******************

  for(unsigned int i=0; i<comp_dom.vector_dh.n_dofs(); ++i)
      {
      comp_dom.map_points(i) -= geom_res(i);
      } 
   
  if (!comp_dom.no_boat && comp_dom.boat_model.is_transom)
     {
     comp_dom.update_support_points();
     
     double transom_draft = fabs(comp_dom.boat_model.CurrentPointCenterTransom(2));
     double transom_aspect_ratio = (fabs(comp_dom.boat_model.CurrentPointLeftTransom(1))+
                                     fabs(comp_dom.boat_model.CurrentPointRightTransom(1)))/transom_draft;

     wind.set_time(100000000.0);
     Vector<double> instantWindValueTinf(dim);
     Point<dim> zero(0,0,0);
     wind.vector_value(zero,instantWindValueTinf);
     Point<dim> VinfTinf;
     for (unsigned int i = 0; i < dim; i++)
       VinfTinf(i) = instantWindValueTinf(i);
     wind.set_time(initial_time);

     double FrT = sqrt(VinfTinf*VinfTinf)/sqrt(9.81*transom_draft);
     double ReT = sqrt(9.81*pow(transom_draft,3.0))/1.307e-6;
     double eta_dry = fmin(0.05*pow(FrT,2.834)*pow(transom_aspect_ratio,0.1352)*pow(ReT,0.01338),1.0);
     double lh = 0.0;
     //if (eta_dry < 1.0) 
        lh = 5.0; 
        //lh = 0.1135*pow(FrT,3.025)*pow(transom_aspect_ratio,0.4603)*pow(ReT,-0.1514);
        //lh = 0.3265*pow(FrT,3.0) - 1.7216*pow(FrT,2.0) + 2.7593*FrT;
     cout<<"FrT: "<<FrT<<"  Tar: "<<transom_aspect_ratio<<"  ReT: "<<ReT<<endl;
     cout<<"****eta_dry: "<<eta_dry<<endl;


     for(unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
        {
        if ( (comp_dom.support_points[i](1) < comp_dom.boat_model.CurrentPointRightTransom(1)) &&
             (comp_dom.support_points[i](1) >= 0.0) &&
             (comp_dom.support_points[i](0) > comp_dom.boat_model.CurrentPointCenterTransom(0)-fabs(comp_dom.boat_model.CurrentPointCenterTransom(2)) ) &&
             (comp_dom.support_points[i](0) < comp_dom.boat_model.CurrentPointCenterTransom(0)+lh*fabs(comp_dom.boat_model.CurrentPointCenterTransom(2)) )   &&
             (comp_dom.support_points[i](2) > 2.0*comp_dom.boat_model.CurrentPointCenterTransom(2))    )
             {
             Point <3> dP0 = comp_dom.support_points[i];
             Point <3> dP; 
         				   //this is the vertical plane
             Handle(Geom_Plane) vertPlane = new Geom_Plane(0.,1.,0.,-dP0(1));
             vertPlane->Transform(comp_dom.boat_model.current_loc.Inverted());
             Handle(Geom_Curve) curve = comp_dom.boat_model.right_transom_bspline;

             GeomAPI_IntCS Intersector(curve, vertPlane);
             int npoints = Intersector.NbPoints();
             AssertThrow((npoints != 0), ExcMessage("Transom curve is not intersecting with vertical plane!"));
             //cout<<"Number of intersections: "<<npoints<<endl;
             double minDistance=1e7;
             for (int j=0; j<npoints;++j)
                 {
                 gp_Pnt int_point = Intersector.Point(j+1);
                 int_point.Transform(comp_dom.boat_model.current_loc);
                 Point<3> inters = Pnt(int_point);

                 if (dP0.distance(inters) < minDistance)
                    {
                    minDistance = dP0.distance(inters);
                    dP = inters;
                    }
                 }

             if ( (comp_dom.support_points[i](0) > dP(0)+comp_dom.min_diameter/20.0) &&
                  (comp_dom.support_points[i](0) < dP(0)+lh*fabs(dP(2))+comp_dom.min_diameter/20.0) )
                {
                //cout<<"££££ "<<dP<<endl;
                double mean_curvature;
                Point<3> normal;
                Point<3> projection;
                comp_dom.boat_model.boat_surface_right->normal_projection_and_diff_forms(projection,
                                                                                         normal,
                                                                                         mean_curvature,
			                                                                 dP);
                //cout<<dP0<<" curve--> "<<dP<<" surface--> "<<projection<<" ("<<dP.distance(projection)<<")"<<endl;
                AssertThrow((dP.distance(projection) < 1e-4*comp_dom.boat_model.boatWetLength), ExcMessage("Normal projection for surface normal evaluation went wrong!"));
                double transom_slope = -normal(0)/normal(2);
                double a = -transom_slope/(lh*fabs(dP(2))) - 1/dP(2)/lh/lh;
                double x = comp_dom.support_points[i](0)-dP(0);
                comp_dom.initial_map_points(3*i+2) = a*x*x + transom_slope*x + dP(2);
                comp_dom.map_points(3*i+2) = comp_dom.initial_map_points(3*i+2);
                //cout<<"****** "<<comp_dom.initial_map_points(3*i+2)<<" TS: "<<transom_slope<<"  a "<<a<<"  x "<<x<<endl;
                }
             }
             else if ( (comp_dom.support_points[i](1) > comp_dom.boat_model.CurrentPointLeftTransom(1)) &&
                       (comp_dom.support_points[i](1) < 0.0) &&
                       (comp_dom.support_points[i](0) > comp_dom.boat_model.CurrentPointCenterTransom(0)-fabs(comp_dom.boat_model.CurrentPointCenterTransom(2))) &&
                       (comp_dom.support_points[i](0) < comp_dom.boat_model.CurrentPointCenterTransom(0)+lh*fabs(comp_dom.boat_model.CurrentPointCenterTransom(2)))  &&
                       (comp_dom.support_points[i](2) > 2.0*comp_dom.boat_model.CurrentPointCenterTransom(2))   )
             {
             Point <3> dP0 = comp_dom.support_points[i];
             Point <3> dP; 
         				   //this is the vertical plane
             Handle(Geom_Plane) vertPlane = new Geom_Plane(0.,1.,0.,-dP0(1));
             vertPlane->Transform(comp_dom.boat_model.current_loc.Inverted());
             Handle(Geom_Curve) curve = comp_dom.boat_model.left_transom_bspline;

             GeomAPI_IntCS Intersector(curve, vertPlane);
             int npoints = Intersector.NbPoints();
             AssertThrow((npoints != 0), ExcMessage("Transom curve is not intersecting with vertical plane!"));
             //cout<<"Number of intersections: "<<npoints<<endl;
             double minDistance=1e7;
             for (int j=0; j<npoints;++j)
                 {
                 gp_Pnt int_point = Intersector.Point(j+1);
                 int_point.Transform(comp_dom.boat_model.current_loc);
                 Point<3> inters = Pnt(int_point);

                 if (dP0.distance(inters) < minDistance)
                    {
                    minDistance = dP0.distance(inters);
                    dP = inters;
                    }
                 }
             if ( (comp_dom.support_points[i](0) > dP(0)+comp_dom.min_diameter/20.0) &&
                  (comp_dom.support_points[i](0) < dP(0)+lh*fabs(dP(2))+comp_dom.min_diameter/20.0) )
                {
                double mean_curvature;
                Point<3> normal;
                Point<3> projection;
                comp_dom.boat_model.boat_surface_left->normal_projection_and_diff_forms(projection,
                                                                                        normal,
                                                                                        mean_curvature,
			                                                                dP);
                AssertThrow((dP.distance(projection) < 1e-4*comp_dom.boat_model.boatWetLength), ExcMessage("Normal projection for surface normal evaluation went wrong!"));
                double transom_slope = -normal(0)/normal(2);
                double a = -transom_slope/(lh*fabs(dP(2))) - 1/dP(2)/lh/lh;
                double x = comp_dom.support_points[i](0)-dP(0);
                comp_dom.initial_map_points(3*i+2) = a*x*x + transom_slope*x + dP(2);
                comp_dom.map_points(3*i+2) = comp_dom.initial_map_points(3*i+2);
                //cout<<"@@@@@@ "<<comp_dom.initial_map_points(3*i+2)<<" TS: "<<transom_slope<<"  a "<<a<<"  x "<<x<<endl;
                }
             }
        }
  comp_dom.vector_constraints.distribute(comp_dom.initial_map_points);
  comp_dom.vector_constraints.distribute(comp_dom.map_points);

  }




  double min_diameter = 1000000000;
  FEValues<dim-1,dim> fe_v(*comp_dom.mapping, comp_dom.fe, *comp_dom.quadrature,
			       update_JxW_values);
  const unsigned int DphiDt_n_q_points = fe_v.n_quadrature_points;
  cell_it
  cell = comp_dom.dh.begin_active(),
  endc = comp_dom.dh.end();

  for (; cell!=endc; ++cell)
      {
      if ((cell->material_id() == comp_dom.free_sur_ID1 ||
          cell->material_id() == comp_dom.free_sur_ID2 ||
          cell->material_id() == comp_dom.free_sur_ID3 ))
          {
          fe_v.reinit(cell);
          double area=0;
          for (unsigned int q=0; q<DphiDt_n_q_points; ++q)
              {
              area += fe_v.JxW(q);
              }
          min_diameter = fmin(min_diameter,2*sqrt(area/3.1415));
          }
       }

  cout<<"Min Diameter Corrected: "<<min_diameter<<endl;

  std::vector<Point<dim> > displaced_support_points(comp_dom.dh.n_dofs());
  DoFTools::map_dofs_to_support_points<dim-1, dim>( *comp_dom.mapping, comp_dom.dh, displaced_support_points);

   for (unsigned int i=0; i<comp_dom.dh.n_dofs(); i++)
      {
      dst(i+comp_dom.vector_dh.n_dofs()) = initial_wave_potential.value(displaced_support_points[i]);
      //std::cout<<i+comp_dom.vector_dh.n_dofs()<<" "<<dst(i+comp_dom.vector_dh.n_dofs())<<std::endl;
      }

   for (unsigned int i=0; i<comp_dom.dh.n_dofs(); i++)
      {
      dst(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()) =
                      initial_norm_potential_grad.value(displaced_support_points[i]);
      //std::cout<<i+comp_dom.vector_dh.n_dofs()<<" "<<dst(i+comp_dom.vector_dh.n_dofs())<<std::endl;
      }

//////////////TEST!/////////////////
   Vector<double> bem_phi(comp_dom.dh.n_dofs());
   Vector<double> bem_dphi_dn(comp_dom.dh.n_dofs());




     //cout<<"AFTER "<<endl;
     //for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
     //    if (constraints.is_constrained(i))
     //       cout<<i<<" "<<src_yy(i+comp_dom.dh.n_dofs()+comp_dom.vector_dh.n_dofs())<<endl;
     

     Vector<double> surface_nodes_backup = comp_dom.surface_nodes;
     Vector<double> other_nodes_backup = comp_dom.other_nodes;

     comp_dom.surface_nodes = 0.0;
     comp_dom.other_nodes = 1.0;

      comp_dom.compute_normals_at_nodes(comp_dom.map_points);
      for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
          {
          if ((comp_dom.flags[i] & water) || (comp_dom.flags[i] & boat))
             {
             bem_dphi_dn(i) = -comp_dom.node_normals[i]*Vinf;
             if ( comp_dom.flags[i] & transom_on_water )
                {
	        comp_dom.surface_nodes(i) = 0;
	        comp_dom.other_nodes(i) = 1;
                std::set<unsigned int> duplicates = comp_dom.double_nodes_set[i];
                duplicates.erase(i); 
                double bem_bc = 0;
                
                for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
                    {
                    bem_bc += bem_dphi_dn(*pos)/duplicates.size();
                    }
                bem_dphi_dn(i) = bem_bc;
                }
             }
          else if 
             (comp_dom.flags[i] & inflow)
             {
             bem_dphi_dn(i) = inflow_norm_potential_grad.value(support_points[i]);
             //cout<<i<<" "<<"   Point: "<<support_points[i]<<"   BC Val: "<<inflow_norm_potential_grad.value(support_points[i])<<endl;
             }
          else
             {
             comp_dom.surface_nodes(i) = 1.0;
             comp_dom.other_nodes(i) = 0.0;
             bem_phi(i) = 0.0;
             } 
          }   
     Vector<double> bem_bc(bem_dphi_dn);


     bem.solve(bem_phi, bem_dphi_dn, bem_bc);

     comp_dom.surface_nodes = surface_nodes_backup;
     comp_dom.other_nodes = other_nodes_backup;

     for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
         {
         dst(i+comp_dom.vector_dh.n_dofs()) = bem_phi(i); 
         dst(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()) = bem_dphi_dn(i); 
         }
     

////////////////////////////////////




//*/
  comp_dom.initial_map_points = comp_dom.map_points;
  vector_constraints.distribute(comp_dom.map_points);


  for(unsigned int i=0; i<comp_dom.vector_dh.n_dofs(); ++i)
      {
      dst(i) = comp_dom.map_points(i);
      //std::cout<<i<<" "<<dst(i)<<std::endl;
      } 


  for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
      {
      if ( ((comp_dom.flags[i] & boat) &&
           !(comp_dom.flags[i] & near_water) ) ||
           (comp_dom.flags[i] & transom_on_water) )
         {
         dst(3*i) -= restart_hull_displacement(0);
         dst(3*i+1) -= restart_hull_displacement(1);
         dst(3*i+2) -= restart_hull_displacement(2);         
         }
      }
  




// on the boat surface the normal potential must respect the boundary conditions (dphi_dn=-Vinf*n)
   std::vector<unsigned int> local_dof_indices (comp_dom.fe.dofs_per_cell);
    
      cell = comp_dom.dh.begin_active(),
      endc = comp_dom.dh.end();
      
      for (; cell!=endc; ++cell)
	{
	  cell->get_dof_indices(local_dof_indices);
	  if (cell->material_id() == comp_dom.wall_sur_ID1 ||
	      cell->material_id() == comp_dom.wall_sur_ID2 ||
	      cell->material_id() == comp_dom.wall_sur_ID3   )
	    for (unsigned int j=0; j<comp_dom.fe.dofs_per_cell; ++j) 
	      {
		unsigned int id=local_dof_indices[j];
                ////////////////////////////////////
                // needs to be changed
                //Point<3> original(displaced_support_points[id](0),fabs(displaced_support_points[id](1)),displaced_support_points[id](2));
                //Point<3> projection;
                //Point<3> normal;
                //double mean_curvature;
                //comp_dom.boat_model.boat_water_line_right->axis_projection_and_diff_forms(projection, normal, mean_curvature, original);
                //double b = displaced_support_points[id](1) < 0.0 ? -1.0 : 1.0; 
                //projection(1)*=b;
                //normal(1)*=b;

                //BoatSurface<3> boat_surface;
                //Point<3> normal2 = boat_surface.HullNormal(displaced_support_points[id]);
                //std::cout<<std::endl;
                //std::cout<<"Point "<<displaced_support_points[id]<<std::endl;
                //std::cout<<"NormalNew "<<normal<<std::endl;
                //std::cout<<"NormalExact "<<normal2<<std::endl;
                ////////////////////////////////////
		//dst(id+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()) = -comp_dom.iges_normals[id]*Vinf;
                Point<3> Vhull(dst(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()),
                               dst(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+1),
                               dst(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+2));
                dst(id+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()) = -comp_dom.node_normals[id]*(Vinf-Vhull); 
	      }
	}
    

   
   max_x_coor_value = 0;
   max_y_coor_value = 0;
   max_z_coor_value = 0;
   for (unsigned int i=0; i < comp_dom.dh.n_dofs(); i++)
       {
       //for printout
       //std::cout<<"Node "<<i<< "["<<support_points(i)<<"] "<<std::endl;
       max_x_coor_value = std::max(max_x_coor_value,std::abs(displaced_support_points[i](0)));
       max_y_coor_value = std::max(max_y_coor_value,std::abs(displaced_support_points[i](1)));
       max_z_coor_value = std::max(max_z_coor_value,std::abs(displaced_support_points[i](2)));
       }

   DphiDt_sys_solution.reinit (comp_dom.dh.n_dofs());
   DphiDt_sys_solution_2.reinit (comp_dom.dh.n_dofs());
   DphiDt_sys_solution_3.reinit (comp_dom.dh.n_dofs());
   break_wave_press.reinit (comp_dom.dh.n_dofs());
   
   vector_sys_solution.reinit (comp_dom.vector_dh.n_dofs());
   vector_sys_solution_2.reinit (comp_dom.vector_dh.n_dofs());

   Vector<double> dummy_sol_dot(dst.size());




  std::string filename = ( output_file_name + "_" +
			   Utilities::int_to_string(0) +
			   ".vtu" );
  //Vector<double> dummy_sol_dot(dst.size());
  output_results(filename, 0, dst, dummy_sol_dot);

  comp_dom.old_map_points = comp_dom.map_points;
  comp_dom.rigid_motion_map_points = 0;     

  last_remesh_time = 0;

  current_sol = dst;
  current_sol_dot = 0.0;
  current_time = 0.0;
  
  

  restart_flag = false;


     // we have to compute the reference transom wet surface with the new mesh here
   if (!comp_dom.no_boat && comp_dom.boat_model.is_transom)
      {
        double transom_wet_surface = 0;
      std::vector<Point<3> > vertices;
      std::vector<CellData<2> > cells;
      for (tria_it elem=comp_dom.tria.begin_active(); elem!= comp_dom.tria.end();++elem)
          {
	  if ((elem->material_id() == comp_dom.wall_sur_ID1 ||
	       elem->material_id() == comp_dom.wall_sur_ID2 ||
	       elem->material_id() == comp_dom.wall_sur_ID3 )) 
	     {
	     if (elem->at_boundary())
                {
                for (unsigned int f=0; f<GeometryInfo<2>::faces_per_cell; ++f)
		    if ( elem->face(f)->boundary_indicator() == 32 ||
                         elem->face(f)->boundary_indicator() == 37 )
                       {
                       unsigned int index_0 = comp_dom.find_point_id(elem->face(f)->vertex(0),comp_dom.ref_points);
                       unsigned int index_1 = comp_dom.find_point_id(elem->face(f)->vertex(1),comp_dom.ref_points);
                       if (comp_dom.ref_points[3*index_1](1) < comp_dom.ref_points[3*index_0](1))
                          {
                          vertices.push_back(comp_dom.ref_points[3*index_0]);
                          vertices.push_back(comp_dom.ref_points[3*index_1]);
                          vertices.push_back(comp_dom.ref_points[3*index_1]+
                                             -1.0*Point<3>(0.0,0.0,comp_dom.ref_points[3*index_1](2)));
                          vertices.push_back(comp_dom.support_points[index_0]+
                                             -1.0*Point<3>(0.0,0.0,comp_dom.ref_points[3*index_0](2)));
                          cells.resize(cells.size()+1);
                          cells[cells.size()-1].vertices[0]=4*(cells.size()-1)+0;
                          cells[cells.size()-1].vertices[1]=4*(cells.size()-1)+1;
                          cells[cells.size()-1].vertices[2]=4*(cells.size()-1)+2;
                          cells[cells.size()-1].vertices[3]=4*(cells.size()-1)+3;
                          }
                       else
                          {
                          vertices.push_back(comp_dom.support_points[index_1]);
                          vertices.push_back(comp_dom.support_points[index_0]);
                          vertices.push_back(comp_dom.support_points[index_0]+
                                             -1.0*Point<3>(0.0,0.0,comp_dom.ref_points[3*index_0](2)));
                          vertices.push_back(comp_dom.support_points[index_1]+
                                             -1.0*Point<3>(0.0,0.0,comp_dom.ref_points[3*index_1](2)));
                          cells.resize(cells.size()+1);
                          cells[cells.size()-1].vertices[0]=4*(cells.size()-1)+0;
                          cells[cells.size()-1].vertices[1]=4*(cells.size()-1)+1;
                          cells[cells.size()-1].vertices[2]=4*(cells.size()-1)+2;
                          cells[cells.size()-1].vertices[3]=4*(cells.size()-1)+3;
                          }                          
                       }
	        }
             }
          }


      SubCellData subcelldata;
      Triangulation<dim-1, dim> transom_tria;
      GridTools::delete_unused_vertices (vertices, cells, subcelldata);
      GridReordering<2,3>::reorder_cells (cells);
      transom_tria.create_triangulation_compatibility(vertices, cells, subcelldata );
      FE_Q<dim-1,dim> transom_fe(1);
      DoFHandler<dim-1,dim> transom_dh(transom_tria);
      transom_dh.distribute_dofs(transom_fe);

      FEValues<dim-1,dim> transom_fe_v(transom_fe, *comp_dom.quadrature,
			               update_values | update_gradients |
			               update_cell_normal_vectors |
			               update_quadrature_points |
			               update_JxW_values);

      const unsigned int transom_n_q_points = transom_fe_v.n_quadrature_points;
   

      std::vector<double> transom_pressure_quad_values(transom_n_q_points);
      for (cell_it cell = transom_dh.begin_active(); cell!=transom_dh.end(); ++cell)
          {
          transom_fe_v.reinit(cell);
     
          for (unsigned int q=0; q<transom_n_q_points; ++q)
              {
              transom_wet_surface += 1 * transom_fe_v.JxW(q);
              }
          }
      ref_transom_wet_surface = transom_wet_surface;
      }

      // here we finally compute the reference cell areas, which are needed
      // to evaluate boat cell deformations and possibly activate smoothing
      //FEValues<dim-1,dim> fe_v(*comp_dom.mapping, comp_dom.fe, *comp_dom.quadrature,
	//		       update_JxW_values);
      //const unsigned int DphiDt_n_q_points = fe_v.n_quadrature_points;
      ref_cell_areas.reinit(comp_dom.tria.n_active_cells());
      cell = comp_dom.dh.begin_active(),
      endc = comp_dom.dh.end();

      unsigned int count=0;
      for (; cell!=endc; ++cell)
          {
          if ((cell->material_id() == comp_dom.wall_sur_ID1 ||
               cell->material_id() == comp_dom.wall_sur_ID2 ||
               cell->material_id() == comp_dom.wall_sur_ID3 ))
             {
             fe_v.reinit(cell);
             for (unsigned int q=0; q<DphiDt_n_q_points; ++q)
                 {
                 ref_cell_areas(count) += 1.0 * fe_v.JxW(q);
                 }
             }
          ++count;
          }

cout<<"££££2 "<<comp_dom.map_points[2067]<<endl;

} 


template <int dim>
unsigned int FreeSurface<dim>:: n_dofs() const{

    return dofs_number;

}

template <int dim>
void FreeSurface<dim>::output_step(Vector<double> & solution,
	                           const unsigned int step_number)
{
   std::cout<<"iteration: "<<step_number<<std::endl;
   std::cout<<std::endl; 
 
   std::string filename = ( output_file_name + "_" +
  		            Utilities::int_to_string(step_number+1) +
			    ".vtu" );    
   output_results(filename, current_time, solution, current_sol_dot);

   ofstream waterLineFile;
   waterLineFile.open ((output_file_name+"_water_line.txt").c_str());
   for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
       if ((comp_dom.flags[i] & water) && (comp_dom.flags[i] & near_boat) && (!(comp_dom.flags[i] & transom_on_water)))
          waterLineFile<<comp_dom.support_points[i]<<endl;
   waterLineFile.close();
              
}


template <int dim>
void FreeSurface<dim>::output_step(Vector<double> & solution,
				   Vector<double> &solution_dot,
				   const double t,
	                           const unsigned int step_number,
		                   const double  h)
{
   std::cout<<"t = "<<t<<"   TS = "<<step_number<<std::endl;
   std::cout<<std::endl; 
 
   std::string filename = ( output_file_name + "_" +
  		            Utilities::int_to_string(step_number+1) +
			    ".vtu" );    
   output_results(filename, t, solution, solution_dot);

   ofstream waterLineFile;
   waterLineFile.open ((output_file_name+"_water_line.txt").c_str());
   for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
       if ((comp_dom.flags[i] & water) && (comp_dom.flags[i] & near_boat) && (!(comp_dom.flags[i] & transom_on_water)))
          waterLineFile<<comp_dom.support_points[i]<<endl;
   waterLineFile.close();
              
}



typedef std::pair<double,unsigned int> mypair;
bool comparator ( const mypair& l, const mypair& r)
   { return l.first < r.first; }

template <int dim>
bool FreeSurface<dim>::solution_check(Vector<double> & solution,
				      Vector<double> &solution_dot,
				      const double t,
	                              const unsigned int step_number,
		                      const double  h)
{



  static unsigned int remeshing_counter = 1;
  static unsigned int dumping_counter = 1;
  
  if( t>=dumping_period*dumping_counter ) 
    {
      std::string fname = output_file_name + "_dump_" +
			  Utilities::int_to_string(dumping_counter);
      ++dumping_counter;
      dump_solution(solution, solution_dot, fname);
    }
		    
  if( t>=remeshing_period*remeshing_counter && comp_dom.dh.n_dofs() < max_number_of_dofs )
    {
      ++remeshing_counter;
      std::cout<<"Checking and updating mesh..."<<std::endl;


      // first thing to do is update the geometry to the current positon
      Point<3> current_hull_displacement;
      Point<3> current_hull_displacement_dot;
      Point<3> current_hull_velocity;
      Point<3> current_hull_velocity_dot;
      Point<3> current_hull_ang_velocity;
      Point<3> current_hull_ang_velocity_dot;
      Point<3> current_hull_quat_vector;
      Point<3> current_hull_quat_vector_dot;
      for (unsigned int k=0; k<3; ++k)
          {
          current_hull_displacement(k) = solution(k+3+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs());
          current_hull_displacement_dot(k) = solution_dot(k+3+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs());
          current_hull_velocity(k) = solution(k+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs());
          current_hull_velocity_dot(k) = solution_dot(k+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs());
          current_hull_ang_velocity(k) = solution(k+6+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs());
          current_hull_ang_velocity_dot(k) = solution_dot(k+6+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs());
          current_hull_quat_vector(k) = solution(k+9+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs());
          current_hull_quat_vector_dot(k) = solution_dot(k+9+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs());
          }
      double current_hull_quat_scalar = solution(12+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs());
      double current_hull_quat_scalar_dot = solution_dot(12+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs());

      restart_hull_displacement(0) = current_hull_displacement(0);
      restart_hull_displacement(1) = current_hull_displacement(1);
      restart_hull_displacement(2) = current_hull_displacement(2);
      restart_hull_quat_vector(0) = current_hull_quat_vector(0);
      restart_hull_quat_vector(1) = current_hull_quat_vector(1);
      restart_hull_quat_vector(2) = current_hull_quat_vector(2);
      restart_hull_quat_scalar = current_hull_quat_scalar;


     // this moves the rigid_motion_map_points vector (target positions) for the internal nodes of the hull 
     // mesh (all nodes except for the free surface ones)
     Vector<double> rigid_motion_velocities(comp_dom.vector_dh.n_dofs());

     // this moves the rigid_motion_map_points vector (target positions) for the internal nodes of the hull 
     // mesh (all nodes except for the free surface ones)
     for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
         {
         if ( ((comp_dom.flags[i] & boat) &&
              !(comp_dom.flags[i] & near_water) ) ||
              (comp_dom.flags[i] & transom_on_water) )
            {
            //if (fabs(t-0.2) <1e-5)
            //cout<<"BEFORE: "<<comp_dom.ref_points[3*i](0)+comp_dom.old_map_points(3*i)<<" "
            //                <<comp_dom.ref_points[3*i](1)+comp_dom.old_map_points(3*i+1)<<" "
            //                <<comp_dom.ref_points[3*i](2)+comp_dom.old_map_points(3*i+2)<<endl;
            gp_Pnt original_boat_mesh_point = Pnt(Point<3>(comp_dom.ref_points[3*i](0)+comp_dom.old_map_points(3*i),
                                                           comp_dom.ref_points[3*i](1)+comp_dom.old_map_points(3*i+1),
                                                           comp_dom.ref_points[3*i](2)+comp_dom.old_map_points(3*i+2)));
            gp_Pnt boat_mesh_point = original_boat_mesh_point;
            // we first take this point (which is in the RESTART hull location) and transform it to be in the
            // REFERENCE configuration 
            boat_mesh_point.Transform(restart_hull_location.Inverted());

            // now we use sacado to compute the residual at this dof, along with its derivatives with respect to 
            // the 7 dofs associated to the rigid linear and angular displacements

            double s_x = restart_hull_displacement(0);
            double s_y = restart_hull_displacement(1);
            double s_z = restart_hull_displacement(2);
            double v_x = restart_hull_quat_vector(0);
            double  v_y = restart_hull_quat_vector(1);
            double v_z = restart_hull_quat_vector(2);
            double s = restart_hull_quat_scalar;

            Point<3> ref_baricenter_pos(comp_dom.boat_model.reference_hull_baricenter(0),
                                        comp_dom.boat_model.reference_hull_baricenter(1),
                                        comp_dom.boat_model.reference_hull_baricenter(2));
   
            Point<3> baricenter_pos(s_x+comp_dom.boat_model.reference_hull_baricenter(0),
                                    s_y+comp_dom.boat_model.reference_hull_baricenter(1),
                                    s_z+comp_dom.boat_model.reference_hull_baricenter(2));

            Point<3> RotMatRow1(1-2*v_y*v_y-2*v_z*v_z, 2*v_x*v_y-2*s*v_z, 2*v_x*v_z+2*s*v_y);
            Point<3> RotMatRow2(2*v_x*v_y+2*s*v_z, 1-2*v_x*v_x-2*v_z*v_z, 2*v_y*v_z-2*s*v_x);
            Point<3> RotMatRow3(2*v_x*v_z-2*s*v_y, 2*v_y*v_z+2*s*v_x, 1-2*v_y*v_y-2*v_x*v_x);

            Point<3> ref_point_pos(boat_mesh_point.X(),boat_mesh_point.Y(),boat_mesh_point.Z());
            Point<3> rigid_lin_displ(s_x,s_y,s_z);
            Point<3> target_point_pos(RotMatRow1*(ref_point_pos+(-1.0)*ref_baricenter_pos),
                                      RotMatRow2*(ref_point_pos+(-1.0)*ref_baricenter_pos),
                                      RotMatRow3*(ref_point_pos+(-1.0)*ref_baricenter_pos));
            target_point_pos += baricenter_pos;
            // now we have the point on the REFERENCE hull, and transform it to go onto the CURRENT hull
            //boat_mesh_point.Transform(reference_to_current_transformation);
            // the rigid motion map points is the difference between the reference point position and the
            // current (rigidly displaced) node position
            comp_dom.rigid_motion_map_points(3*i) = target_point_pos(0)-ref_point_pos(0);
            comp_dom.rigid_motion_map_points(3*i+1) = target_point_pos(1)-ref_point_pos(1);
            comp_dom.rigid_motion_map_points(3*i+2) = target_point_pos(2)-ref_point_pos(2);

            rigid_motion_velocities(3*i) = current_hull_velocity(0)+
                                           current_hull_ang_velocity(1)*(target_point_pos(2)-baricenter_pos(2))-
                                           current_hull_ang_velocity(2)*(target_point_pos(1)-baricenter_pos(1));
            rigid_motion_velocities(3*i+1) = current_hull_velocity(1)+
                                             current_hull_ang_velocity(2)*(target_point_pos(0)-baricenter_pos(0))-
                                             current_hull_ang_velocity(0)*(target_point_pos(2)-baricenter_pos(2));
            rigid_motion_velocities(3*i+2) = current_hull_velocity(2)+
                                             current_hull_ang_velocity(0)*(target_point_pos(1)-baricenter_pos(1))-
                                             current_hull_ang_velocity(1)*(target_point_pos(0)-baricenter_pos(0));
            

            //if (fabs(t-0.2) <1e-5)
            //cout<<"AFTER: "<<Pnt(boat_mesh_point)<<" vs "<<Pnt(original_boat_mesh_point)<<endl;
            //cout<<"RMMP: "<<comp_dom.rigid_motion_map_points(3*i)<<" "
            //              <<comp_dom.rigid_motion_map_points(3*i+1)<<" "
            //              <<comp_dom.rigid_motion_map_points(3*i+2)<<endl;
            //cout<<"NN: "<<comp_dom.rigid_motion_map_points(3*i)+nodes_positions(3*i)+comp_dom.ref_points[3*i](0)<<" "
            //            <<comp_dom.rigid_motion_map_points(3*i+1)+nodes_positions(3*i+1)+comp_dom.ref_points[3*i](1)<<" "
            //            <<comp_dom.rigid_motion_map_points(3*i+2)+nodes_positions(3*i+2)+comp_dom.ref_points[3*i](2)<<endl;
            }
         }

      restart_hull_location = comp_dom.boat_model.set_current_position(restart_hull_displacement,
                                                                       restart_hull_quat_scalar,
                                                                       restart_hull_quat_vector);
      restart_transom_center_point = comp_dom.boat_model.CurrentPointCenterTransom; 
      restart_transom_left_point = comp_dom.boat_model.CurrentPointLeftTransom;
      restart_transom_right_point = comp_dom.boat_model.CurrentPointRightTransom;
      restart_transom_left_tangent = comp_dom.boat_model.current_left_transom_tangent;
      restart_transom_right_tangent = comp_dom.boat_model.current_right_transom_tangent;

      //std::string filename1 = ( "preRemesh.vtu" );
      //output_results(filename1, t, solution, solution_dot);

      Triangulation<dim-1, dim> &tria = comp_dom.tria;
      DoFHandler<dim-1, dim> &dh = comp_dom.dh;
      DoFHandler<dim-1, dim> &vector_dh = comp_dom.vector_dh;
    
      std::vector<Point<dim> > support_points(comp_dom.dh.n_dofs());
      DoFTools::map_dofs_to_support_points<dim-1, dim>( *comp_dom.mapping, comp_dom.dh, support_points);

      Vector <double> disp_vector(comp_dom.vector_dh.n_dofs());
      for (unsigned int i=0; i<comp_dom.dh.n_dofs();++i)
          {
          for (unsigned int j=0; j<dim; ++j)
              disp_vector(i*dim+j) = support_points[i][j]; 
          }

      VectorView<double> Phi(dh.n_dofs(), solution.begin()+vector_dh.n_dofs());
      VectorView<double> Phi_dot(dh.n_dofs(), solution_dot.begin()+vector_dh.n_dofs());
      VectorView<double> dphi_dn(dh.n_dofs(), solution.begin()+vector_dh.n_dofs()+dh.n_dofs());
      VectorView<double> dphi_dn_dot(dh.n_dofs(), solution_dot.begin()+vector_dh.n_dofs()+dh.n_dofs());

      std::vector<Point<3> > old_points(dh.n_dofs());
      old_points = support_points;

      Vector<double> curvatures(vector_dh.n_dofs());
      comp_dom.compute_curvatures(curvatures);
      
    
      Vector<float> estimated_error_per_cell(tria.n_active_cells());

      QGauss<dim-2> quad(2);

      VectorView<double> displacements(vector_dh.n_dofs(), solution.begin());
      Vector <double> elevations(dh.n_dofs());

      for (unsigned int i=0; i<dh.n_dofs();++i)
          if (comp_dom.flags[i] & water)
             elevations(i) = displacements(3*i+2);
//      KellyErrorEstimator<dim-1,dim>::estimate (dh,
// 						quad,
// 						typename FunctionMap<dim>::type(),
// 						elevations,
// 						estimated_error_per_cell);


//     if (t>4.5)
//        comp_dom.n_cycles = 0;

     if (t < 0.0)
         {
         VectorView<double> displacements_dot(vector_dh.n_dofs(), solution_dot.begin());
         KellyErrorEstimator<dim-1,dim>::estimate (vector_dh,
						   quad,
						   typename FunctionMap<dim>::type(),
						   (const Vector<double>&)displacements_dot,
	 					   estimated_error_per_cell);
         }
     else
         {
//         KellyErrorEstimator<dim-1,dim>::estimate (vector_dh,
//						   quad,
//						   typename FunctionMap<dim>::type(),
//						   (const Vector<double>&)displacements,
//	 					   estimated_error_per_cell);
         KellyErrorEstimator<dim-1,dim>::estimate (dh,
 	  					   quad,
 						   typename FunctionMap<dim>::type(),
 						   elevations,
 						   estimated_error_per_cell);        
         }


     double max_boat_error = 0;
     double max_other_error = 0;
     unsigned int counter=0;
     for (cell_it elem=dh.begin_active(); elem!= dh.end();++elem)
	 {
	 if ((elem->material_id() == comp_dom.wall_sur_ID1 ||
	      elem->material_id() == comp_dom.wall_sur_ID2 ||
              elem->material_id() == comp_dom.wall_sur_ID3 ))
             {
             max_boat_error = fmax(max_boat_error,estimated_error_per_cell(counter)); 
             }
         else
             max_other_error = fmax(max_other_error,estimated_error_per_cell(counter));
         ++counter;
         }
 
     counter = 0;
     for (cell_it elem=dh.begin_active(); elem!= dh.end();++elem)
	 {
	 if ((elem->material_id() == comp_dom.wall_sur_ID1 ||
	      elem->material_id() == comp_dom.wall_sur_ID2 ||
	      elem->material_id() == comp_dom.wall_sur_ID3 ))
             estimated_error_per_cell(counter) /= fmax(max_boat_error,1e-6);
         else
             estimated_error_per_cell(counter) /= fmax(max_other_error,1e-6);
         ++counter;
         }


      SolutionTransfer<dim-1, Vector<double>, DoFHandler<dim-1, dim> > soltrans(dh);
 

      Vector<double> positions_x(comp_dom.dh.n_dofs());
      Vector<double> positions_y(comp_dom.dh.n_dofs());
      Vector<double> positions_z(comp_dom.dh.n_dofs());      
      Vector<double> positions_dot_x(comp_dom.dh.n_dofs());
      Vector<double> positions_dot_y(comp_dom.dh.n_dofs());
      Vector<double> positions_dot_z(comp_dom.dh.n_dofs());

      for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
          {
          positions_x(i) = solution(3*i+0)+comp_dom.rigid_motion_map_points(3*i+0)+comp_dom.ref_points[3*i](0);
          positions_y(i) = solution(3*i+1)+comp_dom.rigid_motion_map_points(3*i+1)+comp_dom.ref_points[3*i](1);
          positions_z(i) = solution(3*i+2)+comp_dom.rigid_motion_map_points(3*i+2)+comp_dom.ref_points[3*i](2);
          positions_dot_x(i) = solution_dot(3*i+0)+rigid_motion_velocities(3*i+0);
          positions_dot_y(i) = solution_dot(3*i+1)+rigid_motion_velocities(3*i+1);
          positions_dot_z(i) = solution_dot(3*i+2)+rigid_motion_velocities(3*i+2);
          }          


      std::vector<Vector<double> > all_in;
      all_in.push_back((Vector<double>)Phi);
      all_in.push_back((Vector<double>)Phi_dot);
      all_in.push_back((Vector<double>)dphi_dn);
      all_in.push_back((Vector<double>)dphi_dn_dot);
      all_in.push_back(positions_x);
      all_in.push_back(positions_y);
      all_in.push_back(positions_z);
      all_in.push_back(positions_dot_x);
      all_in.push_back(positions_dot_y);
      all_in.push_back(positions_dot_z); 
 

      GridRefinement::refine_and_coarsen_fixed_number	(tria,
                                                         estimated_error_per_cell,
                                                         refinement_fraction,
							 coarsening_fraction,
							 max_number_of_dofs);
//      GridRefinement::refine_and_coarsen_optimize	(tria,
//                                                         estimated_error_per_cell);



      // here we compute the cell diameters, which are needed
      // to limit cell refinement to bigger cells
      FEValues<dim-1,dim> fe_v(*comp_dom.mapping, comp_dom.fe, *comp_dom.quadrature,
	   		       update_JxW_values);
      const unsigned int n_q_points = fe_v.n_quadrature_points;
      Vector<double> cell_diameters(comp_dom.tria.n_active_cells());
      cell_it
      cell = comp_dom.dh.begin_active(),
      endc = comp_dom.dh.end();

      counter=0;
      for (; cell!=endc; ++cell)
          {
          fe_v.reinit(cell);
          for (unsigned int q=0; q<n_q_points; ++q)
              {
              cell_diameters(counter) += fe_v.JxW(q);
              }
          cell_diameters(counter) = 2*sqrt(cell_diameters(counter)/3.1415);
          ++counter;
          }

      counter=0;
      if(comp_dom.n_cycles == 0)
	{
	for (cell_it elem=dh.begin_active(); elem!= dh.end();++elem)
	    {
	    if(cell_diameters(counter)*8.0/2.0 < comp_dom.min_diameter)
	      elem->clear_refine_flag();

            if(fabs(elem->center()(0)) > 20.0)
	      elem->clear_refine_flag();
	    
	    if ((elem->material_id() != comp_dom.wall_sur_ID1) &&
		(elem->material_id() != comp_dom.wall_sur_ID2) &&
		(elem->material_id() != comp_dom.wall_sur_ID3) &&
                (elem->material_id() != comp_dom.free_sur_ID1) &&
		(elem->material_id() != comp_dom.free_sur_ID2) &&
		(elem->material_id() != comp_dom.free_sur_ID3)) 
	      {
	      elem->clear_refine_flag();
              elem->clear_coarsen_flag();
	      }
	    
	    if ((elem->material_id() == comp_dom.free_sur_ID1 ||
		 elem->material_id() == comp_dom.free_sur_ID2 ||
		 elem->material_id() == comp_dom.free_sur_ID3 )) 
	      {
	      if ((elem->center().distance(Point<3>(0.0,0.0,0.0)) > 10.0) &&
                   elem->at_boundary() )
		 elem->clear_refine_flag();
		 elem->clear_coarsen_flag();
	      }
            counter++;
            }

	  typedef typename Triangulation<dim-1, dim>::active_cell_iterator tria_it;
	  typename std::map<tria_it, tria_it>::iterator it;
	  for(it = comp_dom.boat_to_water_edge_cells.begin();
	      it != comp_dom.boat_to_water_edge_cells.end(); ++it) 
	    {
	      if( it->first->refine_flag_set() &&
		  !it->second->refine_flag_set() )
		{ 
                  for (unsigned int i=0; i<it->second->parent()->n_children(); ++i)
		      //if(!it->second->parent()->child(i).has_children())
                        it->second->parent()->child(i)->clear_coarsen_flag();
		  it->second->set_refine_flag();
		}
	      else if( it->second->refine_flag_set() &&
		       !it->first->refine_flag_set() ) 
		{
                  for (unsigned int i=0; i<it->first->parent()->n_children(); ++i)
                      //if(!it->first->parent()->child(i).has_children())
		        it->first->parent()->child(i)->clear_coarsen_flag();
		  it->first->set_refine_flag();
		}
	      else if( it->first->coarsen_flag_set() &&
		       !it->second->coarsen_flag_set() )
		{
                  for (unsigned int i=0; i<it->first->parent()->n_children(); ++i)
                      //if(!it->first->parent()->child(i).has_children())
		        it->first->parent()->child(i)->clear_coarsen_flag();
		}
	      else if( it->second->coarsen_flag_set() &&
		       !it->first->coarsen_flag_set() ) 
		{
		  for (unsigned int i=0; i<it->second->parent()->n_children(); ++i)
                      //if(!it->second->parent()->child(i).has_children())
		        it->second->parent()->child(i)->clear_coarsen_flag();
		}
	    }
	}
      else
	for (cell_it elem=dh.begin_active(); elem!= dh.end();++elem)
	  {

            if (cell_diameters(counter)/(adaptive_ref_limit) < comp_dom.min_diameter)
               {
	       elem->clear_refine_flag();
               }
            
            if(fabs(elem->center()(0)) > comp_dom.boat_model.boatWetLength*5.0)
	       elem->clear_refine_flag();
	    
	    if ((elem->material_id() == comp_dom.wall_sur_ID1 ||
		 elem->material_id() == comp_dom.wall_sur_ID2 ||
		 elem->material_id() == comp_dom.wall_sur_ID3 )) 
	      {
		elem->clear_refine_flag();
		elem->clear_coarsen_flag();
	      }
	    
	    if ((elem->material_id() == comp_dom.free_sur_ID1 ||
		 elem->material_id() == comp_dom.free_sur_ID2 ||
		 elem->material_id() == comp_dom.free_sur_ID3 )) 
	      {
	      if (elem->at_boundary())
                 {
                 bool clear = true;
                 for(unsigned int f=0; f<GeometryInfo<2>::faces_per_cell; ++f)
		    if( elem->face(f)->boundary_indicator() == 40 ||
                        elem->face(f)->boundary_indicator() == 41 ||
                        elem->face(f)->boundary_indicator() == 26 || //this allows for refinement of waterline cells on water side
                        elem->face(f)->boundary_indicator() == 27 || //this allows for refinement of waterline cells on water side
                        elem->face(f)->boundary_indicator() == 28 || //this allows for refinement of waterline cells on water side
                        elem->face(f)->boundary_indicator() == 29  //this allows for refinement of waterline cells on water side
                      )
                      {
                      clear = false;
                      //if ( ( elem->face(f)->boundary_indicator() == 26 || //this blocks last step of refinement of waterline cells if close to stern
                      //       elem->face(f)->boundary_indicator() == 27 || //this blocks last step of refinement of waterline cells if close to stern
                      //       elem->face(f)->boundary_indicator() == 28 || //this blocks last step of refinement of waterline cells if close to stern
                      //       elem->face(f)->boundary_indicator() == 29 ) && //this blocks last step of refinement of waterline cells if close to stern
                      //     (elem->center()(0) > 1.70) &&
                      //     (elem->diameter()/4.0 < comp_dom.min_diameter) )
                      //   clear = true;
                      }
                 if (clear)
                    { 
		    elem->clear_refine_flag();
		    elem->clear_coarsen_flag();
                    }
                 }
	      }
            else
              {
              elem->clear_refine_flag();
              elem->clear_coarsen_flag();
              }

	  counter++;
          }
 

				       // prepare the triangulation,
      tria.prepare_coarsening_and_refinement();
      soltrans.prepare_for_coarsening_and_refinement(all_in);



      tria.execute_coarsening_and_refinement();


      dh.distribute_dofs(comp_dom.fe);
      vector_dh.distribute_dofs(comp_dom.vector_fe);  
      comp_dom.map_points.reinit(vector_dh.n_dofs());
      comp_dom.smoothing_map_points.reinit(vector_dh.n_dofs());
      comp_dom.old_map_points.reinit(vector_dh.n_dofs());
      comp_dom.rigid_motion_map_points.reinit(vector_dh.n_dofs());
      comp_dom.initial_map_points.reinit(vector_dh.n_dofs());
      comp_dom.ref_points.resize(vector_dh.n_dofs());
      DoFTools::map_dofs_to_support_points<2,3>(StaticMappingQ1<2,3>::mapping,
					      comp_dom.vector_dh, comp_dom.ref_points);
      comp_dom.generate_double_nodes_set();

      
      compute_constraints(constraints, vector_constraints);


      this->dofs_number = vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()+13;

				       //DoFRenumbering::Cuthill_McKee(dh);
				       //DoFRenumbering::Cuthill_McKee(vector_dh);
      std::cout<<"Total number of dofs after refinement: "<<dh.n_dofs()<<std::endl;


      Vector<double> new_curvatures(vector_dh.n_dofs());

      Vector<double> new_Phi(dh.n_dofs());
      Vector<double> new_Phi_dot(dh.n_dofs());

      Vector <double> new_dphi_dn(dh.n_dofs());     
      Vector <double> new_dphi_dn_dot(dh.n_dofs());

      Vector<double> new_positions_x(dh.n_dofs());
      Vector<double> new_positions_y(dh.n_dofs());
      Vector<double> new_positions_z(dh.n_dofs());
      Vector<double> new_positions_dot_x(dh.n_dofs());
      Vector<double> new_positions_dot_y(dh.n_dofs());
      Vector<double> new_positions_dot_z(dh.n_dofs());

      std::vector<Vector<double> > all_out;
      all_out.push_back(new_Phi);
      all_out.push_back(new_Phi_dot);
      all_out.push_back(new_dphi_dn);
      all_out.push_back(new_dphi_dn_dot);
      all_out.push_back(new_positions_x);
      all_out.push_back(new_positions_y);
      all_out.push_back(new_positions_z);
      all_out.push_back(new_positions_dot_x);
      all_out.push_back(new_positions_dot_y);
      all_out.push_back(new_positions_dot_z);
  
      soltrans.interpolate(all_in, all_out);
  
      solution.reinit(dofs_number);
      solution_dot.reinit(dofs_number); 

      constraints.distribute(all_out[0]);
      constraints.distribute(all_out[1]);
      constraints.distribute(all_out[2]);
      constraints.distribute(all_out[3]);
      constraints.distribute(all_out[4]);
      constraints.distribute(all_out[5]);
      constraints.distribute(all_out[6]);
      constraints.distribute(all_out[7]);    
      constraints.distribute(all_out[8]);
      constraints.distribute(all_out[9]);

     // we have to compute rigid_motion_map_points and rigid_motion_velocities on the new mesh
     rigid_motion_velocities.reinit(comp_dom.vector_dh.n_dofs());

     // this moves the rigid_motion_map_points vector (target positions) for the internal nodes of the hull 
     // mesh (all nodes except for the free surface ones)
     for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
         {
         if ( ((comp_dom.flags[i] & boat) &&
              !(comp_dom.flags[i] & near_water) ) ||
              (comp_dom.flags[i] & transom_on_water) )
            {
            //if (fabs(t-0.2) <1e-5)
            //cout<<"BEFORE: "<<comp_dom.ref_points[3*i](0)+comp_dom.old_map_points(3*i)<<" "
            //                <<comp_dom.ref_points[3*i](1)+comp_dom.old_map_points(3*i+1)<<" "
            //                <<comp_dom.ref_points[3*i](2)+comp_dom.old_map_points(3*i+2)<<endl;
            gp_Pnt original_boat_mesh_point = Pnt(Point<3>(comp_dom.ref_points[3*i](0)+comp_dom.old_map_points(3*i),
                                                           comp_dom.ref_points[3*i](1)+comp_dom.old_map_points(3*i+1),
                                                           comp_dom.ref_points[3*i](2)+comp_dom.old_map_points(3*i+2)));
            gp_Pnt boat_mesh_point = original_boat_mesh_point;
            // we first take this point (which is in the RESTART hull location) and transform it to be in the
            // REFERENCE configuration 
            boat_mesh_point.Transform(restart_hull_location.Inverted());

            // now we use sacado to compute the residual at this dof, along with its derivatives with respect to 
            // the 7 dofs associated to the rigid linear and angular displacements

            double s_x = restart_hull_displacement(0);
            double s_y = restart_hull_displacement(1);
            double s_z = restart_hull_displacement(2);
            double v_x = restart_hull_quat_vector(0);
            double  v_y = restart_hull_quat_vector(1);
            double v_z = restart_hull_quat_vector(2);
            double s = restart_hull_quat_scalar;

            Point<3> ref_baricenter_pos(comp_dom.boat_model.reference_hull_baricenter(0),
                                        comp_dom.boat_model.reference_hull_baricenter(1),
                                        comp_dom.boat_model.reference_hull_baricenter(2));
   
            Point<3> baricenter_pos(s_x+comp_dom.boat_model.reference_hull_baricenter(0),
                                    s_y+comp_dom.boat_model.reference_hull_baricenter(1),
                                    s_z+comp_dom.boat_model.reference_hull_baricenter(2));

            Point<3> RotMatRow1(1-2*v_y*v_y-2*v_z*v_z, 2*v_x*v_y-2*s*v_z, 2*v_x*v_z+2*s*v_y);
            Point<3> RotMatRow2(2*v_x*v_y+2*s*v_z, 1-2*v_x*v_x-2*v_z*v_z, 2*v_y*v_z-2*s*v_x);
            Point<3> RotMatRow3(2*v_x*v_z-2*s*v_y, 2*v_y*v_z+2*s*v_x, 1-2*v_y*v_y-2*v_x*v_x);

            Point<3> ref_point_pos(boat_mesh_point.X(),boat_mesh_point.Y(),boat_mesh_point.Z());
            Point<3> rigid_lin_displ(s_x,s_y,s_z);
            Point<3> target_point_pos(RotMatRow1*(ref_point_pos+(-1.0)*ref_baricenter_pos),
                                      RotMatRow2*(ref_point_pos+(-1.0)*ref_baricenter_pos),
                                      RotMatRow3*(ref_point_pos+(-1.0)*ref_baricenter_pos));
            target_point_pos += baricenter_pos;
            // now we have the point on the REFERENCE hull, and transform it to go onto the CURRENT hull
            //boat_mesh_point.Transform(reference_to_current_transformation);
            // the rigid motion map points is the difference between the reference point position and the
            // current (rigidly displaced) node position
            comp_dom.rigid_motion_map_points(3*i) = target_point_pos(0)-ref_point_pos(0);
            comp_dom.rigid_motion_map_points(3*i+1) = target_point_pos(1)-ref_point_pos(1);
            comp_dom.rigid_motion_map_points(3*i+2) = target_point_pos(2)-ref_point_pos(2);

            rigid_motion_velocities(3*i) = current_hull_velocity(0)+
                                           current_hull_ang_velocity(1)*(target_point_pos(2)-baricenter_pos(2))-
                                           current_hull_ang_velocity(2)*(target_point_pos(1)-baricenter_pos(1));
            rigid_motion_velocities(3*i+1) = current_hull_velocity(1)+
                                             current_hull_ang_velocity(2)*(target_point_pos(0)-baricenter_pos(0))-
                                             current_hull_ang_velocity(0)*(target_point_pos(2)-baricenter_pos(2));
            rigid_motion_velocities(3*i+2) = current_hull_velocity(2)+
                                             current_hull_ang_velocity(0)*(target_point_pos(1)-baricenter_pos(1))-
                                             current_hull_ang_velocity(1)*(target_point_pos(0)-baricenter_pos(0));
            

            //if (fabs(t-0.2) <1e-5)
            //cout<<"AFTER: "<<Pnt(boat_mesh_point)<<" vs "<<Pnt(original_boat_mesh_point)<<endl;
            //cout<<"RMMP: "<<comp_dom.rigid_motion_map_points(3*i)<<" "
            //              <<comp_dom.rigid_motion_map_points(3*i+1)<<" "
            //              <<comp_dom.rigid_motion_map_points(3*i+2)<<endl;
            //cout<<"NN: "<<comp_dom.rigid_motion_map_points(3*i)+nodes_positions(3*i)+comp_dom.ref_points[3*i](0)<<" "
            //            <<comp_dom.rigid_motion_map_points(3*i+1)+nodes_positions(3*i+1)+comp_dom.ref_points[3*i](1)<<" "
            //            <<comp_dom.rigid_motion_map_points(3*i+2)+nodes_positions(3*i+2)+comp_dom.ref_points[3*i](2)<<endl;
            }
         }

         
      for (unsigned int i=0; i<dh.n_dofs(); ++i) 
	  {
          solution(3*i+0) = all_out[4](i)-comp_dom.rigid_motion_map_points(3*i)-comp_dom.ref_points[3*i](0);
          solution(3*i+1) = all_out[5](i)-comp_dom.rigid_motion_map_points(3*i+1)-comp_dom.ref_points[3*i](1);
          solution(3*i+2) = all_out[6](i)-comp_dom.rigid_motion_map_points(3*i+2)-comp_dom.ref_points[3*i](2);
          solution_dot(3*i+0) = all_out[7](i)-rigid_motion_velocities(3*i);
          solution_dot(3*i+1) = all_out[8](i)-rigid_motion_velocities(3*i+1);
          solution_dot(3*i+2) = all_out[9](i)-rigid_motion_velocities(3*i+2);
	  solution(i+vector_dh.n_dofs()) = all_out[0](i);
          solution_dot(i+vector_dh.n_dofs()) = all_out[1](i);
          solution(i+vector_dh.n_dofs()+dh.n_dofs()) = all_out[2](i);
          solution_dot(i+vector_dh.n_dofs()+dh.n_dofs()) = all_out[3](i);
	  }


      for (unsigned int k=0; k<3; ++k)
          {
          solution(k+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = current_hull_velocity(k);
          solution_dot(k+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = current_hull_velocity_dot(k);
          solution(k+3+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = current_hull_displacement(k);
          solution_dot(k+3+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = current_hull_displacement_dot(k);
          solution(k+6+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = current_hull_ang_velocity(k);
          solution_dot(k+6+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = current_hull_ang_velocity_dot(k);
          solution(k+9+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = current_hull_quat_vector(k);
          solution_dot(k+9+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = current_hull_quat_vector_dot(k);
          }
      solution(12+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = current_hull_quat_scalar;
      solution_dot(12+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = current_hull_quat_scalar_dot;


      
      for (unsigned int i=0; i<vector_dh.n_dofs(); ++i) 
	  {
	  comp_dom.map_points(i) = solution(i)+comp_dom.rigid_motion_map_points(i);
	  }

   

      DXDt_and_DphiDt_vector.reinit(vector_dh.n_dofs()+dh.n_dofs());

      DphiDt_sparsity_pattern.reinit (comp_dom.dh.n_dofs(),
				      comp_dom.dh.n_dofs(),
				      comp_dom.dh.max_couplings_between_dofs());
      vector_sparsity_pattern.reinit (comp_dom.vector_dh.n_dofs(),
					comp_dom.vector_dh.n_dofs(),
					comp_dom.vector_dh.max_couplings_between_dofs());

      DoFTools::make_sparsity_pattern (comp_dom.dh, DphiDt_sparsity_pattern, constraints);
      DphiDt_sparsity_pattern.compress();

      DoFTools::make_sparsity_pattern (comp_dom.vector_dh, vector_sparsity_pattern, vector_constraints);
      vector_sparsity_pattern.compress();

      working_map_points.reinit (comp_dom.vector_dh.n_dofs());
      working_nodes_velocities.reinit(comp_dom.vector_dh.n_dofs());
      nodes_pos_res.reinit(comp_dom.vector_dh.n_dofs());
      nodes_ref_surf_dist.reinit(comp_dom.vector_dh.n_dofs());
      nodes_diff_jac_x_delta.reinit(dofs_number);
      nodes_alg_jac_x_delta.reinit(dofs_number);
      dae_nonlin_residual.reinit(dofs_number);
      dae_linear_step_residual.reinit(dofs_number); 
      current_sol.reinit(dofs_number);
      current_sol_dot.reinit(dofs_number);
      bem_residual.reinit (comp_dom.dh.n_dofs());
      bem_phi.reinit(comp_dom.dh.n_dofs());
      bem_dphi_dn.reinit(comp_dom.dh.n_dofs());
      temp_src.reinit(comp_dom.vector_dh.n_dofs());
      break_wave_press.reinit(comp_dom.dh.n_dofs());

  
      

      bem.reinit();
      
      support_points.resize(dh.n_dofs());
      DoFTools::map_dofs_to_support_points<2, 3>( *comp_dom.mapping, dh, support_points);
      std::vector<bool> new_boundary_dofs(vector_dh.n_dofs());
      std::vector< bool > comp_sel(3, true);
      DoFTools::extract_boundary_dofs(vector_dh, comp_sel, new_boundary_dofs);
      for (unsigned int i=0; i<dh.n_dofs();++i)
          for (unsigned int j=0; j<old_points.size();++j)
              if (old_points[j].distance(support_points[i]) < 1e-5)
                 {
                 new_boundary_dofs[3*i] = 1;
                 new_boundary_dofs[3*i+1] = 1;
                 new_boundary_dofs[3*i+2] = 1;
                 }


      for (unsigned int i=0; i<comp_dom.vector_dh.n_dofs(); ++i)
          {
          comp_dom.map_points(i) = solution(i)+comp_dom.rigid_motion_map_points(i);
          //cout<<"* "<<i<<" "<<comp_dom.map_points(i)<<endl;
          }


      make_edges_conformal(solution, solution_dot, t, step_number, h);
      make_edges_conformal(solution, solution_dot, t, step_number, h);
      remove_transom_hanging_nodes(solution, solution_dot, t, step_number, h);
      make_edges_conformal(solution, solution_dot, t, step_number, h);
      make_edges_conformal(solution, solution_dot, t, step_number, h);



      for (unsigned int k=0; k<3; ++k)
          {
          solution(k+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = current_hull_velocity(k);
          solution_dot(k+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = current_hull_velocity_dot(k);
          solution(k+3+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = current_hull_displacement(k);
          solution_dot(k+3+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = current_hull_displacement_dot(k);
          solution(k+6+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = current_hull_ang_velocity(k);
          solution_dot(k+6+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = current_hull_ang_velocity_dot(k);
          solution(k+9+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = current_hull_quat_vector(k);
          solution_dot(k+9+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = current_hull_quat_vector_dot(k);
          }
      solution(12+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = current_hull_quat_scalar;
      solution_dot(12+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = current_hull_quat_scalar_dot;

//cout<<"First save "<<endl;
//      std::string filename2 = ( "postRemesh1.vtu" );
//      output_results(filename2, t, solution, solution_dot);
      if (!comp_dom.no_boat)
         comp_dom.evaluate_ref_surf_distances(nodes_ref_surf_dist,false);
      comp_dom.map_points -= nodes_ref_surf_dist;
      comp_dom.update_support_points();
      for (unsigned int i=0; i<comp_dom.vector_dh.n_dofs(); ++i)
          {
          solution(i) = comp_dom.map_points(i)-comp_dom.rigid_motion_map_points(i);
          }



//cout<<"Second save "<<endl;
      //std::string filename20 = ( "postRemesh20.vtu" );
      //output_results(filename20, t, solution, solution_dot);

      // in particular we must get the position of the nodes (in terms of curvilinear length)
      // on the smoothing lines, and the corresponding potential and horizontal velcoity values, in order to
      // interpolate the new values to be assigned at the restart of the simulation

      for (unsigned smooth_id=0; smooth_id<comp_dom.line_smoothers.size(); ++smooth_id)
          {
          Vector<double>  &old_lengths = comp_dom.line_smoothers[smooth_id]->get_lengths_before_smoothing();
          Vector<double>  &new_lengths = comp_dom.line_smoothers[smooth_id]->get_lengths_after_smoothing();
          std::vector<unsigned int> &indices = comp_dom.line_smoothers[smooth_id]->get_node_indices();
          Vector<double> old_potentials(old_lengths.size());
          Vector<double> old_vx(old_lengths.size());
          Vector<double> old_vy(old_lengths.size());
          //Vector<double> new_potentials(old_lengths.size());
          for (unsigned int i=0; i<old_lengths.size();++i)
              {
              old_potentials(i) = solution(indices[i]+vector_dh.n_dofs());
              old_vx(i) = solution_dot(3*indices[i]);
              old_vy(i) = solution_dot(3*indices[i]+1);
              //cout<<i<<"("<<indices[i]<<"->"<<round(indices[i]/3)<<")    "<<old_lengths(i)<<" vs "<<new_lengths(i)<<"  pot: "<<old_potentials(i)<<endl;
              //cout<<indices[i]<<" "<<comp_dom.support_points[indices[i]]<<endl;
              }
          //new_potentials(0) = old_potentials(0);
          //new_potentials(old_lengths.size()-1) = old_potentials(old_lengths.size()-1);
          for (unsigned int i=1; i<old_lengths.size()-1;++i)
              {
              unsigned int jj=1000000;
              for (unsigned int j=1; j<old_lengths.size();++j)
                  {
                  if (new_lengths(i) < old_lengths(j))
                     {
                     jj = j;
                     break;
                     }
                  }
              double fraction = (new_lengths(i)-old_lengths(jj-1))/(old_lengths(jj)-old_lengths(jj-1));
              solution(indices[i]+vector_dh.n_dofs()) = old_potentials(jj-1)+(old_potentials(jj)-old_potentials(jj-1))*fraction;
              //solution_dot(3*indices[i]) = old_vx(jj-1)+(old_vx(jj)-old_vx(jj-1))*fraction;
              //solution_dot(3*indices[i]+1) = old_vy(jj-1)+(old_vy(jj)-old_vy(jj-1))*fraction;
              //cout<<i<<" ---> "<<jj<<" "<<fraction<<" "<<old_potentials(jj-1)<<" "<<old_potentials(jj)<<" "<<new_potentials(i)<<endl;
              }
          }

      for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
          {
          if ( (!(comp_dom.flags[i] & near_boat)==0) &&
               (!(comp_dom.flags[i] & near_water)==0))
             {
             solution_dot(3*i) = 0;
             solution_dot(3*i+1) = 0;
             }
          }



  if (!comp_dom.no_boat && comp_dom.boat_model.is_transom)
     {
     comp_dom.update_support_points();
     double transom_draft = fabs(comp_dom.boat_model.CurrentPointCenterTransom(2));
     //double transom_draft = ref_transom_wet_surface/(fabs(comp_dom.boat_model.PointLeftTransom(1))+
     //                                                fabs(comp_dom.boat_model.PointRightTransom(1)));
     double transom_aspect_ratio = (fabs(comp_dom.boat_model.CurrentPointLeftTransom(1))+
                                     fabs(comp_dom.boat_model.CurrentPointRightTransom(1)))/transom_draft;

     wind.set_time(100000000.0);
     Vector<double> instantWindValueTinf(dim);
     Point<dim> zero(0,0,0);
     wind.vector_value(zero,instantWindValueTinf);
     Point<dim> VinfTinf;
     for (unsigned int i = 0; i < dim; i++)
       VinfTinf(i) = instantWindValueTinf(i);
     wind.set_time(initial_time);

     double FrT = sqrt(VinfTinf*VinfTinf)/sqrt(9.81*transom_draft);
     double ReT = sqrt(9.81*pow(transom_draft,3.0))/1.307e-6;
     double eta_dry = fmin(0.05*pow(FrT,2.834)*pow(transom_aspect_ratio,0.1352)*pow(ReT,0.01338),1.0);
     double lh = 0.0;
     //if (eta_dry < 1.0)
        lh = 5.0; 
        //lh = 0.1135*pow(FrT,3.025)*pow(transom_aspect_ratio,0.4603)*pow(ReT,-0.1514);
        //lh = 0.3265*pow(FrT,3.0) - 1.7216*pow(FrT,2.0) + 2.7593*FrT;
     
     cout<<"****eta_dry: "<<eta_dry<<endl;


     for(unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
        {
        comp_dom.initial_map_points(3*i+2) = 0.0;
        }

     for(unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
        {
        if ( (comp_dom.support_points[i](1) < comp_dom.boat_model.PointRightTransom(1)) &&
             (comp_dom.support_points[i](1) >= 0.0) &&
             (comp_dom.support_points[i](0) > comp_dom.boat_model.PointCenterTransom(0)-fabs(comp_dom.boat_model.PointCenterTransom(2)) ) &&
             (comp_dom.support_points[i](0) < comp_dom.boat_model.PointCenterTransom(0)+lh*fabs(comp_dom.boat_model.PointCenterTransom(2)) )   &&
             (comp_dom.support_points[i](2) > 2.0*comp_dom.boat_model.CurrentPointCenterTransom(2))    )
             {
             Point <3> dP0 = comp_dom.support_points[i];
             Point <3> dP; 
         				   //this is the vertical plane
             Handle(Geom_Plane) vertPlane = new Geom_Plane(0.,1.,0.,-dP0(1));
             vertPlane->Transform(comp_dom.boat_model.current_loc.Inverted());
             Handle(Geom_Curve) curve = comp_dom.boat_model.right_transom_bspline;

             GeomAPI_IntCS Intersector(curve, vertPlane);
             int npoints = Intersector.NbPoints();
             AssertThrow((npoints != 0), ExcMessage("Transom curve is not intersecting with vertical plane!"));
             //cout<<"Number of intersections: "<<npoints<<endl;
             double minDistance=1e7;
             for (int j=0; j<npoints;++j)
                 {
                 gp_Pnt int_point = Intersector.Point(j+1);
                 int_point.Transform(comp_dom.boat_model.current_loc);
                 Point<3> inters = Pnt(int_point);

                 if (dP0.distance(inters) < minDistance)
                    {
                    minDistance = dP0.distance(inters);
                    dP = inters;
                    }
                 }

             if ( (comp_dom.ref_points[i](0) > dP(0)+comp_dom.min_diameter/20.0) &&
                  (comp_dom.ref_points[i](0) < dP(0)+lh*fabs(dP(2))+comp_dom.min_diameter/20.0) )
                {
                double mean_curvature;
                Point<3> normal;
                Point<3> projection;
                comp_dom.boat_model.boat_surface_right->normal_projection_and_diff_forms(projection,
                                                                                        normal,
                                                                                        mean_curvature,
			                                                                dP);
                
                AssertThrow((dP.distance(projection) < 1e-4*comp_dom.boat_model.boatWetLength), ExcMessage("Normal projection for surface normal evaluation went wrong!"));
                double transom_slope = -normal(0)/normal(2);
                double a = -transom_slope/(lh*fabs(dP(2))) - 1/dP(2)/lh/lh;
                double x = comp_dom.support_points[i](0)-dP(0);
                comp_dom.initial_map_points(3*i+2) = a*x*x + transom_slope*x + dP(2);
                //cout<<"a="<<a<<"  t_s="<<transom_slope<<"  dP(2)="<<dP(2)<<endl;
                //cout<<"x="<<x<<"  z_in="<<comp_dom.initial_map_points(3*i+2)<<endl;
                }
             }
             else if ( (comp_dom.support_points[i](1) > comp_dom.boat_model.PointLeftTransom(1)) &&
                       (comp_dom.support_points[i](1) < 0.0) &&
                       (comp_dom.support_points[i](0) > comp_dom.boat_model.PointCenterTransom(0)-fabs(comp_dom.boat_model.PointCenterTransom(2))) &&
                       (comp_dom.support_points[i](0) < comp_dom.boat_model.PointCenterTransom(0)+lh*fabs(comp_dom.boat_model.PointCenterTransom(2)))  &&
                       (comp_dom.support_points[i](2) > 2.0*comp_dom.boat_model.CurrentPointCenterTransom(2))   )
             {
             Point <3> dP0 = comp_dom.support_points[i];
             Point <3> dP; 
         				   //this is the vertical plane
             Handle(Geom_Plane) vertPlane = new Geom_Plane(0.,1.,0.,-dP0(1));
             vertPlane->Transform(comp_dom.boat_model.current_loc.Inverted());
             Handle(Geom_Curve) curve = comp_dom.boat_model.left_transom_bspline;

             GeomAPI_IntCS Intersector(curve, vertPlane);
             int npoints = Intersector.NbPoints();
             AssertThrow((npoints != 0), ExcMessage("Transom curve is not intersecting with vertical plane!"));
             //cout<<"Number of intersections: "<<npoints<<endl;
             double minDistance=1e7;
             for (int j=0; j<npoints;++j)
                 {
                 gp_Pnt int_point = Intersector.Point(j+1);
                 int_point.Transform(comp_dom.boat_model.current_loc);
                 Point<3> inters = Pnt(int_point);

                 if (dP0.distance(inters) < minDistance)
                    {
                    minDistance = dP0.distance(inters);
                    dP = inters;
                    }
                 }
             if ( (comp_dom.ref_points[i](0) > dP(0)+comp_dom.min_diameter/20.0) &&
                  (comp_dom.ref_points[i](0) < dP(0)+lh*fabs(dP(2))+comp_dom.min_diameter/20.0) )
                {
                double mean_curvature;
                Point<3> normal;
                Point<3> projection;
                comp_dom.boat_model.boat_surface_left->normal_projection_and_diff_forms(projection,
                                                                                        normal,
                                                                                        mean_curvature,
			                                                                dP);
                
                AssertThrow((dP.distance(projection) < 1e-4*comp_dom.boat_model.boatWetLength), ExcMessage("Normal projection for surface normal evaluation went wrong!"));
                double transom_slope = -normal(0)/normal(2);
                double a = -transom_slope/(lh*fabs(dP(2))) - 1/dP(2)/lh/lh;
                double x = comp_dom.support_points[i](0)-dP(0);
                comp_dom.initial_map_points(3*i+2) = a*x*x + transom_slope*x + dP(2);
                }
             }
        }

  comp_dom.vector_constraints.distribute(comp_dom.initial_map_points);

  }

      
      for (unsigned int i=0; i<comp_dom.vector_dh.n_dofs(); ++i)
          {
          solution(i) = comp_dom.map_points(i)-comp_dom.rigid_motion_map_points(i);
          } 



      wind.set_time(t);  
      Vector<double> instantWindValue(dim);
      Point<dim> zero(0,0,0);
      wind.vector_value(zero,instantWindValue);
      Point<dim> Vinf;
      for (unsigned int i = 0; i < dim; i++)
          Vinf(i) = instantWindValue(i);
      Point<dim> Vh(0.0,0.0,0.01*2*dealii::numbers::PI*cos(2*dealii::numbers::PI*t));
      
      for (unsigned int i=0; i<dh.n_dofs(); ++i)
	  {
	  if ( comp_dom.flags[i] & boat  )
             //solution(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()) = -comp_dom.iges_normals[i]*Vinf; 
             solution(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()) = comp_dom.node_normals[i]*(Vh-Vinf);
	  }

      differential_components();

      current_time = t;
      current_sol = solution;
      current_sol_dot = solution_dot;
      last_remesh_time = t;
      comp_dom.old_map_points = comp_dom.map_points;
      comp_dom.rigid_motion_map_points = 0;


      prepare_restart(t, solution, solution_dot,true);
      //std::string filename3 = ( "postPostRemesh.vtu" );
      //output_results(filename3, t, solution, solution_dot);
//*/

     // we have to compute the reference transom wet surface with the new mesh here
   if (!comp_dom.no_boat && comp_dom.boat_model.is_transom)
      {
        double transom_wet_surface = 0;
      std::vector<Point<3> > vertices;
      std::vector<CellData<2> > cells;
      for (tria_it elem=comp_dom.tria.begin_active(); elem!= comp_dom.tria.end();++elem)
          {
	  if ((elem->material_id() == comp_dom.wall_sur_ID1 ||
	       elem->material_id() == comp_dom.wall_sur_ID2 ||
	       elem->material_id() == comp_dom.wall_sur_ID3 )) 
	     {
	     if (elem->at_boundary())
                {
                for (unsigned int f=0; f<GeometryInfo<2>::faces_per_cell; ++f)
		    if ( elem->face(f)->boundary_indicator() == 32 ||
                         elem->face(f)->boundary_indicator() == 37 )
                       {
                       unsigned int index_0 = comp_dom.find_point_id(elem->face(f)->vertex(0),comp_dom.ref_points);
                       unsigned int index_1 = comp_dom.find_point_id(elem->face(f)->vertex(1),comp_dom.ref_points);
                       if (comp_dom.ref_points[3*index_1](1) < comp_dom.ref_points[3*index_0](1))
                          {
                          vertices.push_back(comp_dom.ref_points[3*index_0]);
                          vertices.push_back(comp_dom.ref_points[3*index_1]);
                          vertices.push_back(comp_dom.ref_points[3*index_1]+
                                             -1.0*Point<3>(0.0,0.0,comp_dom.ref_points[3*index_1](2)));
                          vertices.push_back(comp_dom.support_points[index_0]+
                                             -1.0*Point<3>(0.0,0.0,comp_dom.ref_points[3*index_0](2)));
                          cells.resize(cells.size()+1);
                          cells[cells.size()-1].vertices[0]=4*(cells.size()-1)+0;
                          cells[cells.size()-1].vertices[1]=4*(cells.size()-1)+1;
                          cells[cells.size()-1].vertices[2]=4*(cells.size()-1)+2;
                          cells[cells.size()-1].vertices[3]=4*(cells.size()-1)+3;
                          }
                       else
                          {
                          vertices.push_back(comp_dom.support_points[index_1]);
                          vertices.push_back(comp_dom.support_points[index_0]);
                          vertices.push_back(comp_dom.support_points[index_0]+
                                             -1.0*Point<3>(0.0,0.0,comp_dom.ref_points[3*index_0](2)));
                          vertices.push_back(comp_dom.support_points[index_1]+
                                             -1.0*Point<3>(0.0,0.0,comp_dom.ref_points[3*index_1](2)));
                          cells.resize(cells.size()+1);
                          cells[cells.size()-1].vertices[0]=4*(cells.size()-1)+0;
                          cells[cells.size()-1].vertices[1]=4*(cells.size()-1)+1;
                          cells[cells.size()-1].vertices[2]=4*(cells.size()-1)+2;
                          cells[cells.size()-1].vertices[3]=4*(cells.size()-1)+3;
                          }                          
                       }
	        }
             }
          }
      

      SubCellData subcelldata;
      Triangulation<dim-1, dim> transom_tria;
      GridTools::delete_unused_vertices (vertices, cells, subcelldata);
      GridReordering<2,3>::reorder_cells (cells);
      transom_tria.create_triangulation_compatibility(vertices, cells, subcelldata );
      FE_Q<dim-1,dim> transom_fe(1);
      DoFHandler<dim-1,dim> transom_dh(transom_tria);
      transom_dh.distribute_dofs(transom_fe);

      FEValues<dim-1,dim> transom_fe_v(transom_fe, *comp_dom.quadrature,
			               update_values | update_gradients |
			               update_cell_normal_vectors |
			               update_quadrature_points |
			               update_JxW_values);

      const unsigned int transom_n_q_points = transom_fe_v.n_quadrature_points;
   

      std::vector<double> transom_pressure_quad_values(transom_n_q_points);
      for (cell_it cell = transom_dh.begin_active(); cell!=transom_dh.end(); ++cell)
          {
          transom_fe_v.reinit(cell);
     
          for (unsigned int q=0; q<transom_n_q_points; ++q)
              {
              transom_wet_surface += 1 * transom_fe_v.JxW(q);
              }
          }
      ref_transom_wet_surface = transom_wet_surface;
      }

      // here we finally compute the reference cell areas, which are needed
      // to evaluate boat cell deformations and possibly activate smoothing

      const unsigned int DphiDt_n_q_points = fe_v.n_quadrature_points;
      ref_cell_areas.reinit(comp_dom.tria.n_active_cells());
      cell = comp_dom.dh.begin_active(),
      endc = comp_dom.dh.end();

      unsigned int count=0;
      for (; cell!=endc; ++cell)
          {
          if ((cell->material_id() == comp_dom.wall_sur_ID1 ||
               cell->material_id() == comp_dom.wall_sur_ID2 ||
               cell->material_id() == comp_dom.wall_sur_ID3 ))
             {
             fe_v.reinit(cell);
             for (unsigned int q=0; q<DphiDt_n_q_points; ++q)
                 {
                 ref_cell_areas(count) += 1.0 * fe_v.JxW(q);
                 }
             }
          ++count;
          }





      std::cout<<"...done checking and updating mesh"<<std::endl;


      return true;

    }

   // this is what happens if max number of nodes has already been reached: nothing
   else if( (t>=remeshing_period*remeshing_counter && comp_dom.dh.n_dofs() > max_number_of_dofs) )
    {
    current_time = t;
    current_sol = solution;
    current_sol_dot = solution_dot;
    return false;
    }
   // in normal situations (no remesh) we just check that the boat mesh is ok (no bad quality cells)
   // and adjust it in case quality is bad. IMPORTANT: the restart procedure has to be improved
   // for this case, bacause the blend_factor used to obtain the last available horizontal displacement
   // field on the free surface is not the one used at restart time. as the difference is very small (1e-4~1e-5)
   // it doesn't seem to harm the restart, but with faster dynamics and ramps restarts might be problematic. thus,
   // we might consider adding nodes horizontal coordinates dofs to algebraic nonlinear restart problem class.
   else
    {
    FEValues<dim-1,dim> fe_v(*comp_dom.mapping, comp_dom.fe, *comp_dom.quadrature,
			     update_JxW_values);
    const unsigned int DphiDt_n_q_points = fe_v.n_quadrature_points;
    Vector<double> cell_areas(comp_dom.tria.n_active_cells());
    cell_areas=0;
    cell_it
    cell = comp_dom.dh.begin_active(),
    endc = comp_dom.dh.end();

    unsigned int count=0;
    bool smooth=false;
    double min_ratio =1.0;
    for (; cell!=endc; ++cell)
        {
        if ((cell->material_id() == comp_dom.wall_sur_ID1 ||
             cell->material_id() == comp_dom.wall_sur_ID2 ||
             cell->material_id() == comp_dom.wall_sur_ID3 ))
           {
           fe_v.reinit(cell);
           for (unsigned int q=0; q<DphiDt_n_q_points; ++q)
               {
               cell_areas(count) += 1.0 * fe_v.JxW(q);
               }
           }
        min_ratio = fmin(min_ratio,cell_areas(count)/ref_cell_areas(count));
        if (cell_areas(count)/ref_cell_areas(count) < 0.5)
           {
           smooth=true;
           break;
           }
        ++count;
        }
    cout<<"Min Areas Ratio For Smoothing: "<<min_ratio<<endl;

    if (smooth)
       {
      // first thing to do is update the geometry to the current positon
      Point<3> current_hull_displacement;
      Point<3> current_hull_displacement_dot;
      Point<3> current_hull_velocity;
      Point<3> current_hull_velocity_dot;
      Point<3> current_hull_ang_velocity;
      Point<3> current_hull_ang_velocity_dot;
      Point<3> current_hull_quat_vector;
      Point<3> current_hull_quat_vector_dot;
      for (unsigned int k=0; k<3; ++k)
          {
          current_hull_displacement(k) = solution(k+3+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs());
          current_hull_displacement_dot(k) = solution_dot(k+3+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs());
          current_hull_velocity(k) = solution(k+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs());
          current_hull_velocity_dot(k) = solution_dot(k+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs());
          current_hull_ang_velocity(k) = solution(k+6+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs());
          current_hull_ang_velocity_dot(k) = solution_dot(k+6+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs());
          current_hull_quat_vector(k) = solution(k+9+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs());
          current_hull_quat_vector_dot(k) = solution_dot(k+9+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs());
          }
      double current_hull_quat_scalar = solution(12+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs());
      double current_hull_quat_scalar_dot = solution_dot(12+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs());

      restart_hull_displacement(0) = current_hull_displacement(0);
      restart_hull_displacement(1) = current_hull_displacement(1);
      restart_hull_displacement(2) = current_hull_displacement(2);
      restart_hull_quat_vector(0) = current_hull_quat_vector(0);
      restart_hull_quat_vector(1) = current_hull_quat_vector(1);
      restart_hull_quat_vector(2) = current_hull_quat_vector(2);
      restart_hull_quat_scalar = current_hull_quat_scalar;

       // this moves the rigid_motion_map_points vector (target positions) for the internal nodes of the hull 
       // mesh (all nodes except for the free surface ones)
       Vector<double> rigid_motion_velocities(comp_dom.vector_dh.n_dofs());
      restart_hull_location = comp_dom.boat_model.set_current_position(restart_hull_displacement,
                                                                       restart_hull_quat_scalar,
                                                                       restart_hull_quat_vector);



      //std::string filename4 = ( "beforeSurfaceRemesh.vtu" );
      //output_results(filename4, t, solution, solution_dot);
       nodes_ref_surf_dist = 0.0;
       if (!comp_dom.no_boat)
          comp_dom.evaluate_ref_surf_distances(nodes_ref_surf_dist,true);
       differential_components();
       for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
           if ( ((comp_dom.flags[i] & boat) &&
                !(comp_dom.flags[i] & near_water) ) ||
                (comp_dom.flags[i] & transom_on_water)  )
              {
              for (unsigned int d=0; d<3; ++d)
                  {
                  comp_dom.map_points(3*i+d) -= nodes_ref_surf_dist(3*i+d);
                  comp_dom.old_map_points(3*i+d) = comp_dom.map_points(3*i+d);
                  }
              if (constraints.is_constrained(i))
                 for (unsigned int d=0; d<3; ++d)
                     solution(i) = comp_dom.map_points(i);
              else
                 {
                 gp_Pnt original_boat_mesh_point = Pnt(Point<3>(comp_dom.ref_points[3*i](0)+comp_dom.old_map_points(3*i),
                                                                comp_dom.ref_points[3*i](1)+comp_dom.old_map_points(3*i+1),
                                                                comp_dom.ref_points[3*i](2)+comp_dom.old_map_points(3*i+2)));
                 gp_Pnt boat_mesh_point = original_boat_mesh_point;
                 // we first take this point (which is in the RESTART hull location) and transform it to be in the
                 // REFERENCE configuration 
                 boat_mesh_point.Transform(restart_hull_location.Inverted());
                 Point<3> rigid_displ = Pnt(boat_mesh_point)+(-1.0*Pnt(original_boat_mesh_point));
                 for (unsigned int d=0; d<3; ++d)
                     comp_dom.rigid_motion_map_points(3*i+d) = rigid_displ(d);
                 for (unsigned int d=0; d<3; ++d)
                     solution(3*i+d) = comp_dom.map_points(3*i+d)-comp_dom.rigid_motion_map_points(3*i+d);
                 }
              //cout<<i<<" "<<nodes_ref_surf_dist(i)<<endl;
              }

      current_time = t;
      current_sol = solution;
      current_sol_dot = solution_dot;
      //std::string filename5 = ( "surfaceRemesh.vtu" );
      //output_results(filename5, t, solution, solution_dot);
       
       prepare_restart(t, solution, solution_dot,false);

       // compute new ref_cell_areas
       ref_cell_areas = 0.0;
       cell = comp_dom.dh.begin_active(),
       endc = comp_dom.dh.end();

       unsigned int count=0;
       for (; cell!=endc; ++cell)
           {
           if ((cell->material_id() == comp_dom.wall_sur_ID1 ||
                cell->material_id() == comp_dom.wall_sur_ID2 ||
                cell->material_id() == comp_dom.wall_sur_ID3 ))
              {
              fe_v.reinit(cell);
              for (unsigned int q=0; q<DphiDt_n_q_points; ++q)
                  {
                  ref_cell_areas(count) += 1.0 * fe_v.JxW(q);
                  }
              }
           ++count;
           }

       return true;
       }
    else
       {
       return false;
       }

    }



    return false;
}


                                      // this routine detects if mesh is not
                                      // conformal at edges (because of double
                                      // nodes) and makes the refinements needed
                                      // to make it conformal. in free surface
                                      // it also has to deal with interpolation
                                      // of the solution

template <int dim>
void FreeSurface<dim>::make_edges_conformal(Vector<double> & solution,
				            Vector<double> &solution_dot,
				            const double t,
	                                    const unsigned int step_number,
		                            const double  h)
{
std::cout<<"Restoring mesh conformity..."<<std::endl;



  Point<3> hull_lin_vel;
  Point<3> hull_lin_displ;
  Point<3> hull_lin_vel_dot;
  Point<3> hull_lin_displ_dot;
  Point<3> hull_ang_vel;
  Point<3> hull_ang_vel_dot;
  Point<3> hull_quat_vect;
  Point<3> hull_quat_vect_dot;
  double hull_quat_scal;
  double hull_quat_scal_dot;
  for (unsigned int d=0; d<3;++d)
      {
      hull_lin_vel(d) = solution(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+d);
      hull_lin_displ(d) = solution(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+3+d);
      hull_lin_vel_dot(d) = solution_dot(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+d);
      hull_lin_displ_dot(d) = solution_dot(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+3+d);
      hull_ang_vel(d) = solution(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+6+d);
      hull_ang_vel_dot(d) = solution_dot(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+6+d);
      hull_quat_vect(d) = solution(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+9+d);
      hull_quat_vect_dot(d) = solution_dot(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+9+d);
      }
  hull_quat_scal = solution(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+12);
  hull_quat_scal_dot = solution_dot(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+12);

  cout<<"*Hull Rigid Displacement: "<<hull_lin_displ<<endl;
  cout<<"*Hull Rigid Velocity: "<<hull_lin_vel<<endl;
  cout<<"*Hull Angualr Velocity: "<<hull_lin_vel<<endl;


     Vector<double> rigid_motion_velocities(comp_dom.vector_dh.n_dofs());

     // this moves the rigid_motion_map_points vector (target positions) for the internal nodes of the hull 
     // mesh (all nodes except for the free surface ones)
     for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
         {
         if ( ((comp_dom.flags[i] & boat) &&
              !(comp_dom.flags[i] & near_water) ) ||
              (comp_dom.flags[i] & transom_on_water) )
            {
            //if (fabs(t-0.2) <1e-5)
            //cout<<"BEFORE: "<<comp_dom.ref_points[3*i](0)+comp_dom.old_map_points(3*i)<<" "
            //                <<comp_dom.ref_points[3*i](1)+comp_dom.old_map_points(3*i+1)<<" "
            //                <<comp_dom.ref_points[3*i](2)+comp_dom.old_map_points(3*i+2)<<endl;
            gp_Pnt original_boat_mesh_point = Pnt(Point<3>(comp_dom.ref_points[3*i](0)+comp_dom.old_map_points(3*i),
                                                           comp_dom.ref_points[3*i](1)+comp_dom.old_map_points(3*i+1),
                                                           comp_dom.ref_points[3*i](2)+comp_dom.old_map_points(3*i+2)));
            gp_Pnt boat_mesh_point = original_boat_mesh_point;
            // we first take this point (which is in the RESTART hull location) and transform it to be in the
            // REFERENCE configuration 
            boat_mesh_point.Transform(restart_hull_location.Inverted());

            // now we use sacado to compute the residual at this dof, along with its derivatives with respect to 
            // the 7 dofs associated to the rigid linear and angular displacements

            double s_x = hull_lin_displ(0);
            double s_y = hull_lin_displ(1);
            double s_z = hull_lin_displ(2);
            double v_x = hull_quat_vect(0);
            double  v_y = hull_quat_vect(1);
            double v_z = hull_quat_vect(2);
            double s = hull_quat_scal;

            Point<3> ref_baricenter_pos(comp_dom.boat_model.reference_hull_baricenter(0),
                                        comp_dom.boat_model.reference_hull_baricenter(1),
                                        comp_dom.boat_model.reference_hull_baricenter(2));
   
            Point<3> baricenter_pos(s_x+comp_dom.boat_model.reference_hull_baricenter(0),
                                    s_y+comp_dom.boat_model.reference_hull_baricenter(1),
                                    s_z+comp_dom.boat_model.reference_hull_baricenter(2));

            Point<3> RotMatRow1(1-2*v_y*v_y-2*v_z*v_z, 2*v_x*v_y-2*s*v_z, 2*v_x*v_z+2*s*v_y);
            Point<3> RotMatRow2(2*v_x*v_y+2*s*v_z, 1-2*v_x*v_x-2*v_z*v_z, 2*v_y*v_z-2*s*v_x);
            Point<3> RotMatRow3(2*v_x*v_z-2*s*v_y, 2*v_y*v_z+2*s*v_x, 1-2*v_y*v_y-2*v_x*v_x);

            Point<3> ref_point_pos(boat_mesh_point.X(),boat_mesh_point.Y(),boat_mesh_point.Z());
            Point<3> rigid_lin_displ(s_x,s_y,s_z);
            Point<3> target_point_pos(RotMatRow1*(ref_point_pos+(-1.0)*ref_baricenter_pos),
                                      RotMatRow2*(ref_point_pos+(-1.0)*ref_baricenter_pos),
                                      RotMatRow3*(ref_point_pos+(-1.0)*ref_baricenter_pos));
            target_point_pos += baricenter_pos;
            // now we have the point on the REFERENCE hull, and transform it to go onto the CURRENT hull
            //boat_mesh_point.Transform(reference_to_current_transformation);
            // the rigid motion map points is the difference between the reference point position and the
            // current (rigidly displaced) node position
            comp_dom.rigid_motion_map_points(3*i) = target_point_pos(0)-ref_point_pos(0);
            comp_dom.rigid_motion_map_points(3*i+1) = target_point_pos(1)-ref_point_pos(1);
            comp_dom.rigid_motion_map_points(3*i+2) = target_point_pos(2)-ref_point_pos(2);

            rigid_motion_velocities(3*i) = hull_lin_vel(0)+
                                           hull_ang_vel(1)*(target_point_pos(2)-baricenter_pos(2))-
                                           hull_ang_vel(2)*(target_point_pos(1)-baricenter_pos(1));
            rigid_motion_velocities(3*i+1) = hull_lin_vel(1)+
                                             hull_ang_vel(2)*(target_point_pos(0)-baricenter_pos(0))-
                                             hull_ang_vel(0)*(target_point_pos(2)-baricenter_pos(2));
            rigid_motion_velocities(3*i+2) = hull_lin_vel(2)+
                                             hull_ang_vel(0)*(target_point_pos(1)-baricenter_pos(1))-
                                             hull_ang_vel(1)*(target_point_pos(0)-baricenter_pos(0));
            
            

            //if (fabs(t-0.2) <1e-5)
            //cout<<"AFTER: "<<Pnt(boat_mesh_point)<<" vs "<<Pnt(original_boat_mesh_point)<<endl;
            //cout<<"RMMP: "<<comp_dom.rigid_motion_map_points(3*i)<<" "
            //              <<comp_dom.rigid_motion_map_points(3*i+1)<<" "
            //              <<comp_dom.rigid_motion_map_points(3*i+2)<<endl;
            //cout<<"NN: "<<comp_dom.rigid_motion_map_points(3*i)+nodes_positions(3*i)+comp_dom.ref_points[3*i](0)<<" "
            //            <<comp_dom.rigid_motion_map_points(3*i+1)+nodes_positions(3*i+1)+comp_dom.ref_points[3*i](1)<<" "
            //            <<comp_dom.rigid_motion_map_points(3*i+2)+nodes_positions(3*i+2)+comp_dom.ref_points[3*i](2)<<endl;
            }
         }
  


      Triangulation<dim-1, dim> &tria = comp_dom.tria;
      DoFHandler<dim-1, dim> &dh = comp_dom.dh;
      DoFHandler<dim-1, dim> &vector_dh = comp_dom.vector_dh;
    

      std::vector<Point<dim> > support_points(dh.n_dofs());
      DoFTools::map_dofs_to_support_points<dim-1, dim>( *comp_dom.mapping, dh, support_points);



      VectorView<double> Phi(dh.n_dofs(), solution.begin()+vector_dh.n_dofs());
      VectorView<double> Phi_dot(dh.n_dofs(), solution_dot.begin()+vector_dh.n_dofs());
      VectorView<double> dphi_dn(dh.n_dofs(), solution.begin()+vector_dh.n_dofs()+dh.n_dofs());
      VectorView<double> dphi_dn_dot(dh.n_dofs(), solution_dot.begin()+vector_dh.n_dofs()+dh.n_dofs());

      std::vector<Point<3> > old_points(dh.n_dofs());
      old_points = support_points;

      Vector<double> curvatures(vector_dh.n_dofs());
      comp_dom.surface_smoother->compute_curvatures(curvatures);

				   // Get support points in the
				   // reference configuration
      std::vector<Point<3> > ref_points(dh.n_dofs());
      DoFTools::map_dofs_to_support_points<2,3>(StaticMappingQ1<2,3>::mapping,
					  dh, ref_points);
     double tol=1e-7;
     for (unsigned int i=0;i<dh.n_dofs();++i)
         {
         if ((comp_dom.flags[i] & edge) &&
             (comp_dom.double_nodes_set[i].size() == 1) )
            { 
            //we identify here the two vertices of the parent cell on the considered side
            //(which is the one with the non conformal node)
            std::vector<Point<3> > nodes(2);
            for (unsigned int kk=0; kk<2;++kk)
                {
                DoFHandler<2,3>::cell_iterator cell = comp_dom.dof_to_elems[i][kk];
                for (unsigned int f=0; f<GeometryInfo<2>::faces_per_cell; ++f)
                    {
	            if (cell->face(f)->at_boundary())
                       {
                       if (ref_points[i].distance(cell->face(f)->vertex(0)) <tol)
                          nodes[kk] = cell->face(f)->vertex(1);
                       else if (ref_points[i].distance(cell->face(f)->vertex(1)) <tol)
                          nodes[kk] = cell->face(f)->vertex(0);
                       }
                    }
                }
            // we can now compute the center of the parent cell face
            Point<3> parent_face_center = 0.5*(nodes[0]+nodes[1]);
            // now we look for the opposite side cell that has a face on an edge, having the same center
            DoFHandler<2,3>::cell_iterator cell1 = comp_dom.dof_to_elems[i][0]->parent();
            for(unsigned int f=0; f<GeometryInfo<2>::faces_per_cell; ++f)
	       if (cell1->face(f)->at_boundary())
                  {
	          for (typename std::set<tria_it>::iterator jt=comp_dom.edge_cells.begin(); jt != comp_dom.edge_cells.end(); ++jt)       
	               for (unsigned int d=0; d<GeometryInfo<2>::faces_per_cell; ++d)
		           if ((*jt)->face(d)->at_boundary())
		              if ( parent_face_center.distance((*jt)->face(d)->center()) < tol)
                                 {
                                 //(*jt)->set_refine_flag();
                                 if ((d==0) || (d==1))
                                       (*jt)->set_refine_flag(RefinementCase<2>::cut_axis(1));
                                    else
                                       (*jt)->set_refine_flag(RefinementCase<2>::cut_axis(0));
                                 }
                 }
            }           
         }

      SolutionTransfer<dim-1, Vector<double>, DoFHandler<dim-1, dim> > soltrans(dh);

      Vector<double> positions_x(comp_dom.dh.n_dofs());
      Vector<double> positions_y(comp_dom.dh.n_dofs());
      Vector<double> positions_z(comp_dom.dh.n_dofs());
      Vector<double> positions_dot_x(comp_dom.dh.n_dofs());
      Vector<double> positions_dot_y(comp_dom.dh.n_dofs());
      Vector<double> positions_dot_z(comp_dom.dh.n_dofs());

      for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
          {
          positions_x(i) = solution(3*i+0)+comp_dom.rigid_motion_map_points(3*i+0)+comp_dom.ref_points[3*i](0);
          positions_y(i) = solution(3*i+1)+comp_dom.rigid_motion_map_points(3*i+1)+comp_dom.ref_points[3*i](1);
          positions_z(i) = solution(3*i+2)+comp_dom.rigid_motion_map_points(3*i+2)+comp_dom.ref_points[3*i](2);
          positions_dot_x(i) = solution_dot(3*i+0);
          positions_dot_y(i) = solution_dot(3*i+1);
          positions_dot_z(i) = solution_dot(3*i+2);
          }          

      std::vector<Vector<double> > all_in;
      all_in.push_back((Vector<double>)Phi);
      all_in.push_back((Vector<double>)Phi_dot);
      all_in.push_back((Vector<double>)dphi_dn);
      all_in.push_back((Vector<double>)dphi_dn_dot);
      all_in.push_back(positions_x);
      all_in.push_back(positions_y);
      all_in.push_back(positions_z);
      all_in.push_back(positions_dot_x);
      all_in.push_back(positions_dot_y);
      all_in.push_back(positions_dot_z); 

      tria.prepare_coarsening_and_refinement();
      soltrans.prepare_for_coarsening_and_refinement(all_in);



      //std::cout << "Refined counter: " << refinedCellCounter << std::endl;
      tria.execute_coarsening_and_refinement();
      dh.distribute_dofs(comp_dom.fe);
      vector_dh.distribute_dofs(comp_dom.vector_fe);  
      comp_dom.map_points.reinit(vector_dh.n_dofs());
      comp_dom.smoothing_map_points.reinit(vector_dh.n_dofs());
      comp_dom.old_map_points.reinit(vector_dh.n_dofs());
      comp_dom.rigid_motion_map_points.reinit(vector_dh.n_dofs());
      comp_dom.initial_map_points.reinit(vector_dh.n_dofs());
      comp_dom.ref_points.resize(vector_dh.n_dofs());
      DoFTools::map_dofs_to_support_points<2,3>(StaticMappingQ1<2,3>::mapping,
					      comp_dom.vector_dh, comp_dom.ref_points);
      comp_dom.generate_double_nodes_set();

      compute_constraints(constraints, vector_constraints);

      dofs_number = vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()+13;

      std::cout<<"Total number of dofs after restoring edges conformity: "<<dh.n_dofs()<<std::endl;


      Vector<double> new_Phi(dh.n_dofs());
      Vector<double> new_Phi_dot(dh.n_dofs());

      Vector <double> new_dphi_dn(dh.n_dofs());     
      Vector <double> new_dphi_dn_dot(dh.n_dofs());
   
      Vector<double> new_positions_x(dh.n_dofs());
      Vector<double> new_positions_y(dh.n_dofs());
      Vector<double> new_positions_z(dh.n_dofs());
      Vector<double> new_positions_dot_x(dh.n_dofs());
      Vector<double> new_positions_dot_y(dh.n_dofs());
      Vector<double> new_positions_dot_z(dh.n_dofs());

      std::vector<Vector<double> > all_out;
      all_out.push_back(new_Phi);
      all_out.push_back(new_Phi_dot);
      all_out.push_back(new_dphi_dn);
      all_out.push_back(new_dphi_dn_dot);
      all_out.push_back(new_positions_x);
      all_out.push_back(new_positions_y);
      all_out.push_back(new_positions_z);
      all_out.push_back(new_positions_dot_x);
      all_out.push_back(new_positions_dot_y);
      all_out.push_back(new_positions_dot_z);

      soltrans.interpolate(all_in, all_out);
  
      solution.reinit(dofs_number);
      solution_dot.reinit(dofs_number); 


      constraints.distribute(all_out[0]);
      constraints.distribute(all_out[1]);
      constraints.distribute(all_out[2]);
      constraints.distribute(all_out[3]);      
      constraints.distribute(all_out[4]);
      constraints.distribute(all_out[5]);
      constraints.distribute(all_out[6]);
      constraints.distribute(all_out[7]);    
      constraints.distribute(all_out[8]);
      constraints.distribute(all_out[9]);
  

     // we have to compute rigid_motion_map_points and rigid_motion_velocities on the new mesh

     rigid_motion_velocities.reinit(comp_dom.vector_dh.n_dofs());


     // this moves the rigid_motion_map_points vector (target positions) for the internal nodes of the hull 
     // mesh (all nodes except for the free surface ones)
     for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
         {
         if ( ((comp_dom.flags[i] & boat) &&
              !(comp_dom.flags[i] & near_water) ) ||
              (comp_dom.flags[i] & transom_on_water) )
            {
            //if (fabs(t-0.2) <1e-5)
            //cout<<"BEFORE: "<<comp_dom.ref_points[3*i](0)+comp_dom.old_map_points(3*i)<<" "
            //                <<comp_dom.ref_points[3*i](1)+comp_dom.old_map_points(3*i+1)<<" "
            //                <<comp_dom.ref_points[3*i](2)+comp_dom.old_map_points(3*i+2)<<endl;
            gp_Pnt original_boat_mesh_point = Pnt(Point<3>(comp_dom.ref_points[3*i](0)+comp_dom.old_map_points(3*i),
                                                           comp_dom.ref_points[3*i](1)+comp_dom.old_map_points(3*i+1),
                                                           comp_dom.ref_points[3*i](2)+comp_dom.old_map_points(3*i+2)));
            gp_Pnt boat_mesh_point = original_boat_mesh_point;
            // we first take this point (which is in the RESTART hull location) and transform it to be in the
            // REFERENCE configuration 
            boat_mesh_point.Transform(restart_hull_location.Inverted());

            // now we use sacado to compute the residual at this dof, along with its derivatives with respect to 
            // the 7 dofs associated to the rigid linear and angular displacements

            double s_x = hull_lin_displ(0);
            double s_y = hull_lin_displ(1);
            double s_z = hull_lin_displ(2);
            double v_x = hull_quat_vect(0);
            double  v_y = hull_quat_vect(1);
            double v_z = hull_quat_vect(2);
            double s = hull_quat_scal;

            Point<3> ref_baricenter_pos(comp_dom.boat_model.reference_hull_baricenter(0),
                                        comp_dom.boat_model.reference_hull_baricenter(1),
                                        comp_dom.boat_model.reference_hull_baricenter(2));
   
            Point<3> baricenter_pos(s_x+comp_dom.boat_model.reference_hull_baricenter(0),
                                    s_y+comp_dom.boat_model.reference_hull_baricenter(1),
                                    s_z+comp_dom.boat_model.reference_hull_baricenter(2));

            Point<3> RotMatRow1(1-2*v_y*v_y-2*v_z*v_z, 2*v_x*v_y-2*s*v_z, 2*v_x*v_z+2*s*v_y);
            Point<3> RotMatRow2(2*v_x*v_y+2*s*v_z, 1-2*v_x*v_x-2*v_z*v_z, 2*v_y*v_z-2*s*v_x);
            Point<3> RotMatRow3(2*v_x*v_z-2*s*v_y, 2*v_y*v_z+2*s*v_x, 1-2*v_y*v_y-2*v_x*v_x);

            Point<3> ref_point_pos(boat_mesh_point.X(),boat_mesh_point.Y(),boat_mesh_point.Z());
            Point<3> rigid_lin_displ(s_x,s_y,s_z);
            Point<3> target_point_pos(RotMatRow1*(ref_point_pos+(-1.0)*ref_baricenter_pos),
                                      RotMatRow2*(ref_point_pos+(-1.0)*ref_baricenter_pos),
                                      RotMatRow3*(ref_point_pos+(-1.0)*ref_baricenter_pos));
            target_point_pos += baricenter_pos;
            // now we have the point on the REFERENCE hull, and transform it to go onto the CURRENT hull
            //boat_mesh_point.Transform(reference_to_current_transformation);
            // the rigid motion map points is the difference between the reference point position and the
            // current (rigidly displaced) node position
            comp_dom.rigid_motion_map_points(3*i) = target_point_pos(0)-ref_point_pos(0);
            comp_dom.rigid_motion_map_points(3*i+1) = target_point_pos(1)-ref_point_pos(1);
            comp_dom.rigid_motion_map_points(3*i+2) = target_point_pos(2)-ref_point_pos(2);

            rigid_motion_velocities(3*i) = hull_lin_vel(0)+
                                           hull_ang_vel(1)*(target_point_pos(2)-baricenter_pos(2))-
                                           hull_ang_vel(2)*(target_point_pos(1)-baricenter_pos(1));
            rigid_motion_velocities(3*i+1) = hull_lin_vel(1)+
                                             hull_ang_vel(2)*(target_point_pos(0)-baricenter_pos(0))-
                                             hull_ang_vel(0)*(target_point_pos(2)-baricenter_pos(2));
            rigid_motion_velocities(3*i+2) = hull_lin_vel(2)+
                                             hull_ang_vel(0)*(target_point_pos(1)-baricenter_pos(1))-
                                             hull_ang_vel(1)*(target_point_pos(0)-baricenter_pos(0));
            
            

            //if (fabs(t-0.2) <1e-5)
            //cout<<"AFTER: "<<Pnt(boat_mesh_point)<<" vs "<<Pnt(original_boat_mesh_point)<<endl;
            //cout<<"RMMP: "<<comp_dom.rigid_motion_map_points(3*i)<<" "
            //              <<comp_dom.rigid_motion_map_points(3*i+1)<<" "
            //              <<comp_dom.rigid_motion_map_points(3*i+2)<<endl;
            //cout<<"NN: "<<comp_dom.rigid_motion_map_points(3*i)+nodes_positions(3*i)+comp_dom.ref_points[3*i](0)<<" "
            //            <<comp_dom.rigid_motion_map_points(3*i+1)+nodes_positions(3*i+1)+comp_dom.ref_points[3*i](1)<<" "
            //            <<comp_dom.rigid_motion_map_points(3*i+2)+nodes_positions(3*i+2)+comp_dom.ref_points[3*i](2)<<endl;
            }
         }


      for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i) 
	  {
          solution(3*i+0) = all_out[4](i)-comp_dom.rigid_motion_map_points(3*i)-comp_dom.ref_points[3*i](0);
          solution(3*i+1) = all_out[5](i)-comp_dom.rigid_motion_map_points(3*i+1)-comp_dom.ref_points[3*i](1);
          solution(3*i+2) = all_out[6](i)-comp_dom.rigid_motion_map_points(3*i+2)-comp_dom.ref_points[3*i](2);
          solution_dot(3*i+0) = all_out[7](i)-rigid_motion_velocities(3*i);
          solution_dot(3*i+1) = all_out[8](i)-rigid_motion_velocities(3*i+1);
          solution_dot(3*i+2) = all_out[9](i)-rigid_motion_velocities(3*i+2);
	  solution(i+vector_dh.n_dofs()) = all_out[0](i);
          solution_dot(i+vector_dh.n_dofs()) = all_out[1](i);
          solution(i+vector_dh.n_dofs()+dh.n_dofs()) = all_out[2](i);
          solution_dot(i+vector_dh.n_dofs()+dh.n_dofs()) = all_out[3](i);
	  }

      for (unsigned int k=0; k<3; ++k)
          {
          solution(k+3+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = hull_lin_displ(k);
          solution_dot(k+3+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = hull_lin_displ_dot(k);
          solution(k+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = hull_lin_vel(k);
          solution_dot(k+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = hull_lin_vel_dot(k);
          solution(k+6+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = hull_ang_vel(k);
          solution_dot(k+6+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = hull_ang_vel_dot(k);
          solution(k+9+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = hull_quat_vect(k);
          solution_dot(k+9+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = hull_quat_vect_dot(k);          
          }
      solution(12+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = hull_quat_scal;
      solution_dot(12+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = hull_quat_scal_dot;  


      DXDt_and_DphiDt_vector.reinit(vector_dh.n_dofs()+dh.n_dofs());

      DphiDt_sparsity_pattern.reinit (dh.n_dofs(),
				      dh.n_dofs(),
				      dh.max_couplings_between_dofs());
      vector_sparsity_pattern.reinit (vector_dh.n_dofs(),
					vector_dh.n_dofs(),
					vector_dh.max_couplings_between_dofs());

      DoFTools::make_sparsity_pattern (dh, DphiDt_sparsity_pattern, constraints);
      DphiDt_sparsity_pattern.compress();

      DoFTools::make_sparsity_pattern (vector_dh, vector_sparsity_pattern, vector_constraints);
      vector_sparsity_pattern.compress();

      working_map_points.reinit(comp_dom.vector_dh.n_dofs());
      working_nodes_velocities.reinit(comp_dom.vector_dh.n_dofs());
      nodes_pos_res.reinit(comp_dom.vector_dh.n_dofs());
      nodes_ref_surf_dist.reinit(comp_dom.vector_dh.n_dofs());
      nodes_diff_jac_x_delta.reinit(dofs_number);
      nodes_alg_jac_x_delta.reinit(dofs_number);
      dae_nonlin_residual.reinit(dofs_number);
      dae_linear_step_residual.reinit(dofs_number);
      current_sol.reinit(dofs_number);
      current_sol_dot.reinit(dofs_number);
      bem_residual.reinit(comp_dom.dh.n_dofs()); 
      bem_phi.reinit(comp_dom.dh.n_dofs());
      bem_dphi_dn.reinit(comp_dom.dh.n_dofs());
      temp_src.reinit(comp_dom.vector_dh.n_dofs());
      break_wave_press.reinit(comp_dom.dh.n_dofs());
  
      

      bem.reinit();

      support_points.resize(dh.n_dofs());
      DoFTools::map_dofs_to_support_points<2, 3>( *comp_dom.mapping, dh, support_points);
      std::vector<bool> new_boundary_dofs(vector_dh.n_dofs());
      std::vector< bool > comp_sel(3, true);
      DoFTools::extract_boundary_dofs(vector_dh, comp_sel, new_boundary_dofs);
      for (unsigned int i=0; i<dh.n_dofs();++i)
          for (unsigned int j=0; j<old_points.size();++j)
              if (old_points[j].distance(support_points[i]) < 1e-5)
                 {
                 new_boundary_dofs[3*i] = 1;
                 new_boundary_dofs[3*i+1] = 1;
                 new_boundary_dofs[3*i+2] = 1;
                 }

      //comp_dom.apply_curvatures(new_curvatures,new_boundary_dofs);
      for (unsigned int i=0; i<comp_dom.vector_dh.n_dofs(); ++i)
          {
          comp_dom.map_points(i) = solution(i)+comp_dom.rigid_motion_map_points(i);
          }
      //std::string filename2 = ( "beforeCrash.vtu" );
      //output_results(filename2, t, solution, solution_dot);

      if (!comp_dom.no_boat)
         comp_dom.evaluate_ref_surf_distances(nodes_ref_surf_dist,false);
      comp_dom.map_points -= nodes_ref_surf_dist;
      for (unsigned int i=0; i<comp_dom.vector_dh.n_dofs(); ++i)
          {
          solution(i) = comp_dom.map_points(i)-comp_dom.rigid_motion_map_points(i);
          }


      // in particular we must get the position of the nodes (in terms of curvilinear length)
      // on the smoothing lines, and the corresponding potential and horizontal velcoity values, in order to
      // interpolate the new values to be assigned at the restart of the simulation

      for (unsigned smooth_id=0; smooth_id<comp_dom.line_smoothers.size(); ++smooth_id)
          {
          Vector<double>  &old_lengths = comp_dom.line_smoothers[smooth_id]->get_lengths_before_smoothing();
          Vector<double>  &new_lengths = comp_dom.line_smoothers[smooth_id]->get_lengths_after_smoothing();
          std::vector<unsigned int> &indices = comp_dom.line_smoothers[smooth_id]->get_node_indices();
          Vector<double> old_potentials(old_lengths.size());
          Vector<double> old_vx(old_lengths.size());
          Vector<double> old_vy(old_lengths.size());
          //Vector<double> new_potentials(old_lengths.size());
          for (unsigned int i=0; i<old_lengths.size();++i)
              {
              old_potentials(i) = solution(indices[i]+vector_dh.n_dofs());
              old_vx(i) = solution_dot(3*indices[i]);
              old_vy(i) = solution_dot(3*indices[i]+1);
              //cout<<i<<"("<<indices[i]<<"->"<<round(indices[i]/3)<<")    "<<old_lengths(i)<<" vs "<<new_lengths(i)<<"  pot: "<<old_potentials(i)<<endl;
              //cout<<indices[i]<<" "<<comp_dom.support_points[indices[i]]<<endl;
              }
          //new_potentials(0) = old_potentials(0);
          //new_potentials(old_lengths.size()-1) = old_potentials(old_lengths.size()-1);
          for (unsigned int i=1; i<old_lengths.size()-1;++i)
              {
              unsigned int jj=1000000;
              for (unsigned int j=1; j<old_lengths.size();++j)
                  {
                  if (new_lengths(i) < old_lengths(j))
                     {
                     jj = j;
                     break;
                     }
                  }
              double fraction = (new_lengths(i)-old_lengths(jj-1))/(old_lengths(jj)-old_lengths(jj-1));
              solution(indices[i]+vector_dh.n_dofs()) = old_potentials(jj-1)+(old_potentials(jj)-old_potentials(jj-1))*fraction;
              //solution_dot(3*indices[i]) = old_vx(jj-1)+(old_vx(jj)-old_vx(jj-1))*fraction;
              //solution_dot(3*indices[i]+1) = old_vy(jj-1)+(old_vy(jj)-old_vy(jj-1))*fraction;
              //cout<<i<<" ---> "<<jj<<" "<<fraction<<" "<<old_potentials(jj-1)<<" "<<old_potentials(jj)<<" "<<new_potentials(i)<<endl;
              }
          }

      for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
          {
          if ( (!(comp_dom.flags[i] & near_boat)==0) &&
               (!(comp_dom.flags[i] & near_water)==0))
             {
             solution_dot(3*i) = 0;
             solution_dot(3*i+1) = 0;
             }
          }

      //std::string filename2 = ( "postPostPostRemesh.vtu" );
      //output_results(filename2, t, solution, solution_dot);

std::cout<<"...Done restoring mesh conformity"<<std::endl;
}



template <int dim>
void FreeSurface<dim>::remove_transom_hanging_nodes(Vector<double> & solution,
				                    Vector<double> &solution_dot,
				                    const double t,
	                                            const unsigned int step_number,
		                                    const double  h)
{
std::cout<<"Removing hanging nodes from transom stern..."<<std::endl;

  Point<3> hull_lin_vel;
  Point<3> hull_lin_displ;
  Point<3> hull_lin_vel_dot;
  Point<3> hull_lin_displ_dot;
  Point<3> hull_ang_vel;
  Point<3> hull_ang_vel_dot;
  Point<3> hull_quat_vect;
  Point<3> hull_quat_vect_dot;
  double hull_quat_scal;
  double hull_quat_scal_dot;
  for (unsigned int d=0; d<3;++d)
      {
      hull_lin_vel(d) = solution(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+d);
      hull_lin_displ(d) = solution(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+3+d);
      hull_lin_vel_dot(d) = solution_dot(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+d);
      hull_lin_displ_dot(d) = solution_dot(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+3+d);
      hull_ang_vel(d) = solution(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+6+d);
      hull_ang_vel_dot(d) = solution_dot(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+6+d);
      hull_quat_vect(d) = solution(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+9+d);
      hull_quat_vect_dot(d) = solution_dot(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+9+d);
      }
  hull_quat_scal = solution(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+12);
  hull_quat_scal_dot = solution_dot(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+12);

  cout<<"*Hull Rigid Displacement: "<<hull_lin_displ<<endl;
  cout<<"*Hull Rigid Velocity: "<<hull_lin_vel<<endl;
  cout<<"*Hull Angualr Velocity: "<<hull_lin_vel<<endl;


     Vector<double> rigid_motion_velocities(comp_dom.vector_dh.n_dofs());

     // this moves the rigid_motion_map_points vector (target positions) for the internal nodes of the hull 
     // mesh (all nodes except for the free surface ones)
     for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
         {
         if ( ((comp_dom.flags[i] & boat) &&
              !(comp_dom.flags[i] & near_water) ) ||
              (comp_dom.flags[i] & transom_on_water) )
            {
            //if (fabs(t-0.2) <1e-5)
            //cout<<"BEFORE: "<<comp_dom.ref_points[3*i](0)+comp_dom.old_map_points(3*i)<<" "
            //                <<comp_dom.ref_points[3*i](1)+comp_dom.old_map_points(3*i+1)<<" "
            //                <<comp_dom.ref_points[3*i](2)+comp_dom.old_map_points(3*i+2)<<endl;
            gp_Pnt original_boat_mesh_point = Pnt(Point<3>(comp_dom.ref_points[3*i](0)+comp_dom.old_map_points(3*i),
                                                           comp_dom.ref_points[3*i](1)+comp_dom.old_map_points(3*i+1),
                                                           comp_dom.ref_points[3*i](2)+comp_dom.old_map_points(3*i+2)));
            gp_Pnt boat_mesh_point = original_boat_mesh_point;
            // we first take this point (which is in the RESTART hull location) and transform it to be in the
            // REFERENCE configuration 
            boat_mesh_point.Transform(restart_hull_location.Inverted());

            // now we use sacado to compute the residual at this dof, along with its derivatives with respect to 
            // the 7 dofs associated to the rigid linear and angular displacements

            double s_x = hull_lin_displ(0);
            double s_y = hull_lin_displ(1);
            double s_z = hull_lin_displ(2);
            double v_x = hull_quat_vect(0);
            double  v_y = hull_quat_vect(1);
            double v_z = hull_quat_vect(2);
            double s = hull_quat_scal;

            Point<3> ref_baricenter_pos(comp_dom.boat_model.reference_hull_baricenter(0),
                                        comp_dom.boat_model.reference_hull_baricenter(1),
                                        comp_dom.boat_model.reference_hull_baricenter(2));
   
            Point<3> baricenter_pos(s_x+comp_dom.boat_model.reference_hull_baricenter(0),
                                    s_y+comp_dom.boat_model.reference_hull_baricenter(1),
                                    s_z+comp_dom.boat_model.reference_hull_baricenter(2));

            Point<3> RotMatRow1(1-2*v_y*v_y-2*v_z*v_z, 2*v_x*v_y-2*s*v_z, 2*v_x*v_z+2*s*v_y);
            Point<3> RotMatRow2(2*v_x*v_y+2*s*v_z, 1-2*v_x*v_x-2*v_z*v_z, 2*v_y*v_z-2*s*v_x);
            Point<3> RotMatRow3(2*v_x*v_z-2*s*v_y, 2*v_y*v_z+2*s*v_x, 1-2*v_y*v_y-2*v_x*v_x);

            Point<3> ref_point_pos(boat_mesh_point.X(),boat_mesh_point.Y(),boat_mesh_point.Z());
            Point<3> rigid_lin_displ(s_x,s_y,s_z);
            Point<3> target_point_pos(RotMatRow1*(ref_point_pos+(-1.0)*ref_baricenter_pos),
                                      RotMatRow2*(ref_point_pos+(-1.0)*ref_baricenter_pos),
                                      RotMatRow3*(ref_point_pos+(-1.0)*ref_baricenter_pos));
            target_point_pos += baricenter_pos;
            // now we have the point on the REFERENCE hull, and transform it to go onto the CURRENT hull
            //boat_mesh_point.Transform(reference_to_current_transformation);
            // the rigid motion map points is the difference between the reference point position and the
            // current (rigidly displaced) node position
            comp_dom.rigid_motion_map_points(3*i) = target_point_pos(0)-ref_point_pos(0);
            comp_dom.rigid_motion_map_points(3*i+1) = target_point_pos(1)-ref_point_pos(1);
            comp_dom.rigid_motion_map_points(3*i+2) = target_point_pos(2)-ref_point_pos(2);

            rigid_motion_velocities(3*i) = hull_lin_vel(0)+
                                           hull_ang_vel(1)*(target_point_pos(2)-baricenter_pos(2))-
                                           hull_ang_vel(2)*(target_point_pos(1)-baricenter_pos(1));
            rigid_motion_velocities(3*i+1) = hull_lin_vel(1)+
                                             hull_ang_vel(2)*(target_point_pos(0)-baricenter_pos(0))-
                                             hull_ang_vel(0)*(target_point_pos(2)-baricenter_pos(2));
            rigid_motion_velocities(3*i+2) = hull_lin_vel(2)+
                                             hull_ang_vel(0)*(target_point_pos(1)-baricenter_pos(1))-
                                             hull_ang_vel(1)*(target_point_pos(0)-baricenter_pos(0));
            
            

            //if (fabs(t-0.2) <1e-5)
            //cout<<"AFTER: "<<Pnt(boat_mesh_point)<<" vs "<<Pnt(original_boat_mesh_point)<<endl;
            //cout<<"RMMP: "<<comp_dom.rigid_motion_map_points(3*i)<<" "
            //              <<comp_dom.rigid_motion_map_points(3*i+1)<<" "
            //              <<comp_dom.rigid_motion_map_points(3*i+2)<<endl;
            //cout<<"NN: "<<comp_dom.rigid_motion_map_points(3*i)+nodes_positions(3*i)+comp_dom.ref_points[3*i](0)<<" "
            //            <<comp_dom.rigid_motion_map_points(3*i+1)+nodes_positions(3*i+1)+comp_dom.ref_points[3*i](1)<<" "
            //            <<comp_dom.rigid_motion_map_points(3*i+2)+nodes_positions(3*i+2)+comp_dom.ref_points[3*i](2)<<endl;
            }
         }


      Triangulation<dim-1, dim> &tria = comp_dom.tria;
      DoFHandler<dim-1, dim> &dh = comp_dom.dh;
      DoFHandler<dim-1, dim> &vector_dh = comp_dom.vector_dh;
    

      std::vector<Point<dim> > support_points(dh.n_dofs());
      DoFTools::map_dofs_to_support_points<dim-1, dim>( *comp_dom.mapping, dh, support_points);



      VectorView<double> Phi(dh.n_dofs(), solution.begin()+vector_dh.n_dofs());
      VectorView<double> Phi_dot(dh.n_dofs(), solution_dot.begin()+vector_dh.n_dofs());
      VectorView<double> dphi_dn(dh.n_dofs(), solution.begin()+vector_dh.n_dofs()+dh.n_dofs());
      VectorView<double> dphi_dn_dot(dh.n_dofs(), solution_dot.begin()+vector_dh.n_dofs()+dh.n_dofs());

      std::vector<Point<3> > old_points(dh.n_dofs());
      old_points = support_points;

      Vector<double> curvatures(vector_dh.n_dofs());
      comp_dom.surface_smoother->compute_curvatures(curvatures);

				   // Get support points in the
				   // reference configuration
      std::vector<Point<3> > ref_points(dh.n_dofs());
      DoFTools::map_dofs_to_support_points<2,3>(StaticMappingQ1<2,3>::mapping,
					  dh, ref_points);
    unsigned int refinedCellCounter = 1;

    while(refinedCellCounter)
     {
     refinedCellCounter = 0;
     for (unsigned int i=0;i<comp_dom.dh.n_dofs();++i)
         {
         if ((comp_dom.flags[i] & transom_on_water) )
         //if ((comp_dom.flags[i] & water) && (comp_dom.flags[i] & near_boat))
            {
            //cout<<i<<": "<<support_points[i]<<endl;
            std::vector<cell_it>  cells = comp_dom.dof_to_elems[i];
            for (unsigned int k=0; k<cells.size(); ++k)
                {
                //cout<<k<<":  "<<cells[k]<<"   ("<<cells[k]->center()<<")"<<endl;
                for (unsigned int j=0; j<GeometryInfo<2>::faces_per_cell; ++j)
                    {
                    //cout<<"j: "<<j<<"  nb: "<<cells[k]->neighbor_index(j)<<"  ("<<endl;
                    if (cells[k]->neighbor_index(j) != -1)
                       if ( cells[k]->neighbor(j)->at_boundary() &&
                            (cells[k]->neighbor_is_coarser(j) || cells[k]->refine_flag_set() ) &&
                            !(cells[k]->neighbor(j)->refine_flag_set()) )
                          {
                          //cout<<"FOUND: "<<cells[k]->neighbor(j)<<" ("<<cells[k]->neighbor(j)->center()<<")"<<endl;
                          cells[k]->neighbor(j)->set_refine_flag();
                          refinedCellCounter++;
                          }
                    }
                }
            }
         }

     }

      SolutionTransfer<dim-1, Vector<double>, DoFHandler<dim-1, dim> > soltrans(dh);

      Vector<double> positions_x(comp_dom.dh.n_dofs());
      Vector<double> positions_y(comp_dom.dh.n_dofs());
      Vector<double> positions_z(comp_dom.dh.n_dofs());
      Vector<double> positions_dot_x(comp_dom.dh.n_dofs());
      Vector<double> positions_dot_y(comp_dom.dh.n_dofs());
      Vector<double> positions_dot_z(comp_dom.dh.n_dofs());

      for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
          {
          positions_x(i) = solution(3*i+0)+comp_dom.rigid_motion_map_points(3*i+0)+comp_dom.ref_points[3*i](0);
          positions_y(i) = solution(3*i+1)+comp_dom.rigid_motion_map_points(3*i+1)+comp_dom.ref_points[3*i](1);
          positions_z(i) = solution(3*i+2)+comp_dom.rigid_motion_map_points(3*i+2)+comp_dom.ref_points[3*i](2);
          positions_dot_x(i) = solution_dot(3*i+0);
          positions_dot_y(i) = solution_dot(3*i+1);
          positions_dot_z(i) = solution_dot(3*i+2);
          }          

      std::vector<Vector<double> > all_in;
      all_in.push_back((Vector<double>)Phi);
      all_in.push_back((Vector<double>)Phi_dot);
      all_in.push_back((Vector<double>)dphi_dn);
      all_in.push_back((Vector<double>)dphi_dn_dot);
      all_in.push_back(positions_x);
      all_in.push_back(positions_y);
      all_in.push_back(positions_z);
      all_in.push_back(positions_dot_x);
      all_in.push_back(positions_dot_y);
      all_in.push_back(positions_dot_z); 

      tria.prepare_coarsening_and_refinement();
      soltrans.prepare_for_coarsening_and_refinement(all_in);



      //std::cout << "Refined counter: " << refinedCellCounter << std::endl;
      tria.execute_coarsening_and_refinement();
      dh.distribute_dofs(comp_dom.fe);
      vector_dh.distribute_dofs(comp_dom.vector_fe);  
      comp_dom.map_points.reinit(vector_dh.n_dofs());
      comp_dom.smoothing_map_points.reinit(vector_dh.n_dofs());
      comp_dom.old_map_points.reinit(vector_dh.n_dofs());
      comp_dom.rigid_motion_map_points.reinit(vector_dh.n_dofs());
      comp_dom.initial_map_points.reinit(vector_dh.n_dofs());
      comp_dom.ref_points.resize(vector_dh.n_dofs());
      DoFTools::map_dofs_to_support_points<2,3>(StaticMappingQ1<2,3>::mapping,
					      comp_dom.vector_dh, comp_dom.ref_points);
      comp_dom.generate_double_nodes_set();

      compute_constraints(constraints, vector_constraints);

      dofs_number = vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()+13;

      std::cout<<"Total number of dofs after fixing transom stern: "<<dh.n_dofs()<<std::endl;


      Vector<double> new_Phi(dh.n_dofs());
      Vector<double> new_Phi_dot(dh.n_dofs());

      Vector <double> new_dphi_dn(dh.n_dofs());     
      Vector <double> new_dphi_dn_dot(dh.n_dofs());
   
      Vector<double> new_positions_x(dh.n_dofs());
      Vector<double> new_positions_y(dh.n_dofs());
      Vector<double> new_positions_z(dh.n_dofs());
      Vector<double> new_positions_dot_x(dh.n_dofs());
      Vector<double> new_positions_dot_y(dh.n_dofs());
      Vector<double> new_positions_dot_z(dh.n_dofs());

      std::vector<Vector<double> > all_out;
      all_out.push_back(new_Phi);
      all_out.push_back(new_Phi_dot);
      all_out.push_back(new_dphi_dn);
      all_out.push_back(new_dphi_dn_dot);
      all_out.push_back(new_positions_x);
      all_out.push_back(new_positions_y);
      all_out.push_back(new_positions_z);
      all_out.push_back(new_positions_dot_x);
      all_out.push_back(new_positions_dot_y);
      all_out.push_back(new_positions_dot_z);

      soltrans.interpolate(all_in, all_out);
  
      solution.reinit(dofs_number);
      solution_dot.reinit(dofs_number); 


      constraints.distribute(all_out[0]);
      constraints.distribute(all_out[1]);
      constraints.distribute(all_out[2]);
      constraints.distribute(all_out[3]);      
      constraints.distribute(all_out[4]);
      constraints.distribute(all_out[5]);
      constraints.distribute(all_out[6]);
      constraints.distribute(all_out[7]);    
      constraints.distribute(all_out[8]);
      constraints.distribute(all_out[9]);

      // we have to compute rigid_motion_map_points and rigid_motion_velocities on the new mesh

     rigid_motion_velocities.reinit(comp_dom.vector_dh.n_dofs());

     // this moves the rigid_motion_map_points vector (target positions) for the internal nodes of the hull 
     // mesh (all nodes except for the free surface ones)
     for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
         {
         if ( ((comp_dom.flags[i] & boat) &&
              !(comp_dom.flags[i] & near_water) ) ||
              (comp_dom.flags[i] & transom_on_water) )
            {
            //if (fabs(t-0.2) <1e-5)
            //cout<<"BEFORE: "<<comp_dom.ref_points[3*i](0)+comp_dom.old_map_points(3*i)<<" "
            //                <<comp_dom.ref_points[3*i](1)+comp_dom.old_map_points(3*i+1)<<" "
            //                <<comp_dom.ref_points[3*i](2)+comp_dom.old_map_points(3*i+2)<<endl;
            gp_Pnt original_boat_mesh_point = Pnt(Point<3>(comp_dom.ref_points[3*i](0)+comp_dom.old_map_points(3*i),
                                                           comp_dom.ref_points[3*i](1)+comp_dom.old_map_points(3*i+1),
                                                           comp_dom.ref_points[3*i](2)+comp_dom.old_map_points(3*i+2)));
            gp_Pnt boat_mesh_point = original_boat_mesh_point;
            // we first take this point (which is in the RESTART hull location) and transform it to be in the
            // REFERENCE configuration 
            boat_mesh_point.Transform(restart_hull_location.Inverted());

            // now we use sacado to compute the residual at this dof, along with its derivatives with respect to 
            // the 7 dofs associated to the rigid linear and angular displacements

            double s_x = hull_lin_displ(0);
            double s_y = hull_lin_displ(1);
            double s_z = hull_lin_displ(2);
            double v_x = hull_quat_vect(0);
            double  v_y = hull_quat_vect(1);
            double v_z = hull_quat_vect(2);
            double s = hull_quat_scal;

            Point<3> ref_baricenter_pos(comp_dom.boat_model.reference_hull_baricenter(0),
                                        comp_dom.boat_model.reference_hull_baricenter(1),
                                        comp_dom.boat_model.reference_hull_baricenter(2));
   
            Point<3> baricenter_pos(s_x+comp_dom.boat_model.reference_hull_baricenter(0),
                                    s_y+comp_dom.boat_model.reference_hull_baricenter(1),
                                    s_z+comp_dom.boat_model.reference_hull_baricenter(2));

            Point<3> RotMatRow1(1-2*v_y*v_y-2*v_z*v_z, 2*v_x*v_y-2*s*v_z, 2*v_x*v_z+2*s*v_y);
            Point<3> RotMatRow2(2*v_x*v_y+2*s*v_z, 1-2*v_x*v_x-2*v_z*v_z, 2*v_y*v_z-2*s*v_x);
            Point<3> RotMatRow3(2*v_x*v_z-2*s*v_y, 2*v_y*v_z+2*s*v_x, 1-2*v_y*v_y-2*v_x*v_x);

            Point<3> ref_point_pos(boat_mesh_point.X(),boat_mesh_point.Y(),boat_mesh_point.Z());
            Point<3> rigid_lin_displ(s_x,s_y,s_z);
            Point<3> target_point_pos(RotMatRow1*(ref_point_pos+(-1.0)*ref_baricenter_pos),
                                      RotMatRow2*(ref_point_pos+(-1.0)*ref_baricenter_pos),
                                      RotMatRow3*(ref_point_pos+(-1.0)*ref_baricenter_pos));
            target_point_pos += baricenter_pos;
            // now we have the point on the REFERENCE hull, and transform it to go onto the CURRENT hull
            //boat_mesh_point.Transform(reference_to_current_transformation);
            // the rigid motion map points is the difference between the reference point position and the
            // current (rigidly displaced) node position
            comp_dom.rigid_motion_map_points(3*i) = target_point_pos(0)-ref_point_pos(0);
            comp_dom.rigid_motion_map_points(3*i+1) = target_point_pos(1)-ref_point_pos(1);
            comp_dom.rigid_motion_map_points(3*i+2) = target_point_pos(2)-ref_point_pos(2);

            rigid_motion_velocities(3*i) = hull_lin_vel(0)+
                                           hull_ang_vel(1)*(target_point_pos(2)-baricenter_pos(2))-
                                           hull_ang_vel(2)*(target_point_pos(1)-baricenter_pos(1));
            rigid_motion_velocities(3*i+1) = hull_lin_vel(1)+
                                             hull_ang_vel(2)*(target_point_pos(0)-baricenter_pos(0))-
                                             hull_ang_vel(0)*(target_point_pos(2)-baricenter_pos(2));
            rigid_motion_velocities(3*i+2) = hull_lin_vel(2)+
                                             hull_ang_vel(0)*(target_point_pos(1)-baricenter_pos(1))-
                                             hull_ang_vel(1)*(target_point_pos(0)-baricenter_pos(0));
            
            

            //if (fabs(t-0.2) <1e-5)
            //cout<<"AFTER: "<<Pnt(boat_mesh_point)<<" vs "<<Pnt(original_boat_mesh_point)<<endl;
            //cout<<"RMMP: "<<comp_dom.rigid_motion_map_points(3*i)<<" "
            //              <<comp_dom.rigid_motion_map_points(3*i+1)<<" "
            //              <<comp_dom.rigid_motion_map_points(3*i+2)<<endl;
            //cout<<"NN: "<<comp_dom.rigid_motion_map_points(3*i)+nodes_positions(3*i)+comp_dom.ref_points[3*i](0)<<" "
            //            <<comp_dom.rigid_motion_map_points(3*i+1)+nodes_positions(3*i+1)+comp_dom.ref_points[3*i](1)<<" "
            //            <<comp_dom.rigid_motion_map_points(3*i+2)+nodes_positions(3*i+2)+comp_dom.ref_points[3*i](2)<<endl;
            }
         }

      for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i) 
	  {
          solution(3*i+0) = all_out[4](i)-comp_dom.rigid_motion_map_points(3*i)-comp_dom.ref_points[3*i](0);
          solution(3*i+1) = all_out[5](i)-comp_dom.rigid_motion_map_points(3*i+1)-comp_dom.ref_points[3*i](1);
          solution(3*i+2) = all_out[6](i)-comp_dom.rigid_motion_map_points(3*i+2)-comp_dom.ref_points[3*i](2);
          solution_dot(3*i+0) = all_out[7](i);
          solution_dot(3*i+1) = all_out[8](i);
          solution_dot(3*i+2) = all_out[9](i);
	  solution(i+vector_dh.n_dofs()) = all_out[0](i);
          solution_dot(i+vector_dh.n_dofs()) = all_out[1](i);
          solution(i+vector_dh.n_dofs()+dh.n_dofs()) = all_out[2](i);
          solution_dot(i+vector_dh.n_dofs()+dh.n_dofs()) = all_out[3](i);
	  }

      for (unsigned int k=0; k<3; ++k)
          {
          solution(k+3+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = hull_lin_displ(k);
          solution_dot(k+3+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = hull_lin_displ_dot(k);
          solution(k+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = hull_lin_vel(k);
          solution_dot(k+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = hull_lin_vel_dot(k);
          solution(k+6+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = hull_ang_vel(k);
          solution_dot(k+6+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = hull_ang_vel_dot(k);
          solution(k+9+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = hull_quat_vect(k);
          solution_dot(k+9+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = hull_quat_vect_dot(k);          
          }
      solution(12+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = hull_quat_scal;
      solution_dot(12+vector_dh.n_dofs()+dh.n_dofs()+dh.n_dofs()) = hull_quat_scal_dot; 



      DXDt_and_DphiDt_vector.reinit(vector_dh.n_dofs()+dh.n_dofs());

      DphiDt_sparsity_pattern.reinit (dh.n_dofs(),
				      dh.n_dofs(),
				      dh.max_couplings_between_dofs());
      vector_sparsity_pattern.reinit (vector_dh.n_dofs(),
					vector_dh.n_dofs(),
					vector_dh.max_couplings_between_dofs());

      DoFTools::make_sparsity_pattern (dh, DphiDt_sparsity_pattern, constraints);
      DphiDt_sparsity_pattern.compress();

      DoFTools::make_sparsity_pattern (vector_dh, vector_sparsity_pattern, vector_constraints);
      vector_sparsity_pattern.compress();

      working_map_points.reinit(comp_dom.vector_dh.n_dofs());
      working_nodes_velocities.reinit(comp_dom.vector_dh.n_dofs());
      nodes_pos_res.reinit(comp_dom.vector_dh.n_dofs());
      nodes_ref_surf_dist.reinit(comp_dom.vector_dh.n_dofs());
      nodes_diff_jac_x_delta.reinit(dofs_number);
      nodes_alg_jac_x_delta.reinit(dofs_number);
      dae_nonlin_residual.reinit(dofs_number);
      dae_linear_step_residual.reinit(dofs_number);
      current_sol.reinit(dofs_number);
      current_sol_dot.reinit(dofs_number);
      bem_residual.reinit(comp_dom.dh.n_dofs()); 
      bem_phi.reinit(comp_dom.dh.n_dofs());
      bem_dphi_dn.reinit(comp_dom.dh.n_dofs());
      temp_src.reinit(comp_dom.vector_dh.n_dofs());
      break_wave_press.reinit(comp_dom.dh.n_dofs());
  
      

      bem.reinit();

      support_points.resize(dh.n_dofs());
      DoFTools::map_dofs_to_support_points<2, 3>( *comp_dom.mapping, dh, support_points);
      std::vector<bool> new_boundary_dofs(vector_dh.n_dofs());
      std::vector< bool > comp_sel(3, true);
      DoFTools::extract_boundary_dofs(vector_dh, comp_sel, new_boundary_dofs);
      for (unsigned int i=0; i<dh.n_dofs();++i)
          for (unsigned int j=0; j<old_points.size();++j)
              if (old_points[j].distance(support_points[i]) < 1e-5)
                 {
                 new_boundary_dofs[3*i] = 1;
                 new_boundary_dofs[3*i+1] = 1;
                 new_boundary_dofs[3*i+2] = 1;
                 }

      for (unsigned int i=0; i<comp_dom.vector_dh.n_dofs(); ++i)
          {
          comp_dom.map_points(i) = solution(i)+comp_dom.rigid_motion_map_points(i);
          }
       
      //std::string filename2 = ( "beforeCrash.vtu" );
      //output_results(filename2, t, solution, solution_dot);

      if (!comp_dom.no_boat)
         comp_dom.evaluate_ref_surf_distances(nodes_ref_surf_dist,false);
      comp_dom.map_points -= nodes_ref_surf_dist;
      for (unsigned int i=0; i<comp_dom.vector_dh.n_dofs(); ++i)
          {
          solution(i) = comp_dom.map_points(i)-comp_dom.rigid_motion_map_points(i);
          }




      // in particular we must get the position of the nodes (in terms of curvilinear length)
      // on the smoothing lines, and the corresponding potential and horizontal velcoity values, in order to
      // interpolate the new values to be assigned at the restart of the simulation

      for (unsigned smooth_id=0; smooth_id<comp_dom.line_smoothers.size(); ++smooth_id)
          {
          Vector<double>  &old_lengths = comp_dom.line_smoothers[smooth_id]->get_lengths_before_smoothing();
          Vector<double>  &new_lengths = comp_dom.line_smoothers[smooth_id]->get_lengths_after_smoothing();
          std::vector<unsigned int> &indices = comp_dom.line_smoothers[smooth_id]->get_node_indices();
          Vector<double> old_potentials(old_lengths.size());
          Vector<double> old_vx(old_lengths.size());
          Vector<double> old_vy(old_lengths.size());
          //Vector<double> new_potentials(old_lengths.size());
          for (unsigned int i=0; i<old_lengths.size();++i)
              {
              old_potentials(i) = solution(indices[i]+vector_dh.n_dofs());
              old_vx(i) = solution_dot(3*indices[i]);
              old_vy(i) = solution_dot(3*indices[i]+1);
              //cout<<i<<"("<<indices[i]<<"->"<<round(indices[i]/3)<<")    "<<old_lengths(i)<<" vs "<<new_lengths(i)<<"  pot: "<<old_potentials(i)<<endl;
              //cout<<indices[i]<<" "<<comp_dom.support_points[indices[i]]<<endl;
              }
          //new_potentials(0) = old_potentials(0);
          //new_potentials(old_lengths.size()-1) = old_potentials(old_lengths.size()-1);
          for (unsigned int i=1; i<old_lengths.size()-1;++i)
              {
              unsigned int jj=1000000;
              for (unsigned int j=1; j<old_lengths.size();++j)
                  {
                  if (new_lengths(i) < old_lengths(j))
                     {
                     jj = j;
                     break;
                     }
                  }
              double fraction = (new_lengths(i)-old_lengths(jj-1))/(old_lengths(jj)-old_lengths(jj-1));
              solution(indices[i]+vector_dh.n_dofs()) = old_potentials(jj-1)+(old_potentials(jj)-old_potentials(jj-1))*fraction;
              //solution_dot(3*indices[i]) = old_vx(jj-1)+(old_vx(jj)-old_vx(jj-1))*fraction;
              //solution_dot(3*indices[i]+1) = old_vy(jj-1)+(old_vy(jj)-old_vy(jj-1))*fraction;
              //cout<<i<<" ---> "<<jj<<" "<<fraction<<" "<<old_potentials(jj-1)<<" "<<old_potentials(jj)<<" "<<new_potentials(i)<<endl;
              }
          }

      for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
          {
          if ( (!(comp_dom.flags[i] & near_boat)==0) &&
               (!(comp_dom.flags[i] & near_water)==0))
             {
             solution_dot(3*i) = 0;
             solution_dot(3*i+1) = 0;
             }
          }

      //std::string filename2 = ( "postPostPostRemesh.vtu" );
      //output_results(filename2, t, solution, solution_dot);

std::cout<<"...Done removing hanging nodes from transom stern"<<std::endl;
}


template <int dim>
void FreeSurface<dim>::prepare_restart(const double t, Vector<double> &y, Vector<double> &yp, bool restart_flag)
{
std::cout<<"Preparing interpolated solution for restart"<<std::endl;



     Vector<double> res;
     res.reinit(sys_comp.size());
     //yp.reinit(sys_comp.size());
//cout<<"CHECK1"<<endl;
//residual(t,res,y,yp);
//setup_jacobian_prec(t,y,yp,0.0);





  // first of all, all yp for coordinates is to be set to zero: at restarts all is idle  
  if (restart_flag)
     for (unsigned int i=0; i<comp_dom.vector_dh.n_dofs(); ++i)
         yp(i) = 0;




  const VectorView<double> nodes_positions(comp_dom.vector_dh.n_dofs(),y.begin());
  const VectorView<double> nodes_velocities(comp_dom.vector_dh.n_dofs(),yp.begin());
  const VectorView<double> phi(comp_dom.dh.n_dofs(),y.begin()+comp_dom.vector_dh.n_dofs());
  const VectorView<double> phi_time_derivs(comp_dom.dh.n_dofs(),yp.begin()+comp_dom.vector_dh.n_dofs());
  const VectorView<double> dphi_dn(comp_dom.dh.n_dofs(),y.begin()+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs());
  const VectorView<double> dphi_dn_time_derivs(comp_dom.dh.n_dofs(),yp.begin()+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs());

  Point<3> hull_lin_vel;
  Point<3> hull_lin_displ;
  Point<3> hull_lin_vel_dot;
  Point<3> hull_lin_displ_dot;
  Point<3> hull_ang_vel;
  Point<3> hull_quat_vector;
  Point<3> hull_ang_vel_dot;
  Point<3> hull_quat_vector_dot;
  for (unsigned int d=0; d<3;++d)
      {
      hull_lin_vel(d) = y(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+d);
      hull_lin_displ(d) = y(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+3+d);
      hull_lin_vel_dot(d) = yp(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+d);
      hull_lin_displ_dot(d) = yp(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+3+d);
      hull_ang_vel(d) = y(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+6+d);
      hull_quat_vector(d) = y(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+9+d);
      hull_ang_vel_dot(d) = yp(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+6+d);
      hull_quat_vector_dot(d) = yp(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+9+d);
      }
  double hull_quat_scalar = y(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+12);
  double hull_quat_scalar_dot = yp(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+12);

  cout<<"*Hull Rigid Displacement: "<<hull_lin_displ<<endl;
  cout<<"*Hull Rigid Velocity: "<<hull_lin_vel<<endl;
  cout<<"*Hull Angular Velocity: "<<hull_ang_vel<<endl;   
  //let's start with rigid modes variables
  for (unsigned int d=0; d<3;++d)
      {
      yp(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+3+d) = hull_lin_vel(d);
      }

//  for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
//         cout<<i<<"  Poss: "<<" "<<nodes_positions(3*i)<<" "<<nodes_positions(3*i+1)<<" "<<nodes_positions(3*i+2)<<endl;  
//  for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
//         cout<<i<<"  Vels: "<<" "<<nodes_velocities(3*i)<<" "<<nodes_velocities(3*i+1)<<" "<<nodes_velocities(3*i+2)<<endl;
//  for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
//         cout<<i<<"  Pot: "<<" "<<phi(i)<<"  Pod_dot: "<<phi_time_derivs(i)<<"  dphi_dn: "<<dphi_dn(i)<<endl;
  // TO BE REMOVED
  DphiDt_sys_rhs.reinit (comp_dom.dh.n_dofs()); 
  DphiDt_sys_rhs_2.reinit (comp_dom.dh.n_dofs());
  DphiDt_sys_rhs_3.reinit (comp_dom.dh.n_dofs()); 
  DphiDt_sys_rhs_4.reinit (comp_dom.dh.n_dofs()); 
  DphiDt_sys_matrix.reinit(DphiDt_sparsity_pattern);
  DphiDt_sys_matrix_2.reinit(DphiDt_sparsity_pattern);
  
  
  wind.set_time(t);
  Vector<double> instantWindValue(dim);
  Point<dim> zero(0,0,0);
  wind.vector_value(zero,instantWindValue);
  std::cout<<std::endl<<"Simulation time= "<<t<<"   Vinf= ";
  instantWindValue.print(std::cout,4,false,true);
  std::cout << "Ndofs ODE= " << y.size();
  std::cout << "  Ndofs BEM= " << comp_dom.dh.n_dofs();
  std::cout<<std::endl;  
  wind.set_time(t);
// we'll need Vinf
  Vector<double> wind_value(dim);
  wind.vector_value(Point<3>(0.0,0.0,0.0),wind_value);
  Point<dim> Vinf;
  for (unsigned int i = 0; i < dim; i++)
      Vinf(i) = wind_value(i);

  cout<<"Hull Rigid Displacement: "<<hull_lin_displ<<endl;
  cout<<"Hull Rigid Velocity: "<<hull_lin_vel<<endl;

//////////////////////////////////////////////////////////////////////
  // here we take care of the geometric part of the variables
//////////////////////////////////////////////////////////////////////
     // this moves the rigid_motion_map_points vector (target positions) for the internal nodes of the hull 
     // mesh (all nodes except for the free surface ones)
     for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
         {
         if ( ((comp_dom.flags[i] & boat) &&
              !(comp_dom.flags[i] & near_water) ) ||
              (comp_dom.flags[i] & transom_on_water) )
            {
            //if (fabs(t-0.2) <1e-5)
            //cout<<"BEFORE: "<<comp_dom.ref_points[3*i](0)+comp_dom.old_map_points(3*i)<<" "
            //                <<comp_dom.ref_points[3*i](1)+comp_dom.old_map_points(3*i+1)<<" "
            //                <<comp_dom.ref_points[3*i](2)+comp_dom.old_map_points(3*i+2)<<endl;
            gp_Pnt original_boat_mesh_point = Pnt(Point<3>(comp_dom.ref_points[3*i](0)+comp_dom.old_map_points(3*i),
                                                           comp_dom.ref_points[3*i](1)+comp_dom.old_map_points(3*i+1),
                                                           comp_dom.ref_points[3*i](2)+comp_dom.old_map_points(3*i+2)));
            gp_Pnt boat_mesh_point = original_boat_mesh_point;
            // we first take this point (which is in the RESTART hull location) and transform it to be in the
            // REFERENCE configuration 
            boat_mesh_point.Transform(restart_hull_location.Inverted());

            // now we use sacado to compute the residual at this dof, along with its derivatives with respect to 
            // the 7 dofs associated to the rigid linear and angular displacements
            double s_x,s_y,s_z,v_x,v_y,v_z,s;

            s_x = hull_lin_displ(0);
            s_y = hull_lin_displ(1);
            s_z = hull_lin_displ(2);
            v_x = hull_quat_vector(0);
            v_y = hull_quat_vector(1);
            v_z = hull_quat_vector(2);
            s = hull_quat_scalar;

            Point<3> ref_baricenter_pos(comp_dom.boat_model.reference_hull_baricenter(0),
                                        comp_dom.boat_model.reference_hull_baricenter(1),
                                        comp_dom.boat_model.reference_hull_baricenter(2));
   
            Point<3> baricenter_pos(s_x+comp_dom.boat_model.reference_hull_baricenter(0),
                                    s_y+comp_dom.boat_model.reference_hull_baricenter(1),
                                    s_z+comp_dom.boat_model.reference_hull_baricenter(2));

            Point<3> RotMatRow1(1-2*v_y*v_y-2*v_z*v_z, 2*v_x*v_y-2*s*v_z, 2*v_x*v_z+2*s*v_y);
            Point<3> RotMatRow2(2*v_x*v_y+2*s*v_z, 1-2*v_x*v_x-2*v_z*v_z, 2*v_y*v_z-2*s*v_x);
            Point<3> RotMatRow3(2*v_x*v_z-2*s*v_y, 2*v_y*v_z+2*s*v_x, 1-2*v_y*v_y-2*v_x*v_x);

            Point<3> ref_point_pos(boat_mesh_point.X(),boat_mesh_point.Y(),boat_mesh_point.Z());
            Point<3> rigid_lin_displ(s_x,s_y,s_z);
            Point<3> target_point_pos(RotMatRow1*(ref_point_pos+(-1.0)*ref_baricenter_pos),
                                      RotMatRow2*(ref_point_pos+(-1.0)*ref_baricenter_pos),
                                      RotMatRow3*(ref_point_pos+(-1.0)*ref_baricenter_pos));
            target_point_pos += baricenter_pos;
            // now we have the point on the REFERENCE hull, and transform it to go onto the CURRENT hull
            //boat_mesh_point.Transform(reference_to_current_transformation);
            // the rigid motion map points is the difference between the reference point position and the
            // current (rigidly displaced) node position
            comp_dom.rigid_motion_map_points(3*i) = target_point_pos(0)-ref_point_pos(0);
            comp_dom.rigid_motion_map_points(3*i+1) = target_point_pos(1)-ref_point_pos(1);
            comp_dom.rigid_motion_map_points(3*i+2) = target_point_pos(2)-ref_point_pos(2);

            //working_map_points(3*i) = ref_point_pos(0).val()-comp_dom.ref_points[3*i](0)+comp_dom.rigid_motion_map_points(3*i);
            //working_map_points(3*i+1) = ref_point_pos(1).val()-comp_dom.ref_points[3*i](1)+comp_dom.rigid_motion_map_points(3*i+1);
            //working_map_points(3*i+2) = ref_point_pos(2).val()-comp_dom.ref_points[3*i](2)+comp_dom.rigid_motion_map_points(3*i+2);


            //if (fabs(t-0.2) <1e-5)
            //cout<<"AFTER: "<<Pnt(boat_mesh_point)<<" vs "<<Pnt(original_boat_mesh_point)<<endl;
            //cout<<"RMMP: "<<comp_dom.rigid_motion_map_points(3*i)<<" "
            //              <<comp_dom.rigid_motion_map_points(3*i+1)<<" "
            //              <<comp_dom.rigid_motion_map_points(3*i+2)<<endl;
            //cout<<"NN: "<<comp_dom.rigid_motion_map_points(3*i)+nodes_positions(3*i)+comp_dom.ref_points[3*i](0)<<" "
            //            <<comp_dom.rigid_motion_map_points(3*i+1)+nodes_positions(3*i+1)+comp_dom.ref_points[3*i](1)<<" "
            //            <<comp_dom.rigid_motion_map_points(3*i+2)+nodes_positions(3*i+2)+comp_dom.ref_points[3*i](2)<<endl;
            }
         }
  
  //   }
  // INEFFICIENTE!!!!!!!!!!
  


 Vector<double> full_map_points(comp_dom.rigid_motion_map_points);
  full_map_points.add(nodes_positions);
  comp_dom.update_mapping(full_map_points);
  comp_dom.update_support_points();

  
  //we work on a local COPY of map_points
  working_map_points = comp_dom.map_points;
  nodes_pos_res = 0;



  // as for now x and y coordinates of nodes are not moved (no surface smoothing)
  // so on x and y coordinates (except for water nodes) we only need to put a one
  // on the matrix diagonal
  for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
      {
      if ( !(comp_dom.flags[i] & near_boat) &&
           (comp_dom.flags[i] & water) &&
           (comp_dom.flags[i] & edge) &&
           !(comp_dom.flags[i] & transom_on_water) &&
           (!constraints.is_constrained(i)))
         {
         working_map_points(3*i) = comp_dom.old_map_points(3*i);
         working_map_points(3*i+1) = comp_dom.old_map_points(3*i+1);
         if (comp_dom.flags[i] & near_inflow)
            {
            working_map_points(3*i+2) = comp_dom.old_map_points(3*i+2);
            }
         }
      else if ( !(comp_dom.flags[i] & near_water) &&
                (comp_dom.flags[i] & boat) &&
               (!constraints.is_constrained(i)))
         {
         working_map_points(3*i) = comp_dom.old_map_points(3*i);
         working_map_points(3*i+1) = comp_dom.old_map_points(3*i+1);
         working_map_points(3*i+2) = comp_dom.old_map_points(3*i+2);
         //cout<<"&& "<<3*i<<" "<<3*i+1<<" "<<3*i+2<<endl;
         }
      else if ( !(comp_dom.flags[i] & near_water) &&
                !(comp_dom.flags[i] & water) &&
                !(comp_dom.flags[i] & boat) &&
               (!constraints.is_constrained(i)))
         {
         working_map_points(3*i) = comp_dom.old_map_points(3*i);
         working_map_points(3*i+1) = comp_dom.old_map_points(3*i+1);
         working_map_points(3*i+2) = comp_dom.old_map_points(3*i+2);
         //cout<<"&& "<<3*i<<" "<<3*i+1<<" "<<3*i+2<<endl;         
         }
      else if ((comp_dom.flags[i] & transom_on_water) )
         {
         working_map_points(3*i) = comp_dom.old_map_points(3*i);
         working_map_points(3*i+1) = comp_dom.old_map_points(3*i+1);
         working_map_points(3*i+2) = comp_dom.old_map_points(3*i+2);
         }          
      }

  // blending factor is needed to avoid that right after restart
  // the free surface mesh smoothing causes infinite horizontal nodes velocity:
  // here we are at restart, so blending factor is ZERO 
  double blend_factor = 0.0;
  if (restart_flag)
     {
     }
  else
     {

     if (t - last_remesh_time < 0.5*remeshing_period)
        blend_factor = sin(3.141592654*(t-last_remesh_time)/remeshing_period);
     else
        blend_factor = 1.0;
     //std::cout<<"t "<<t<<"  last_remesh_time "<<last_remesh_time<<"  remeshing_period/5 "<<remeshing_period/5<<std::endl;
     }

  std::cout<<"blend_factor = "<<blend_factor<<std::endl;


     //this takes care of the right water line nodes projection (without smoothing)
     for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
         { 
         if ( (comp_dom.flags[i] & water) &&
              (comp_dom.flags[i] & near_boat) &&
              (comp_dom.flags[i] & right_side) &&
              !(comp_dom.flags[i] & transom_on_water) &&
              (comp_dom.moving_point_ids[3] != i) &&
              (comp_dom.moving_point_ids[4] != i) &&
              (comp_dom.moving_point_ids[5] != i) &&
              (comp_dom.moving_point_ids[6] != i) ) // to avoid the bow and stern node
            {//cout<<"**** "<<i<<endl;
            Point<3> proj_node;
            double iges_curvature;
            Point<3> direction(comp_dom.iges_normals[i](0),comp_dom.iges_normals[i](1),0.0);
            //cout<<3*i+1<<"   "<<comp_dom.support_points[i]<<"   ("<<comp_dom.iges_normals[i]<<")"<<endl;
            if (fabs(comp_dom.old_iges_normals[i](0))<sqrt(3)/3*fabs(comp_dom.old_iges_normals[i](1)))
               direction(0) = 0.0;
            else
               direction(1) = 0.0;
            //if (fabs(comp_dom.old_iges_normals[i](0)) > 1e-3)
               //cout<<3*i<<"  dir:    ("<<direction<<")   |    n_i("<<comp_dom.old_iges_normals[i]<<")"<<endl;
            comp_dom.boat_model.boat_water_line_right->assigned_axis_projection_and_diff_forms(proj_node,
                                                                                               comp_dom.iges_normals[i],
                                                                                               comp_dom.iges_mean_curvatures[i],
                                                                                               comp_dom.support_points[i],
                                                                                               direction);  // hor normal dir projection
            //comp_dom.boat_model.boat_water_line_right->axis_projection_and_diff_forms(proj_node,
            //                                                                          comp_dom.iges_normals[i],
            //                                                                          iges_curvature,
            //                                                                          comp_dom.support_points[i]);  // y axis projection
            //working_map_points(3*i) = comp_dom.old_map_points(3*i); // x of the node must not change
            if (fabs(comp_dom.old_iges_normals[i](0))<sqrt(3)/3*fabs(comp_dom.old_iges_normals[i](1)))
               {
               working_map_points(3*i) = comp_dom.old_map_points(3*i); // x of the node must not change
               working_map_points(3*i+1) = proj_node(1) - comp_dom.ref_points[3*i](1);
               }
            else
               {
               working_map_points(3*i) = proj_node(0) - comp_dom.ref_points[3*i](0);
               working_map_points(3*i+1) = comp_dom.old_map_points(3*i+1); // y of the node must not change
               }            
            working_map_points(3*i+2) = proj_node(2) - comp_dom.ref_points[3*i](2);
            //cout<<3*i+1<<"   "<<working_map_points(3*i+1)-comp_dom.map_points(3*i+1)<<"   ("<<iges_normal<<")"<<endl;
            //if (fabs(comp_dom.old_iges_normals[i](0))>sqrt(3)/3*fabs(comp_dom.old_iges_normals[i](1)))
            //   {
            //   cout<<i<<"   "<<proj_node<<"   ("<<comp_dom.support_points[i]<<")"<<endl;
            //   cout<<direction<<" | "<<working_map_points(3*i)<<" "<<working_map_points(3*i+1)<<endl;
            //   }
            }              
         }

     //this takes care of the left water line nodes projection (without smoothing)
     for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
         { 
         if ( (comp_dom.flags[i] & water) &&
              (comp_dom.flags[i] & near_boat) &&
              (comp_dom.flags[i] & left_side) &&
              !(comp_dom.flags[i] & transom_on_water) &&
              (comp_dom.moving_point_ids[3] != i) &&
              (comp_dom.moving_point_ids[4] != i) &&
              (comp_dom.moving_point_ids[5] != i) &&
              (comp_dom.moving_point_ids[6] != i) ) // to avoid the bow and stern node
            {//cout<<"**** "<<i<<endl;
            Point<3> proj_node;
            double iges_curvature;
            Point<3> direction(comp_dom.iges_normals[i](0),comp_dom.iges_normals[i](1),0.0);
            if (fabs(comp_dom.old_iges_normals[i](0))<sqrt(3)/3*fabs(comp_dom.old_iges_normals[i](1)))
               direction(0) = 0.0;
            else
               direction(1) = 0.0;
            //if (fabs(comp_dom.iges_normals[i](0))<0.001)
            //   cout<<3*i<<"  dir:    ("<<direction<<")"<<endl;
            comp_dom.boat_model.boat_water_line_left->assigned_axis_projection_and_diff_forms(proj_node,
                                                                                              comp_dom.iges_normals[i],
                                                                                              comp_dom.iges_mean_curvatures[i],
                                                                                              comp_dom.support_points[i],
                                                                                              direction);  // hor normal dir projection
            //comp_dom.boat_model.boat_water_line_left->axis_projection_and_diff_forms(proj_node,
            //                                                                         comp_dom.iges_normals[i],
            //                                                                         iges_curvature,
            //                                                                         comp_dom.support_points[i]);  // y axis projection
            //working_map_points(3*i) = comp_dom.old_map_points(3*i); // x of the node must not change
            if (fabs(comp_dom.old_iges_normals[i](0))<sqrt(3)/3*fabs(comp_dom.old_iges_normals[i](1)))
               {
               working_map_points(3*i) = comp_dom.old_map_points(3*i); // x of the node must not change
               working_map_points(3*i+1) = proj_node(1) - comp_dom.ref_points[3*i](1);
               }
            else
               {
               working_map_points(3*i) = proj_node(0) - comp_dom.ref_points[3*i](0);
               working_map_points(3*i+1) = comp_dom.old_map_points(3*i+1); // y of the node must not change
               }            
            working_map_points(3*i+2) = proj_node(2) - comp_dom.ref_points[3*i](2);
            //cout<<i<<"   "<<temp_src(3*i+1)<<"   ("<<comp_dom.iges_normals[i]<<")"<<endl;
            //if (fabs(comp_dom.old_iges_normals[i](0))>sqrt(3)/3*fabs(comp_dom.old_iges_normals[i](1)))
            //   {
            //   cout<<i<<"   "<<proj_node<<"   ("<<comp_dom.support_points[i]<<")"<<endl;
            //   cout<<direction<<" | "<<working_map_points(3*i)<<" "<<working_map_points(3*i+1)<<endl;
            //   }
            }              
         }

     //this takes care of the bow and stern nodes 
     if (!comp_dom.no_boat)
        for (unsigned int k=3; k<7; ++k)
            { 
            unsigned int i = comp_dom.moving_point_ids[k];
            {
            Point <3> dP0 = comp_dom.support_points[i];
            Point <3> dP;
         				   //this is the horizontal plane
            Handle(Geom_Plane) horPlane = new Geom_Plane(0.,0.,1.,-dP0(2));
            Handle(Geom_Curve) curve;
            TopLoc_Location L = comp_dom.boat_model.current_loc;
            TopLoc_Location L_inv = L.Inverted();
            horPlane->Transform(L_inv.Transformation());
            if (comp_dom.boat_model.is_transom)
               {
               if (k==3 || k==4)
                  curve = comp_dom.boat_model.equiv_keel_bspline;
               else if (k == 6)
                  curve = comp_dom.boat_model.left_transom_bspline;
               else
                  curve = comp_dom.boat_model.right_transom_bspline;
               }
            else
               {
               curve = comp_dom.boat_model.equiv_keel_bspline;
               }

            TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(curve);
            edge.Location(L);
            BRepAdaptor_Curve AC(edge);
            gp_Pnt P;
            gp_Vec V1;
            GeomAPI_IntCS Intersector(curve, horPlane);
            int npoints = Intersector.NbPoints();

            AssertThrow((npoints != 0), ExcMessage("Keel or transom curve is not intersecting with horizontal plane!"));
            double minDistance=1e7;
            double t,u,v;
            for (int j=0; j<npoints;++j)
                {
                gp_Pnt int_point = Intersector.Point(j+1);
                int_point.Transform(L.Transformation());
                Point<3> inters = Pnt(int_point);
                Intersector.Parameters(j+1,u,v,t);
                if (dP0.distance(inters) < minDistance)
                   {
                   minDistance = dP0.distance(inters);
                   dP = inters;
                   AC.D1(t,P,V1);
                   }
                }
            //cout<<"Check plane-curve intersection:"<<endl;
            //cout<<"Origin: "<<dP0<<"   Proj: "<<dP<<"  dist: "<<minDistance<<endl;
            //cout<<Pnt(P)<<endl;
            /*
            // here temporarily for kcs hull tests
            if (minDistance > 0.5*comp_dom.boat_model.boatWetLength)
               {
               Standard_Real First = curve->FirstParameter();
               Standard_Real Last = curve->LastParameter();
               gp_Pnt PIn(0.0,0.0,0.0);
               gp_Pnt PFin(0.0,0.0,0.0);
               gp_Vec VIn;
               gp_Vec VFin;
               curve->D1(First,PIn,VIn);
               curve->D1(Last,PFin,VFin);
               cout<<"New part one: "<<Pnt(PIn)<<" | "<<Pnt(PFin)<<endl;
               if (dP0.distance(Pnt(PIn)) < dP0.distance(Pnt(PFin)))
                  {
                  double delta_z = dP0(2) - PIn.Z();
                  dP = Point<3>(PIn.X()+delta_z*VIn.X()/VIn.Z(),PIn.Y()+delta_z*VIn.Y()/VIn.Z(),dP0(2));
                  V1 = VIn;
                  }
               else
                  {
                  double delta_z = dP0(2) - PFin.Z();
                  dP = Point<3>(PFin.X()+delta_z*VFin.X()/VFin.Z(),PIn.Y()+delta_z*VFin.Y()/VFin.Z(),dP0(2));
                  V1 = VFin;
                  }
               cout<<"New part two: "<<dP<<" | "<<V1.X()<<" "<<V1.Y()<<" "<<V1.Z()<<" | "<<dP0<<endl;
               }
            */
            //cout<<k<<"("<<i<<") ---> ("<<dP0<<") vs ("<<dP<<")"<<endl;
            working_map_points(3*i) = dP(0)-comp_dom.ref_points[3*i](0);
            working_map_points(3*i+1) = dP(1)-comp_dom.ref_points[3*i](1);
            working_map_points(3*i+2) = dP(2)-comp_dom.ref_points[3*i](2);
            comp_dom.edges_tangents[3*i] = V1.X();
            comp_dom.edges_tangents[3*i+1] = V1.Y();
            comp_dom.edges_tangents[3*i+2] = V1.Z();
            //cout<<i<<" (point) "<<comp_dom.support_points[i]<<" vs "<<dP<<endl;
            //cout<<i<<" (edges_tangents) "<<comp_dom.edges_tangents(3*i)<<","<<comp_dom.edges_tangents(3*i+1)<<","<<comp_dom.edges_tangents(3*i+2)<<endl;
            }              
            }

     // this cycle hooks the boat and far field double nodes
     // to their water twins that have been moved
     for (unsigned int i=0; i<comp_dom.vector_dh.n_dofs(); ++i)
         {
         if ( (comp_dom.vector_flags[i] & water) &&
              (comp_dom.vector_flags[i] & edge)  )
            {
            std::set<unsigned int> duplicates = comp_dom.vector_double_nodes_set[i];
            duplicates.erase(i); 
            for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
                {//cout<<*pos<<"("<<i<<") "<<comp_dom.map_points(i)<<" vs "<<comp_dom.map_points(*pos)<<"   diff "<<comp_dom.map_points(i)-comp_dom.map_points(*pos)<<endl;
                working_map_points(*pos) = working_map_points(i);
                }
            }
         }

     //we enforce constraint on the new geometry assigned by the DAE solver
     vector_constraints.distribute(working_map_points);

     // mesh is ready, now we copy on map points the new nodes displacements
     comp_dom.map_points = working_map_points;
     // and update the support points
     comp_dom.update_support_points();
     // and compute the node normals (needed by bem b.c.)
     comp_dom.compute_normals_at_nodes(comp_dom.map_points);

     for (unsigned int i=0; i<comp_dom.vector_dh.n_dofs(); ++i)
         {
         y(i) = comp_dom.map_points(i)-comp_dom.rigid_motion_map_points(i);
         }    

//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// here we prepare the initial_map_points to be used in case of transom stern
  if (!comp_dom.no_boat && comp_dom.boat_model.is_transom)
     {
     comp_dom.update_support_points();
     
     double transom_draft = fabs(comp_dom.boat_model.CurrentPointCenterTransom(2));
     double transom_aspect_ratio = (fabs(comp_dom.boat_model.CurrentPointLeftTransom(1))+
                                     fabs(comp_dom.boat_model.CurrentPointRightTransom(1)))/transom_draft;

     wind.set_time(100000000.0);
     Vector<double> instantWindValueTinf(dim);
     Point<dim> zero(0,0,0);
     wind.vector_value(zero,instantWindValueTinf);
     Point<dim> VinfTinf;
     for (unsigned int i = 0; i < dim; i++)
       VinfTinf(i) = instantWindValueTinf(i);
     wind.set_time(initial_time);

     double FrT = sqrt(VinfTinf*VinfTinf)/sqrt(9.81*transom_draft);
     double ReT = sqrt(9.81*pow(transom_draft,3.0))/1.307e-6;
     double eta_dry = fmin(0.05*pow(FrT,2.834)*pow(transom_aspect_ratio,0.1352)*pow(ReT,0.01338),1.0);
     double lh = 0.0;
     //if (eta_dry < 1.0) 
        lh = 5.0; 
        //lh = 0.1135*pow(FrT,3.025)*pow(transom_aspect_ratio,0.4603)*pow(ReT,-0.1514);
        //lh = 0.3265*pow(FrT,3.0) - 1.7216*pow(FrT,2.0) + 2.7593*FrT;
     cout<<FrT<<" "<<transom_aspect_ratio<<" "<<ReT<<endl;
     cout<<"****eta_dry: "<<eta_dry<<endl;

     for(unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
        {
        comp_dom.initial_map_points(3*i+2) = 0.0;
        }

     for(unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
        {
        if ( (comp_dom.support_points[i](1) < comp_dom.boat_model.CurrentPointRightTransom(1)) &&
             (comp_dom.support_points[i](1) >= 0.0) &&
             (comp_dom.support_points[i](0) > comp_dom.boat_model.CurrentPointCenterTransom(0)-fabs(comp_dom.boat_model.CurrentPointCenterTransom(2)) ) &&
             (comp_dom.support_points[3*i](0) < comp_dom.boat_model.CurrentPointCenterTransom(0)+lh*fabs(comp_dom.boat_model.CurrentPointCenterTransom(2)) )   &&
             (comp_dom.support_points[i](2) > 2.0*comp_dom.boat_model.CurrentPointCenterTransom(2))    )
             {
             Point <3> dP0 = comp_dom.support_points[i];
             Point <3> dP; 
         				   //this is the vertical plane
             Handle(Geom_Plane) vertPlane = new Geom_Plane(0.,1.,0.,-dP0(1));
             vertPlane->Transform(comp_dom.boat_model.current_loc.Inverted());
             Handle(Geom_Curve) curve = comp_dom.boat_model.right_transom_bspline;

             GeomAPI_IntCS Intersector(curve, vertPlane);
             int npoints = Intersector.NbPoints();
             AssertThrow((npoints != 0), ExcMessage("Transom curve is not intersecting with vertical plane!"));
             //cout<<"Number of intersections: "<<npoints<<endl;
             double minDistance=1e7;
             for (int j=0; j<npoints;++j)
                 {
                 gp_Pnt int_point = Intersector.Point(j+1);
                 int_point.Transform(comp_dom.boat_model.current_loc);
                 Point<3> inters = Pnt(int_point);

                 if (dP0.distance(inters) < minDistance)
                    {
                    minDistance = dP0.distance(inters);
                    dP = inters;
                    }
                 }

             if ( (comp_dom.support_points[i](0) > dP(0)+comp_dom.min_diameter/20.0) &&
                  (comp_dom.support_points[i](0) < dP(0)+lh*fabs(dP(2))+comp_dom.min_diameter/20.0) )
                {
                double mean_curvature;
                Point<3> normal;
                Point<3> projection;
                comp_dom.boat_model.boat_surface_right->normal_projection_and_diff_forms(projection,
                                                                                        normal,
                                                                                        mean_curvature,
			                                                                dP);
                AssertThrow((dP.distance(projection) < 1e-4*comp_dom.boat_model.boatWetLength), ExcMessage("Normal projection for surface normal evaluation went wrong!"));
                double transom_slope = -normal(0)/normal(2);
                double a = -transom_slope/(lh*fabs(dP(2))) - 1/dP(2)/lh/lh;
                double x = comp_dom.support_points[i](0)-dP(0);
                comp_dom.initial_map_points(3*i+2) = a*x*x + transom_slope*x + dP(2);
                }
             }
             else if ( (comp_dom.support_points[i](1) > comp_dom.boat_model.CurrentPointLeftTransom(1)) &&
                       (comp_dom.support_points[i](1) < 0.0) &&
                       (comp_dom.support_points[i](0) > comp_dom.boat_model.CurrentPointCenterTransom(0)-fabs(comp_dom.boat_model.CurrentPointCenterTransom(2))) &&
                       (comp_dom.support_points[i](0) < comp_dom.boat_model.CurrentPointCenterTransom(0)+lh*fabs(comp_dom.boat_model.CurrentPointCenterTransom(2)))  &&
                       (comp_dom.support_points[i](2) > 2.0*comp_dom.boat_model.CurrentPointCenterTransom(2))   )
             {
             Point <3> dP0 = comp_dom.support_points[i];
             Point <3> dP; 
         				   //this is the vertical plane
             Handle(Geom_Plane) vertPlane = new Geom_Plane(0.,1.,0.,-dP0(1));
             vertPlane->Transform(comp_dom.boat_model.current_loc.Inverted());
             Handle(Geom_Curve) curve = comp_dom.boat_model.left_transom_bspline;

             GeomAPI_IntCS Intersector(curve, vertPlane);
             int npoints = Intersector.NbPoints();
             AssertThrow((npoints != 0), ExcMessage("Transom curve is not intersecting with vertical plane!"));
             //cout<<"Number of intersections: "<<npoints<<endl;
             double minDistance=1e7;
             for (int j=0; j<npoints;++j)
                 {
                 gp_Pnt int_point = Intersector.Point(j+1);
                 int_point.Transform(comp_dom.boat_model.current_loc);
                 Point<3> inters = Pnt(int_point);

                 if (dP0.distance(inters) < minDistance)
                    {
                    minDistance = dP0.distance(inters);
                    dP = inters;
                    }
                 }
             if ( (comp_dom.support_points[i](0) > dP(0)+comp_dom.min_diameter/20.0) &&
                  (comp_dom.support_points[i](0) < dP(0)+lh*fabs(dP(2))+comp_dom.min_diameter/20.0) )
                {
                double mean_curvature;
                Point<3> normal;
                Point<3> projection;
                comp_dom.boat_model.boat_surface_left->normal_projection_and_diff_forms(projection,
                                                                                        normal,
                                                                                        mean_curvature,
			                                                                dP);
                AssertThrow((dP.distance(projection) < 1e-4*comp_dom.boat_model.boatWetLength), ExcMessage("Normal projection for surface normal evaluation went wrong!"));
                double transom_slope = -normal(0)/normal(2);
                double a = -transom_slope/(lh*fabs(dP(2))) - 1/dP(2)/lh/lh;
                double x = comp_dom.support_points[i](0)-dP(0);
                comp_dom.initial_map_points(3*i+2) = a*x*x + transom_slope*x + dP(2);
                }
             }
        }
  comp_dom.vector_constraints.distribute(comp_dom.initial_map_points);
  }






//////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////
  // here we take care of the free surface boundary condition (differential)
  // part of the variables, and the boat neumann condition variables
////////////////////////////////////////////////////////////////////// 
   

// we first need to fix the nodes horiziontal velocities on the boat (they must be parallel to the boat
// and respect the differential equation)

  Vector<double> complete_potential_gradients(comp_dom.vector_dh.n_dofs()); 
  compute_potential_gradients(complete_potential_gradients,phi,dphi_dn);

  for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
      {
      if ( (comp_dom.flags[i] & water) &&
           (comp_dom.flags[i] & near_boat) &&
           !(comp_dom.flags[i] & transom_on_water) &&
           (comp_dom.moving_point_ids[3] != i) &&
           (comp_dom.moving_point_ids[4] != i) &&
           (comp_dom.moving_point_ids[5] != i) &&
           (comp_dom.moving_point_ids[6] != i) )
         {
         Point<3> phi_gradient(complete_potential_gradients(3*i),complete_potential_gradients(3*i+1),complete_potential_gradients(3*i+2));
         Point<3> eta_gradient(-comp_dom.node_normals[i](0)/comp_dom.node_normals[i](2),-comp_dom.node_normals[i](1)/comp_dom.node_normals[i](2),0.0);
         double eta_dot = (phi_gradient(2)-eta_gradient*(Vinf+phi_gradient))/(1-eta_gradient(1)*comp_dom.iges_normals[i](2)/comp_dom.iges_normals[i](1));
         //yp(3*i+1) = -eta_dot*comp_dom.iges_normals[i](2)/comp_dom.iges_normals[i](1);
         std::set <unsigned int> duplicates = comp_dom.vector_double_nodes_set[3*i+1];
         for (std::set <unsigned int>::iterator pos = duplicates.begin(); pos != duplicates.end(); ++pos)
             {
             yp(*pos) = -eta_dot*comp_dom.iges_normals[i](2)/comp_dom.iges_normals[i](1);
             }
         //cout<<"WL "<<i<<" "<<eta_dot<<" "<<yp(3*i+1)<<"   ("<<comp_dom.iges_normals[i]<<")"<<endl;
         }
      }

     //this takes care of the bow and stern nodes 
     if (!comp_dom.no_boat)
        for (unsigned int k=3; k<7; ++k)
        { 
        unsigned int i = comp_dom.moving_point_ids[k];
        Point<3> phi_gradient(complete_potential_gradients(3*i),complete_potential_gradients(3*i+1),complete_potential_gradients(3*i+2));
        Point<3> eta_gradient(-comp_dom.node_normals[i](0)/comp_dom.node_normals[i](2),-comp_dom.node_normals[i](1)/comp_dom.node_normals[i](2),0.0);
        Point<3> t(comp_dom.edges_tangents[3*i],comp_dom.edges_tangents[3*i+1],comp_dom.edges_tangents[3*i+2]);     
        double eta_dot = (phi_gradient(2)-eta_gradient*(Vinf+phi_gradient))/(1.0-t(0)/t(2)-t(1)/t(2));
        std::set <unsigned int> duplicates = comp_dom.vector_double_nodes_set[3*i];
        for (std::set <unsigned int>::iterator pos = duplicates.begin(); pos != duplicates.end(); ++pos)
            {
            yp(*pos) = eta_dot*t(0)/t(2);
            }
        duplicates = comp_dom.vector_double_nodes_set[3*i+1];
        for (std::set <unsigned int>::iterator pos = duplicates.begin(); pos != duplicates.end(); ++pos)
            {
            yp(*pos) = eta_dot*t(1)/t(2);
            }

        //cout<<"KT "<<i<<" "<<eta_dot<<" "<<yp(3*i)<<" "<<yp(3*i+1)<<"   (";
        //cout<<comp_dom.edges_tangents(3*i)<<" "<<comp_dom.edges_tangents(3*i+1)<<" "<<comp_dom.edges_tangents(3*i+2)<<")"<<endl;
        }              



// building reference cell
  Triangulation<2,3> ref_triangulation;

  std::vector<Point<3> > ref_vertices;
  std::vector<CellData<2> > ref_cells;
  SubCellData ref_subcelldata;

  ref_vertices.resize(4);
  ref_vertices[0](0)=-1.0; ref_vertices[0](1)=-1.0; ref_vertices[0](2)=0.0;
  ref_vertices[1](0)= 1.0; ref_vertices[1](1)=-1.0; ref_vertices[1](2)=0.0;
  ref_vertices[2](0)=-1.0; ref_vertices[2](1)= 1.0; ref_vertices[2](2)=0.0;
  ref_vertices[3](0)= 1.0; ref_vertices[3](1)= 1.0; ref_vertices[3](2)=0.0;

  ref_cells.resize(1);

  ref_cells[0].vertices[0]=0; ref_cells[0].vertices[1]=1; ref_cells[0].vertices[2]=3; ref_cells[0].vertices[3]=2;
  ref_cells[0].material_id = 1;

  GridTools::delete_unused_vertices (ref_vertices, ref_cells, ref_subcelldata);
  GridReordering<2,3>::reorder_cells (ref_cells);

  ref_triangulation.create_triangulation_compatibility(ref_vertices, ref_cells, ref_subcelldata );

  FE_Q<2,3> fe(1);
  DoFHandler<2,3> ref_dh(ref_triangulation);
  ref_dh.distribute_dofs(fe);

  FEValues<2,3> ref_fe_v(StaticMappingQ1<2,3>::mapping, fe, *comp_dom.quadrature,
   		         update_values | update_gradients |
		         update_cell_normal_vectors |
		         update_quadrature_points |
		         update_JxW_values);

  const unsigned int n_q_points = ref_fe_v.n_quadrature_points;
  const unsigned int  dofs_per_cell   = fe.dofs_per_cell;

  cell_it ref_cell = ref_dh.begin_active();
  ref_fe_v.reinit(ref_cell);  

///////////////////////////////////

// test: let's try assemble actual mass matrix

  DphiDt_sys_matrix = 0;
  DphiDt_sys_solution = 0;
  DphiDt_sys_solution_2 = 0;
  DphiDt_sys_solution_3 = 0;
  DphiDt_sys_rhs = 0;
  DphiDt_sys_rhs_2 = 0;
  DphiDt_sys_rhs_3 = 0;

//let's build the residual on the free surface cells (differential components)
  std::vector<double> eta_res(comp_dom.dh.n_dofs(),0.0);
  std::vector<double> phi_res(comp_dom.dh.n_dofs(),0.0);
  std::vector<double> dphi_dn_res(comp_dom.dh.n_dofs(),0.0);
  std::vector<double> x_smoothing_res(comp_dom.dh.n_dofs(),0.0);
  std::vector<double> y_smoothing_res(comp_dom.dh.n_dofs(),0.0);


  double g = 9.81;

  FullMatrix<double>   local_DphiDt_matrix (dofs_per_cell, dofs_per_cell);
  FullMatrix<double>   local_DphiDt_matrix_2 (dofs_per_cell, dofs_per_cell);
  Vector<double>       local_DphiDt_rhs (dofs_per_cell);
  Vector<double>       local_DphiDt_rhs_2 (dofs_per_cell);
  Vector<double>       local_DphiDt_rhs_3 (dofs_per_cell);

  cell_it
  cell = comp_dom.dh.begin_active(),
  endc = comp_dom.dh.end();

  std::vector<fad_double> coors(3*dofs_per_cell,0.0);
  std::vector<fad_double> phis(dofs_per_cell,0.0);
  std::vector<fad_double> dphi_dns(dofs_per_cell,0.0);
  std::vector<fad_double> coors_dot(3*dofs_per_cell,0.0);
  std::vector<fad_double> phis_dot(dofs_per_cell,0.0);
  std::vector<fad_double> dphi_dns_dot(dofs_per_cell,0.0);
  std::vector<fad_double> x_displs(dofs_per_cell,0.0);
  std::vector<fad_double> y_displs(dofs_per_cell,0.0);
  
  std::vector<fad_double> loc_eta_res(dofs_per_cell,0.0);
  std::vector<fad_double> loc_phi_res(dofs_per_cell,0.0);
  std::vector<fad_double> loc_dphi_dn_res(dofs_per_cell,0.0);
  std::vector<fad_double> loc_x_smooth_res(dofs_per_cell,0.0);
  std::vector<fad_double> loc_y_smooth_res(dofs_per_cell,0.0);
 
  std::vector< std::vector<fad_double> > loc_stiffness_matrix(dofs_per_cell,std::vector<fad_double>(dofs_per_cell));
  std::vector< std::vector<fad_double> > loc_mass_matrix(dofs_per_cell,std::vector<fad_double>(dofs_per_cell));
  std::vector< std::vector<fad_double> > loc_supg_mass_matrix(dofs_per_cell,std::vector<fad_double>(dofs_per_cell));
 

  std::vector<unsigned int> local_dof_indices(dofs_per_cell);
   

  //for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
  //       cout<<i<<"  VelsAgain: "<<" "<<nodes_velocities(3*i)<<" "<<nodes_velocities(3*i+1)<<" "<<nodes_velocities(3*i+2)<<endl;

  for (; cell!=endc; ++cell)
      {
      //if (cell->material_id() == comp_dom.free_sur_ID1 ||
      //    cell->material_id() == comp_dom.free_sur_ID2 ||
      //    cell->material_id() == comp_dom.free_sur_ID3 ||
      //    cell->material_id() == comp_dom.wall_sur_ID1 ||
      //    cell->material_id() == comp_dom.wall_sur_ID2 ||
      //    cell->material_id() == comp_dom.wall_sur_ID3 )
         {
         //std::cout<<std::endl;
         //std::cout<<"Cell: "<<cell<<"  Center: "<<cell->center()<<" --- "<<dofs_per_cell<<std::endl;

         local_DphiDt_matrix = 0;
         local_DphiDt_matrix_2 = 0;
         local_DphiDt_rhs = 0;
         local_DphiDt_rhs_2 = 0;
         local_DphiDt_rhs_3 = 0; 

         cell->get_dof_indices(local_dof_indices);
         for (unsigned int i=0; i<dofs_per_cell; ++i)
             {
             loc_eta_res[i] = 0;
             loc_phi_res[i] = 0;
             loc_dphi_dn_res[i] = 0;
             loc_x_smooth_res[i] = 0;
             loc_y_smooth_res[i] = 0;
             for (unsigned int j=0; j<dofs_per_cell; ++j)
                 {
                 loc_mass_matrix[i][j] = 0;
                 loc_supg_mass_matrix[i][j] = 0;
                 loc_stiffness_matrix[i][j] = 0;
                 }
             //cout<<local_dof_indices[i]<<" ";
             }

         for (unsigned int i=0; i<dofs_per_cell; ++i)
             {
             for (unsigned int j=0; j<3; ++j)
                 {
                 //cout<<3*local_dof_indices[i]+j<<endl;
                 coors[3*i+j] = comp_dom.support_points[local_dof_indices[i]](j);
                 coors[3*i+j].diff(3*i+j,10*dofs_per_cell);
                 coors_dot[3*i+j] = nodes_velocities(3*local_dof_indices[i]+j);
                 coors_dot[3*i+j].diff(3*i+j+5*dofs_per_cell,10*dofs_per_cell);
                 //std::cout<<i<<"-------> "<<coors_dot[3*i+j].val()<<" vs "<<nodes_velocities(3*local_dof_indices[i]+j)<<endl;
                 }
             //std::cout<<i<<"-------> "<<coors_dot[3*i].val()<<" "<<coors_dot[3*i+1].val()<<" "<<coors_dot[3*i+2].val()<<" "<<endl;
             }
         for (unsigned int i=0; i<dofs_per_cell; ++i)
             {
             phis[i] = phi(local_dof_indices[i]);
             phis[i].diff(i+3*dofs_per_cell,10*dofs_per_cell);
             phis_dot[i] = phi_time_derivs(local_dof_indices[i]);
             phis_dot[i].diff(i+8*dofs_per_cell,10*dofs_per_cell);
             dphi_dns[i] = dphi_dn(local_dof_indices[i]);
             dphi_dns[i].diff(i+4*dofs_per_cell,10*dofs_per_cell);
             dphi_dns_dot[i] = dphi_dn_time_derivs(local_dof_indices[i]);
             dphi_dns_dot[i].diff(i+9*dofs_per_cell,10*dofs_per_cell);
             //std::cout<<i<<"--> "<<local_dof_indices[i]<<"--------->"<<bem_phi(local_dof_indices[i])<<"  "<<bem_dphi_dn(local_dof_indices[i])<<endl;
             }

         // computation of displacements
         for (unsigned int i=0; i<dofs_per_cell; ++i)
             {
             x_displs[i] = coors[3*i] - comp_dom.ref_points[3*local_dof_indices[i]](0);
             y_displs[i] = coors[3*i+1] - comp_dom.ref_points[3*local_dof_indices[i]](1);
             //cout<<i<<" "<<x_displs[i].val()<<" "<<y_displs[i].val()<<endl;
             //cout<<i<<" "<<coors[3*i].val()<<" "<<comp_dom.ref_points[3*local_dof_indices[i]](0)<<" "<<comp_dom.support_points[local_dof_indices[i]](0)<<endl;
             //cout<<i<<" "<<coors[3*i+1].val()<<" "<<comp_dom.ref_points[3*local_dof_indices[i]](1)<<" "<<comp_dom.support_points[local_dof_indices[i]](1)<<endl;
             }
         // computation of cell center
         Point<3,fad_double> center(0.0,0.0,0.0);
         for (unsigned int i=0; i<dofs_per_cell; ++i)
             {
             center += (Point<3,fad_double>(coors[3*i],coors[3*i+1],coors[3*i+2]))/dofs_per_cell;
             }

         std::vector<fad_double> eta_dot_rhs_fun(n_q_points);
         std::vector<fad_double> phi_dot_rhs_fun(n_q_points);
         std::vector< Point<3,fad_double> > fluid_vel(n_q_points);
         std::vector<fad_double> q_JxW(n_q_points);
   
   
         for (unsigned int q=0; q<n_q_points; ++q)
             {
             Point<3,fad_double> q_point(0.0,0.0,0.0);
             Point<3,fad_double> u_deriv_pos(0.0,0.0,0.0);
             Point<3,fad_double> v_deriv_pos(0.0,0.0,0.0);
             fad_double u_deriv_phi = 0;
             fad_double v_deriv_phi = 0;
             fad_double q_dphi_dn = 0;
             fad_double q_x_dot = 0;
             fad_double q_y_dot = 0;
             fad_double q_z_dot = 0;
             fad_double q_eta = 0;
             for (unsigned int i=0; i<dofs_per_cell; ++i)
                 {
                 unsigned int index = local_dof_indices[i];
                 q_point += fad_double(ref_fe_v.shape_value(i,q))*Point<3,fad_double>(coors[3*i],coors[3*i+1],coors[3*i+2]);
                 u_deriv_pos += fad_double(ref_fe_v.shape_grad(i,q)[0])*Point<3,fad_double>(coors[3*i],coors[3*i+1],coors[3*i+2]);
                 v_deriv_pos += fad_double(ref_fe_v.shape_grad(i,q)[1])*Point<3,fad_double>(coors[3*i],coors[3*i+1],coors[3*i+2]);
                 u_deriv_phi += fad_double(ref_fe_v.shape_grad(i,q)[0])*phis[i];
                 v_deriv_phi += fad_double(ref_fe_v.shape_grad(i,q)[1])*phis[i];
                 q_dphi_dn += fad_double(ref_fe_v.shape_value(i,q))*dphi_dns[i];
                 q_x_dot += fad_double(ref_fe_v.shape_value(i,q))*coors_dot[3*i];
                 q_y_dot += fad_double(ref_fe_v.shape_value(i,q))*coors_dot[3*i+1];
                 q_z_dot += fad_double(ref_fe_v.shape_value(i,q))*coors_dot[3*i+2];
                 q_eta +=  fad_double(ref_fe_v.shape_value(i,q))*coors[3*i+2];
                 //std::cout<<i<<"-------> "<<u_deriv_pos<<" "<<v_deriv_pos<<" "<<u_deriv_phi<<" "<<v_deriv_phi<<endl;
                 //std::cout<<i<<"-------> "<<coors[3*i].val()<<" "<<coors[3*i+1]<<" "<<coors[3*i+2]<<" "<<endl;
                 }
             Point<3,fad_double> q_normal(u_deriv_pos(1)*v_deriv_pos(2)-u_deriv_pos(2)*v_deriv_pos(1),
                                          u_deriv_pos(2)*v_deriv_pos(0)-u_deriv_pos(0)*v_deriv_pos(2),
                                          u_deriv_pos(0)*v_deriv_pos(1)-u_deriv_pos(1)*v_deriv_pos(0));
             //std::cout<<"q_normal="<<q_normal<<std::endl;
             //std::cout<<"q_y_dot="<<q_y_dot<<std::endl;
             fad_double q_jac_det = q_normal.norm();
             q_normal/=q_jac_det;
             fad_double a = 1.0/((u_deriv_pos*u_deriv_pos)*(v_deriv_pos*v_deriv_pos)-(u_deriv_pos*v_deriv_pos)*(u_deriv_pos*v_deriv_pos));
             fad_double d11 = a*(u_deriv_pos(0)*v_deriv_pos*v_deriv_pos-v_deriv_pos(0)*u_deriv_pos*v_deriv_pos);
             fad_double d21 = a*(u_deriv_pos(1)*v_deriv_pos*v_deriv_pos-v_deriv_pos(1)*u_deriv_pos*v_deriv_pos);
             fad_double d31 = a*(u_deriv_pos(2)*v_deriv_pos*v_deriv_pos-v_deriv_pos(2)*u_deriv_pos*v_deriv_pos);
             fad_double d12 = a*(v_deriv_pos(0)*u_deriv_pos*u_deriv_pos-u_deriv_pos(0)*u_deriv_pos*v_deriv_pos);
             fad_double d22 = a*(v_deriv_pos(1)*u_deriv_pos*u_deriv_pos-u_deriv_pos(1)*u_deriv_pos*v_deriv_pos);
             fad_double d32 = a*(v_deriv_pos(2)*u_deriv_pos*u_deriv_pos-u_deriv_pos(2)*u_deriv_pos*v_deriv_pos);
             Point<3,fad_double> phi_surf_grad(d11*u_deriv_phi+d12*v_deriv_phi,
                                               d21*u_deriv_phi+d22*v_deriv_phi,
                                               d31*u_deriv_phi+d32*v_deriv_phi);
             Point<3,fad_double> phi_surf_grad_corrected(phi_surf_grad(0) - phi_surf_grad(2)*q_normal(0)/q_normal(2),
                                                         phi_surf_grad(1) - phi_surf_grad(2)*q_normal(1)/q_normal(2),
                                                         0.0);
             Point<3,fad_double> phi_grad = phi_surf_grad + q_normal*q_dphi_dn;

             //std::cout<<"q_point="<<q_point<<"   q_normal="<<q_normal<<"   q_dphi_dn="<<q_dphi_dn<<std::endl;
             //cout<<q<<" phi_grad("<<phi_grad<<")  phi_surf_grad("<<phi_surf_grad<<")"<<endl;
             Point<3,fad_double> eta_grad(-q_normal(0)/q_normal(2),-q_normal(1)/q_normal(2),0.0);
             Point<3,fad_double> q_nodes_vel(q_x_dot,q_y_dot,q_z_dot);
             fluid_vel[q] = Point<3,fad_double>(fad_double(Vinf(0)),fad_double(Vinf(1)),fad_double(Vinf(2))) + phi_grad;
             fad_double fluid_vel_norm = fluid_vel[q].norm();
             if (fluid_vel_norm < 1e-3)
                fluid_vel_norm = -8.0e+05*pow(fluid_vel_norm,3.0) + 1.7e+03*pow(fluid_vel_norm,2.0) + 0.0001;
             fad_double cell_diameter;
             for (unsigned int i=0; i<dofs_per_cell; ++i)
                 {
                 cell_diameter += pow(fluid_vel[q]*(Point<3,fad_double>(coors[3*i],coors[3*i+1],coors[3*i+2])-center),2.0)/dofs_per_cell;
                 }
             cell_diameter = sqrt(cell_diameter)*2;

             eta_dot_rhs_fun[q] = phi_grad(2) + eta_grad*(q_nodes_vel-fluid_vel[q]);
             phi_dot_rhs_fun[q] = phi_grad*phi_grad/2 - g*q_eta + phi_surf_grad_corrected*(q_nodes_vel-fluid_vel[q]);
             q_JxW[q] = q_jac_det*ref_fe_v.JxW(q);

             //cout<<q<<" fvel("<<fluid_vel[q].val()<<")  fvel_norm="<<fluid_vel_norm<<"   q_JxW="<<q_JxW[q]<<endl;
             //cout<<q<<" erhs("<<eta_dot_rhs_fun[q]<<")  prhs("<<phi_dot_rhs_fun[q]<<")"<<endl;
             //cout<<q<<" phi_grad("<<phi_grad<<")  phi_surf_grad("<<phi_surf_grad<<")"<<endl;
            // cout<<q<<"   "<<phi_dot_rhs_fun[q].val()<<endl;//" "<<phi_dot_rhs_fun[q].val()<<endl;
             //cout<<q<<"   "<<phi_grad(0).val()<<" "<<phi_grad(1).val()<<" "<<phi_grad(2).val()<<endl;
             
             if (cell->material_id() == comp_dom.free_sur_ID1 ||
                 cell->material_id() == comp_dom.free_sur_ID2 ||
                 cell->material_id() == comp_dom.free_sur_ID3 )
                {
                for (unsigned int i=0;i<dofs_per_cell;++i)
                    {
                    Point<3,fad_double> N_i_surf_grad(d11*ref_fe_v.shape_grad(i,q)[0]+d12*ref_fe_v.shape_grad(i,q)[1],
                                                     d21*ref_fe_v.shape_grad(i,q)[0]+d22*ref_fe_v.shape_grad(i,q)[1],
                                                     d31*ref_fe_v.shape_grad(i,q)[0]+d32*ref_fe_v.shape_grad(i,q)[1]);
                    fad_double N_i_supg = fad_double(ref_fe_v.shape_value(i,q)) +
                                          N_i_surf_grad*fluid_vel[q]/fluid_vel_norm*cell_diameter/sqrt(2);
                    loc_eta_res[i] -= eta_dot_rhs_fun[q]*N_i_supg*q_JxW[q];
                    loc_phi_res[i] -= phi_dot_rhs_fun[q]*N_i_supg*q_JxW[q];
                    loc_x_smooth_res[i] -= blend_factor*0*fad_double(ref_fe_v.shape_value(i,q))*q_JxW[q];
                    loc_y_smooth_res[i] -= blend_factor*0*fad_double(ref_fe_v.shape_value(i,q))*q_JxW[q];
                    local_DphiDt_rhs(i) += (phi_dot_rhs_fun[q]*N_i_supg*q_JxW[q]).val();
                    local_DphiDt_rhs_2(i) += (eta_dot_rhs_fun[q]*N_i_supg*q_JxW[q]).val();
               //     cout<<q<<"   "<<endl; //<<local_DphiDt_rhs(i)<<" "<<local_DphiDt_rhs_2(i)<<endl;
               //     if (q_point(0).val()>36 && q_point(0).val()<70 && abs(q_point(1).val())< 6 && abs(q_point(2).val())< 1)             
               // {
               // std::cout<<"q_point=("<<q_point(0).val()<<","<<q_point(1).val()<<","<<q_point(2).val()<<")  q_dphi_dn="<<q_dphi_dn.val()<<endl;
               // std::cout<<"phi_grad=("<<phi_grad(0).val()<<","<<phi_grad(1).val()<<","<<phi_grad(2).val()<<")"<<endl;
               // std::cout<<"fluid_vel=("<<fluid_vel[q](0).val()<<","<<fluid_vel[q](1).val()<<","<<fluid_vel[q](2).val()<<")"<<endl;
               // std::cout<<"q_nodes_vel=("<<q_nodes_vel(0).val()<<","<<q_nodes_vel(1).val()<<","<<q_nodes_vel(2).val()<<")"<<endl;
               // std::cout<<"phi_surf_grad_corrected=("<<phi_surf_grad_corrected(0).val()<<","<<phi_surf_grad_corrected(1).val()<<","<<phi_surf_grad_corrected(2).val()<<")"<<endl;
               // cout<<q<<" erhs("<<eta_dot_rhs_fun[q].val()<<")  prhs("<<phi_dot_rhs_fun[q].val()<<")"<<endl;
              //  std::cout<<"local_DphiDt_rhs_2(i)="<<local_DphiDt_rhs_2(i)<<"local_DphiDt_rhs_2(i) ="<<local_DphiDt_rhs_2(i)<<")"<<endl;
              //  }
                
                    //cout<<q<<"  "<<i<<"   "<<phi_grad(2)<<"    "<<eta_grad<<"    "<<q_nodes_vel-fluid_vel[q]<<endl;
                    //cout<<q<<"  "<<i<<"   "<<N_i_surf_grad<<"    "<<fluid_vel[q]/fluid_vel_norm<<"   "<<cell_diameter/sqrt(2)<<endl;
                    //cout<<q<<"  "<<i<<"   "<<N_i_supg.val()<<"   "<<phi_dot_rhs_fun[q].val()<<"   "<<q_JxW[q].val()<<endl;
                    for (unsigned int j=0;j<dofs_per_cell;++j)
                        {
                        //loc_eta_res[i] += fad_double(ref_fe_v.shape_value(j,q))*coors_dot[3*j+2]*N_i_supg*q_JxW[q];
                        //loc_phi_res[i] += fad_double(ref_fe_v.shape_value(j,q))*phis_dot[j]*N_i_supg*q_JxW[q];
                        //local_DphiDt_matrix.add(i,j,ref_fe_v.shape_value(j,q)*(N_i_supg*q_JxW[q]).val());
                        loc_supg_mass_matrix[i][j] += fad_double(ref_fe_v.shape_value(j,q))*N_i_supg*q_JxW[q];
                        Point<3,fad_double> N_j_surf_grad(d11*ref_fe_v.shape_grad(j,q)[0]+d12*ref_fe_v.shape_grad(j,q)[1],
                                                          d21*ref_fe_v.shape_grad(j,q)[0]+d22*ref_fe_v.shape_grad(j,q)[1],
                                                          d31*ref_fe_v.shape_grad(j,q)[0]+d32*ref_fe_v.shape_grad(j,q)[1]);
                        loc_stiffness_matrix[i][j] += N_i_surf_grad*N_j_surf_grad*q_JxW[q];
                        }
                    //if (fmax(abs(loc_eta_res[i].val()),abs(loc_phi_res[i].val()))>1e-6)
                    //   cout<<q<<"  "<<i<<"   "<<loc_eta_res[i].val()<<"("<<coors_dot[3*i+2].val()<<")  "<<loc_phi_res[i].val()<<"("<<phis_dot[i].val()<<")  "<<endl;   
                    }
                
                }

             if (cell->material_id() == comp_dom.wall_sur_ID1 ||
                 cell->material_id() == comp_dom.wall_sur_ID2 ||
                 cell->material_id() == comp_dom.wall_sur_ID3 )
                {
                for (unsigned int i=0;i<dofs_per_cell;++i)
                    {
                    loc_dphi_dn_res[i] -= -(q_normal*Point<3,fad_double>(fad_double(Vinf(0)),fad_double(Vinf(1)),fad_double(Vinf(2))))*
                                           fad_double(ref_fe_v.shape_value(i,q))*q_JxW[q];
                    //cout<<q<<"  "<<i<<"   "<<-(q_normal*Point<3,fad_double>(fad_double(Vinf(0)),fad_double(Vinf(1)),fad_double(Vinf(2)))).val()<<"    "<<cell->center()<<"    "<<endl;
                    //cout<<q<<"  "<<i<<"   "<<N_i_surf_grad<<"    "<<fluid_vel[q]/fluid_vel_norm<<"   "<<cell_diameter/sqrt(2)<<endl;
                    //cout<<q<<"  "<<i<<"   "<<N_i_supg.val()<<"   "<<phi_dot_res_fun[q].val()<<"   "<<q_JxW[q].val()<<endl;
                    local_DphiDt_rhs_3(i) += (-(q_normal*Point<3,fad_double>(fad_double(Vinf(0)),fad_double(Vinf(1)),fad_double(Vinf(2))))*
                                           fad_double(ref_fe_v.shape_value(i,q))*q_JxW[q]).val();
                    //cout<<"**** "<<loc_dphi_dn_res[i].val()<<" "<<local_DphiDt_rhs_3(i)<<" "<<loc_dphi_dn_res[i].val()+local_DphiDt_rhs_3(i)<<endl;
                    for (unsigned int j=0;j<dofs_per_cell;++j)
                        {
                        //loc_dphi_dn_res[i] += fad_double(ref_fe_v.shape_value(i,q))*fad_double(ref_fe_v.shape_value(j,q))*dphi_dns[j]*q_JxW[q];
                        loc_mass_matrix[i][j] += fad_double(ref_fe_v.shape_value(j,q))*fad_double(ref_fe_v.shape_value(i,q))*q_JxW[q]; 
                        }
                    //if (abs(loc_dphi_dn_res[i].val())>1e-7)
                    //   cout<<q<<"  "<<i<<"   "<<loc_dphi_dn_res[i].val()<<endl;   
                    }
                }

             if (cell->material_id() != comp_dom.wall_sur_ID1 &&
                 cell->material_id() != comp_dom.wall_sur_ID2 &&
                 cell->material_id() != comp_dom.wall_sur_ID3 &&
                 cell->material_id() != comp_dom.free_sur_ID1 &&
                 cell->material_id() != comp_dom.free_sur_ID2 &&
                 cell->material_id() != comp_dom.free_sur_ID3)
                {
                for (unsigned int i=0;i<dofs_per_cell;++i)
                    {
                    loc_dphi_dn_res[i] -= 0;
                    //cout<<q<<"  "<<i<<"   "<<-(q_normal*Point<3,fad_double>(fad_double(Vinf(0)),fad_double(Vinf(1)),fad_double(Vinf(2)))).val()<<"    "<<cell->center()<<"    "<<endl;
                    //cout<<q<<"  "<<i<<"   "<<N_i_surf_grad<<"    "<<fluid_vel[q]/fluid_vel_norm<<"   "<<cell_diameter/sqrt(2)<<endl;
                    //cout<<q<<"  "<<i<<"   "<<N_i_supg.val()<<"   "<<phi_dot_res_fun[q].val()<<"   "<<q_JxW[q].val()<<endl;
                    local_DphiDt_rhs_3(i) += 0;
                    //cout<<"**** "<<loc_dphi_dn_res[i].val()<<" "<<local_DphiDt_rhs_3(i)<<" "<<loc_dphi_dn_res[i].val()+local_DphiDt_rhs_3(i)<<endl;
                    for (unsigned int j=0;j<dofs_per_cell;++j)
                        {
                        //loc_dphi_dn_res[i] += fad_double(ref_fe_v.shape_value(i,q))*fad_double(ref_fe_v.shape_value(j,q))*dphi_dns[j]*q_JxW[q];
                        loc_mass_matrix[i][j] += fad_double(ref_fe_v.shape_value(j,q))*fad_double(ref_fe_v.shape_value(i,q))*q_JxW[q]; 
                        }
                    //if (abs(loc_dphi_dn_res[i].val())>1e-7)
                    //   cout<<q<<"  "<<i<<"   "<<loc_dphi_dn_res[i].val()<<endl;   
                    }
                }

             for (unsigned int i=0;i<dofs_per_cell;++i)
                 {
                 Point<3,fad_double> N_i_surf_grad(d11*ref_fe_v.shape_grad(i,q)[0]+d12*ref_fe_v.shape_grad(i,q)[1],
                                                   d21*ref_fe_v.shape_grad(i,q)[0]+d22*ref_fe_v.shape_grad(i,q)[1],
                                                   d31*ref_fe_v.shape_grad(i,q)[0]+d32*ref_fe_v.shape_grad(i,q)[1]);
                 fad_double N_i_supg = fad_double(ref_fe_v.shape_value(i,q)) +
                                       N_i_surf_grad*fluid_vel[q]/fluid_vel_norm*cell_diameter/sqrt(2);
                 for (unsigned int j=0;j<dofs_per_cell;++j)
                     {
                     local_DphiDt_matrix.add(i,j,ref_fe_v.shape_value(j,q)*(N_i_supg*q_JxW[q]).val());
                     local_DphiDt_matrix_2.add(i,j,ref_fe_v.shape_value(j,q)*ref_fe_v.shape_value(i,q)*(q_JxW[q]).val());
                     //loc_supg_mass_matrix[i][j] += fad_double(ref_fe_v.shape_value(j,q))*N_i_supg*q_JxW[q];
                     //loc_mass_matrix[i][j] += fad_double(ref_fe_v.shape_value(j,q))*fad_double(ref_fe_v.shape_value(i,q))*q_JxW[q]; 
                     }
                 } 
  	     }

         for (unsigned int i=0;i<dofs_per_cell;++i)
                 {
                 DphiDt_sys_rhs(local_dof_indices[i]) += local_DphiDt_rhs(i);
                 DphiDt_sys_rhs_2(local_dof_indices[i]) += local_DphiDt_rhs_2(i);
                 DphiDt_sys_rhs_3(local_dof_indices[i]) += local_DphiDt_rhs_3(i);
                 for (unsigned int j=0;j<dofs_per_cell;++j)
                     {
                     DphiDt_sys_matrix.add(local_dof_indices[i],local_dof_indices[j],local_DphiDt_matrix(i,j));
                     DphiDt_sys_matrix_2.add(local_dof_indices[i],local_dof_indices[j],local_DphiDt_matrix_2(i,j));
                     }
                 } 

         if (cell->material_id() == comp_dom.free_sur_ID1 ||
             cell->material_id() == comp_dom.free_sur_ID2 ||
             cell->material_id() == comp_dom.free_sur_ID3 )
             {
             for (unsigned int i=0;i<dofs_per_cell;++i)
                 {
                 for (unsigned int j=0;j<dofs_per_cell;++j)
                     {
                     loc_eta_res[i] += loc_supg_mass_matrix[i][j]*coors_dot[3*j+2];
                     loc_phi_res[i] += loc_supg_mass_matrix[i][j]*phis_dot[j];
                     loc_x_smooth_res[i] += loc_stiffness_matrix[i][j]*(x_displs[j]-(1-blend_factor)*comp_dom.old_map_points(3*local_dof_indices[j]));
                     loc_y_smooth_res[i] += loc_stiffness_matrix[i][j]*(y_displs[j]-(1-blend_factor)*comp_dom.old_map_points(3*local_dof_indices[j]+1));
                     }
                 if ( !constraints.is_constrained(local_dof_indices[i]) &&
                      !(comp_dom.flags[local_dof_indices[i]] & transom_on_water) )
                    {
                    unsigned int ii = local_dof_indices[i];
                    eta_res[ii] += loc_eta_res[i].val();
                    phi_res[ii] += loc_phi_res[i].val();
                                         
                    }
                 if ( !(constraints.is_constrained(local_dof_indices[i])) &&
                      !(comp_dom.flags[local_dof_indices[i]] & edge) )
                    {
                    unsigned int ii = local_dof_indices[i];
                    x_smoothing_res[ii] += loc_x_smooth_res[i].val();
                    y_smoothing_res[ii] += loc_y_smooth_res[i].val();
                                         
                    }
                 }
             }
         if (cell->material_id() != comp_dom.free_sur_ID1 &&
             cell->material_id() != comp_dom.free_sur_ID2 &&
             cell->material_id() != comp_dom.free_sur_ID3 )
             {
             for (unsigned int i=0;i<dofs_per_cell;++i)
                 {
                 unsigned int ii = local_dof_indices[i];
                 for (unsigned int j=0;j<dofs_per_cell;++j)
                     {
                     loc_dphi_dn_res[i] += loc_mass_matrix[i][j]*dphi_dns[j];
                     }
                 if (!constraints.is_constrained(ii))
                    {
                    dphi_dn_res[ii] += loc_dphi_dn_res[i].val();
                    }
                 }
             }
         }
      }
//for (unsigned int i=0; i<comp_dom.dh.n_dofs();++i)
//    cout<<i<<"--->  "<<eta_res[i]<<endl;;

     SparseMatrix<double> DphiDt_sys_matrix_copy;
     DphiDt_sys_matrix_copy.reinit(DphiDt_sparsity_pattern);
     DphiDt_sys_matrix_copy.copy_from(DphiDt_sys_matrix);

     //for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
     //    {
     //    cout<<"*** "<<i<<" "<<DphiDt_sys_rhs(i)<<" "<<DphiDt_sys_rhs_2(i)<<" "<<DphiDt_sys_rhs_3(i)<<endl;
     //    }

     constraints.condense(DphiDt_sys_matrix,DphiDt_sys_rhs);
     constraints.condense(DphiDt_sys_matrix_copy,DphiDt_sys_rhs_2);
     constraints.condense(DphiDt_sys_matrix_2,DphiDt_sys_rhs_3);
 

     SparseDirectUMFPACK DphiDt_direct;
     SparseDirectUMFPACK DphiDt_direct_copy;
     SparseDirectUMFPACK DphiDt_direct_2;
     
     DphiDt_direct.initialize(DphiDt_sys_matrix);
     DphiDt_direct_copy.initialize(DphiDt_sys_matrix_copy);
     DphiDt_direct_2.initialize(DphiDt_sys_matrix_2);

     DphiDt_direct.vmult(DphiDt_sys_solution, DphiDt_sys_rhs); // solving for phi_dot
     constraints.distribute(DphiDt_sys_solution);
     DphiDt_direct_copy.vmult(DphiDt_sys_solution_2, DphiDt_sys_rhs_2); // solving for eta_dot
     constraints.distribute(DphiDt_sys_solution_2);
     DphiDt_direct_2.vmult(DphiDt_sys_solution_3, DphiDt_sys_rhs_3); // solving for dphi_dn
     constraints.distribute(DphiDt_sys_solution_3);

     Vector<double> RES(DphiDt_sys_solution.size());
     DphiDt_sys_matrix.vmult(RES,DphiDt_sys_solution_2);
     RES*=-1.0;
     RES.add(DphiDt_sys_rhs_2);
     RES*=-1.0;

     //for (unsigned int i=0; i<comp_dom.dh.n_dofs();i++)
     //    {
     //    if (constraints.is_constrained(i) == 0)
     //       cout<<"eta_dot("<<i<<") "<<DphiDt_sys_solution_2(i)<<"   res("<<i<<") "<<RES(i)<<"  eta_res("<<i<<") "<<eta_res[i]<<endl;
     //    }

     for (unsigned int i=0; i<comp_dom.dh.n_dofs();i++)
         {
         if ( (comp_dom.flags[i] & water) &&
              !(comp_dom.flags[i] & transom_on_water) ) 
            {
            yp(3*i+2) = DphiDt_sys_solution_2(i);
     //       cout<<3*i+2<<" -> "<<yp(3*i+2)<<endl;
            yp(i+comp_dom.vector_dh.n_dofs()) = DphiDt_sys_solution(i);
     //       cout<<i+comp_dom.vector_dh.n_dofs()<<" -> "<<yp(i+comp_dom.vector_dh.n_dofs())<<endl;
            }
         else
            {
            y(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()) = DphiDt_sys_solution_3(i);
            }
         }


     RestartNonlinearProblemAlg rest_nonlin_prob_alg(*this,comp_dom,t,y,yp,jacobian_matrix);
     NewtonSolver restart_solver_alg(rest_nonlin_prob_alg);

     std::map<unsigned int,unsigned int> &map_alg = rest_nonlin_prob_alg.indices_map;
     Vector<double> restart_prob_solution_alg(rest_nonlin_prob_alg.n_dofs());
     for (std::map<unsigned int, unsigned int>::iterator it = map_alg.begin(); it != map_alg.end(); ++it)
         restart_prob_solution_alg(it->second) = y(it->first);
     restart_solver_alg.solve(restart_prob_solution_alg,5);

     //for (unsigned int i=0; i<restart_prob_solution.size(); ++i)
     //    cout<<i<<" "<<restart_prob_solution(i)<<endl;

     for (std::map<unsigned int, unsigned int>::iterator it = map_alg.begin(); it != map_alg.end(); ++it)
         {
         y(it->first) = restart_prob_solution_alg(it->second);
         //cout<<it->first<<" "<<yp(it->first)<<"("<<restart_prob_solution(it->second)<<")"<<endl;
         }
//////////////////////////////////////////////////////////////////////
  // here we take care of the bem part of the variables
////////////////////////////////////////////////////////////////////// 


     bem_phi = (const Vector<double> &)phi;
     constraints.distribute(bem_phi);

     bem_dphi_dn = (const Vector<double> &)dphi_dn;
     //constraints.distribute(bem_dphi_dn); 

     Vector<double> bem_bc(comp_dom.dh.n_dofs());
     for (unsigned int i=0; i<comp_dom.dh.n_dofs(); i++)
         {
         if (comp_dom.flags[i] & water)
            bem_bc(i) = bem_phi(i);
         else
            {
            if (comp_dom.flags[i] & boat) 
               bem_bc(i) = bem_dphi_dn(i);
            else
               bem_bc(i) = 0.0;
            }
         }

      // trying a fix for transom stern nodes
      for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
          {
          if ( comp_dom.flags[i] & transom_on_water )
             {
	     comp_dom.surface_nodes(i) = 0;
	     comp_dom.other_nodes(i) = 1;
             std::set<unsigned int> duplicates = comp_dom.double_nodes_set[i];
             duplicates.erase(i); 
             bem_bc(i) = 0;
             jacobian_matrix.add(i,i,-1.0);
             for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
                 {
                 bem_bc(i) += bem_dphi_dn(*pos)/duplicates.size();
                 }
             bem_dphi_dn(i) = bem_bc(i);
             }
          } 

       // this is to enforce constraints in a more strict way
       // pure normals might not be respecting constraints
       //constraints.distribute(bem_phi);
       //constraints.distribute(bem_dphi_dn);
       //constraints.distribute(bem_bc); 

        bem.solve(bem_phi, bem_dphi_dn, bem_bc);

       // this is to enforce constraints in a more strict way
       // the vector given back by bem has constraints
       // imposed up to the bem GMRES tolerance
       //constraints.distribute(bem_phi);
       //constraints.distribute(bem_dphi_dn);

       // finally this imposes that the potential on water side of 
       // transom line equals that on boat side
       for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
           {
           if (comp_dom.flags[i] & transom_on_water)
              {
              std::set<unsigned int> duplicates = comp_dom.double_nodes_set[i];
              duplicates.erase(i);
              bem_phi(i) = bem_phi(*duplicates.begin());
              }
           }


       // now bem_phi and bem_dphi_dn are correct, we copy them on phi and dphi_dn
       for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
           { 
           y(i+comp_dom.vector_dh.n_dofs()) = bem_phi(i);
           y(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()) = bem_dphi_dn(i);
           //if (comp_dom.flags[i] & water)
           //   {
           //   cout<<"bem_phi("<<i<<") "<<y(i+comp_dom.vector_dh.n_dofs())<<endl;
           //   cout<<"bem_dphi_dn("<<i<<") "<<y(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs())<<endl;
           //   }
           }

//////////////////////////////////////////////////////////////////////




  //residual(t,res,y,yp);
  //setup_jacobian_prec(t,y,yp,0.0);

  RestartNonlinearProblemDiff rest_nonlin_prob_diff(*this,comp_dom,t,y,yp,jacobian_dot_matrix);
  std::map<unsigned int,unsigned int> &map_diff = rest_nonlin_prob_diff.indices_map;

/*
  // these lines test the correctness of the jacobian for the
  // restart (reduced) nonlinear problem

  Vector<double> restart_prob_solution(rest_nonlin_prob_diff.n_dofs());
  Vector<double> restart_prob_residual(rest_nonlin_prob_diff.n_dofs());
  for (std::map<unsigned int, unsigned int>::iterator it = map_diff.begin(); it != map_diff.end(); ++it)
      restart_prob_solution(it->second) = yp(it->first);

     Vector<double> delta_y(rest_nonlin_prob_diff.n_dofs());
     Vector<double> delta_res(rest_nonlin_prob_diff.n_dofs());
     //delta_y.add(1e-8);
     //delta_y(974) = 1e-8;
     for (unsigned int i=0; i<rest_nonlin_prob_diff.n_dofs();++i)
         {
         double f = (double)rand()/RAND_MAX;
         delta_y(i) = -1e-8 + f * (2e-8);
         cout<<i<<" "<<delta_y(i)<<endl;
         }
     rest_nonlin_prob_diff.residual(restart_prob_residual,restart_prob_solution);
     rest_nonlin_prob_diff.jacobian(delta_res,restart_prob_solution,delta_y);
     restart_prob_solution.add(delta_y);
     //yp.add(delta_y);
     delta_res.add(restart_prob_residual);
     rest_nonlin_prob_diff.residual(restart_prob_residual,restart_prob_solution);
     cout<<"----------Test---------"<<endl;
     for (unsigned int i=0; i<rest_nonlin_prob_diff.n_dofs(); ++i)
         if (fabs(restart_prob_residual(i)-delta_res(i)) > 1e-10)
         cout<<i<<"  "<<delta_res(i)<<" vs "<<restart_prob_residual(i)<<"      err "<<restart_prob_residual(i)-delta_res(i)<<"   "<<sys_comp(i)<<endl;
     delta_res*=-1;
     delta_res.add(restart_prob_residual);
     cout<<"Absolute error norm: "<<delta_res.l2_norm()<<endl;
     cout<<"Relative error norm: "<<delta_res.l2_norm()/delta_y.l2_norm()<<endl;
     cout<<"----------Done---------"<<endl;
     for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
         cout<<i<<" "<<comp_dom.support_points[i]<<"     sn "<<comp_dom.surface_nodes(i)<<"    ic "<<constraints.is_constrained(i)<<endl;
     
//*/



  NewtonSolver restart_solver_diff(rest_nonlin_prob_diff);

  Vector<double> restart_prob_solution_diff(rest_nonlin_prob_diff.n_dofs());
  for (std::map<unsigned int, unsigned int>::iterator it = map_diff.begin(); it != map_diff.end(); ++it)
      restart_prob_solution_diff(it->second) = yp(it->first);
  restart_solver_diff.solve(restart_prob_solution_diff,8);

  //for (unsigned int i=0; i<restart_prob_solution.size(); ++i)
  //    cout<<i<<" "<<restart_prob_solution(i)<<endl;

  for (std::map<unsigned int, unsigned int>::iterator it = map_diff.begin(); it != map_diff.end(); ++it)
      {
      yp(it->first) = restart_prob_solution_diff(it->second);
      //cout<<it->first<<" "<<yp(it->first)<<"("<<restart_prob_solution_diff(it->second)<<")"<<endl;
      }


     //std::string filename1 = ( "post_restart_mesh.vtu" );
     //output_results(filename1, t, y, yp);


   
     // these lines test the jacobian of the DAE system
/* 
     Vector<double> delta_y(this->n_dofs());
     Vector<double> delta_res(this->n_dofs());

     for (unsigned int i=0; i<this->n_dofs();++i)
         {
         double f = (double)rand()/RAND_MAX;
         delta_y(i) = -1e-8 + f * (2e-8);
         cout<<i<<" "<<delta_y(i)<<endl;
         }
     Vector<double> delta_y_dot(delta_y);
     delta_y_dot*=00000001.0; 
     residual(t,res,y,yp);
     jacobian(t,delta_res,y,yp,delta_y,00000001.0);
     y.add(delta_y);
     yp.add(delta_y_dot);
     delta_res.add(res);
     residual(t,res,y,yp);
     cout<<"----------Test---------"<<endl;
     for (unsigned int i=0; i<this->n_dofs(); ++i)
         if (fabs(res(i)-delta_res(i)) > 1e-15)
         //if (fabs(res(i)) > 1e-20)
         cout<<i<<"  "<<delta_res(i)<<" vs "<<res(i)<<"      err "<<res(i)-delta_res(i)<<"   "<<sys_comp(i)<<endl;
     delta_res*=-1;
     delta_res.add(res);
     cout<<"Absolute error norm: "<<delta_res.l2_norm()<<endl;
     cout<<"Relative error norm: "<<delta_res.l2_norm()/delta_y.l2_norm()<<endl;
     cout<<"----------Done---------"<<endl;
     for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
         cout<<i<<" "<<comp_dom.support_points[i]<<"     sn "<<comp_dom.surface_nodes(i)<<"    ic "<<constraints.is_constrained(i)<<endl;
  //*/
std::cout<<"... Done preparing interpolated solution for restart"<<std::endl;


}




template <int dim>
Vector<double>& FreeSurface<dim>::get_diameters()
{
  if(diameters.size() != n_dofs())
    { 
      diameters.reinit(dofs_number);
      diameters.add(1000000);

      cell_it
	vector_cell = comp_dom.vector_dh.begin_active(),
	vector_endc = comp_dom.vector_dh.end();

      cell_it
	phi_cell = comp_dom.dh.begin_active(),
	phi_endc = comp_dom.dh.end();

      FEValues<dim-1,dim> fe_v(*comp_dom.mapping, comp_dom.fe, *comp_dom.quadrature,
			    update_JxW_values);
      const unsigned int n_q_points = fe_v.n_quadrature_points;
      std::vector<unsigned int> vector_local_dof_indices(comp_dom.vector_fe.dofs_per_cell);    
      std::vector<unsigned int> phi_local_dof_indices(comp_dom.fe.dofs_per_cell);  

      for (; phi_cell!=phi_endc,vector_cell!=vector_endc; ++phi_cell,++vector_cell)
	{
	  Assert(phi_cell->index() == vector_cell->index(), ExcInternalError());

	  phi_cell->get_dof_indices(phi_local_dof_indices);
	  vector_cell->get_dof_indices(vector_local_dof_indices);
	  for (unsigned int i=0; i<comp_dom.vector_fe.dofs_per_cell; ++i)
	    { 
	      diameters(vector_local_dof_indices[i]) =
		fmin(diameters(vector_local_dof_indices[i]),vector_cell->diameter());          
	    }
	  for (unsigned int i=0; i<comp_dom.fe.dofs_per_cell; ++i)
	    { 
	      diameters(phi_local_dof_indices[i]+comp_dom.vector_dh.n_dofs()) =
		fmin(diameters(phi_local_dof_indices[i]+comp_dom.vector_dh.n_dofs()),phi_cell->diameter());
              diameters(phi_local_dof_indices[i]+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()) =
		fmin(diameters(phi_local_dof_indices[i]+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()),phi_cell->diameter()); 
	    }

	}
    }
  return diameters;
}

template <int dim>
int FreeSurface<dim>::residual(Vector<double> &dst,  
	                       const Vector<double> &src_yy)
{
  double alpha = 0;
  // here nodes_alg_jac_x_delta won't be used, it's just to feed it with something
  // withoud having to generate a (big) dummy vector
  nodes_alg_jac_x_delta=0;
  residual_and_jacobian(current_time,dst,current_sol,src_yy,nodes_alg_jac_x_delta,0.0,false);


  dae_nonlin_residual = dst;
  dae_nonlin_residual*=-1;
  

  return 0;
}


template <int dim>
int FreeSurface<dim>::residual(const double t, 
			       Vector<double> &dst,  
			       const Vector<double> &src_yy,
			       const Vector<double> &src_yp)
{
  double alpha = 0;
  // here nodes_alg_jac_x_delta won't be used, it's just to feed it with something
  // withoud having to generate a (big) dummy vector
  residual_and_jacobian(t,dst,src_yy,src_yp,nodes_alg_jac_x_delta,alpha,false);


  dae_nonlin_residual = dst;
  dae_nonlin_residual*=-1;

  //for (unsigned int i=0; i<dst.size(); ++i)
  //    {
  //    cout<<i<<" -> "<<dst(i)<<endl;
  //    }
  

  return 0;
}



template <int dim>
void FreeSurface<dim>::vmult(Vector<double> &dst, const Vector<double> &src) const
{

  dst = 0;

  Vector<double> pphi(comp_dom.dh.n_dofs());
  Vector<double> dpphi_dn(comp_dom.dh.n_dofs());
  Vector<double> bem_pphi(comp_dom.dh.n_dofs());
  Vector<double> bem_dpphi_dn(comp_dom.dh.n_dofs());    

  for (unsigned int i=0;i<comp_dom.dh.n_dofs();++i)
      {
      pphi(i) = src(i+comp_dom.vector_dh.n_dofs());
      bem_pphi(i) = src(i+comp_dom.vector_dh.n_dofs());
      dpphi_dn(i) = src(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs());
      bem_dpphi_dn(i) = src(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs());
      }

//  const VectorView<double> phi(comp_dom.dh.n_dofs(),src.begin()+comp_dom.vector_dh.n_dofs());
//  const VectorView<double> dphi_dn(comp_dom.dh.n_dofs(),src.begin()+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs());

  //bem_phi = pphi; //(const Vector<double> &)phi;
  constraints.distribute(bem_pphi);

  //bem_dphi_dn = dpphi_dn;//(const Vector<double> &)dphi_dn;
  constraints.distribute(bem_dpphi_dn); 

  Vector<double> bem_bc(comp_dom.dh.n_dofs());
  for (unsigned int i=0; i<comp_dom.dh.n_dofs(); i++)
      {
      if (comp_dom.flags[i] & water)
         bem_bc(i) = bem_pphi(i);
      else
         bem_bc(i) = bem_dpphi_dn(i);
      }

  bem.solve_system(bem_pphi, bem_dpphi_dn, bem_bc);

  jacobian_dot_matrix.vmult(dst,src);
  dst*=alpha;
  jacobian_matrix.vmult_add(dst,src);

  for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
      {
      if ( //(!constraints.is_constrained(i)) &&
           (comp_dom.flags[i] & water) )
         {
         dst(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()) = dpphi_dn(i) - bem_dphi_dn(i);
         }
      else if ( //(!constraints.is_constrained(i)) &&
                (!(comp_dom.flags[i] & water)) )
         {
         dst(i+comp_dom.vector_dh.n_dofs()) = pphi(i) - bem_phi(i);
         }
      }





}


template <int dim>
int FreeSurface<dim>::jacobian(Vector<double> &dst,  
	                       const Vector<double> &src_yy,
		               const Vector<double> &src)
{

Vector<double> src_copy(src);
jacobian(current_time,dst,current_sol,src_yy,src_copy,1e7);

return 0;
}


template <int dim>
int FreeSurface<dim>::jacobian(const double t,
		       Vector<double> &dst,  
		       const Vector<double> &src_yy,
		       const Vector<double> &src_yp,
		       const Vector<double> &src,
		       const double alpha)
{
  
  //for (unsigned int i=0; i<dst.size(); ++i)
  //    cout<<"src("<<i<<") = "<<src(i)<<endl;
  //residual_and_jacobian(t,dst,src_yy,src_yp,src,alpha,true);

  dst = 0;

  const VectorView<double> phi(comp_dom.dh.n_dofs(),src.begin()+comp_dom.vector_dh.n_dofs());
  const VectorView<double> dphi_dn(comp_dom.dh.n_dofs(),src.begin()+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs());

  bem_phi = (const Vector<double> &)phi;
  constraints.distribute(bem_phi);

  bem_dphi_dn = (const Vector<double> &)dphi_dn;
  constraints.distribute(bem_dphi_dn); 

  Vector<double> bem_bc(comp_dom.dh.n_dofs());
  for (unsigned int i=0; i<comp_dom.dh.n_dofs(); i++)
      {
      if (comp_dom.flags[i] & water ||
          comp_dom.flags[i] & pressure)
         bem_bc(i) = phi(i);
      else
         bem_bc(i) = dphi_dn(i);
      }

  // trying a fix for transom stern nodes
  for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
      {
      if ( comp_dom.flags[i] & transom_on_water )
         {
         comp_dom.surface_nodes(i) = 0;
         comp_dom.other_nodes(i) = 1;
         std::set<unsigned int> duplicates = comp_dom.double_nodes_set[i];
         duplicates.erase(i); 
         bem_bc(i) = 0;
         for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
             {
             bem_bc(i) += bem_dphi_dn(*pos)/duplicates.size();
             }
         bem_dphi_dn(i) = bem_bc(i);
         }
      }

     // trying a fix for water/pressure nodes (water side)
      for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
          {
          if ( (comp_dom.flags[i] & water) &&  (comp_dom.flags[i] & near_pressure) )
             {
             std::set<unsigned int> duplicates = comp_dom.double_nodes_set[i];
             duplicates.erase(i);
             unsigned int count=0;
             bem_bc(i) = 0;
             for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
                 if (comp_dom.flags[*pos] & pressure)
                    count++;
                 
             for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
                 {
                 if (comp_dom.flags[*pos] & pressure)
                    {
                    bem_bc(i) += phi(*pos)/count;
                    }
                 }
             bem_phi(i) = bem_bc(i);
             }
          } 

  bem.solve_system(bem_phi, bem_dphi_dn, bem_bc);
  
/*
   cout<<"44!!!!!!"<<endl;
    for (SparsityPattern::iterator col=jacobian_sparsity_pattern.begin(1484); col!=jacobian_sparsity_pattern.end(1484); ++col)
        {
        unsigned int j = col->column(); 
        cout<<j<<"("<<jacobian_matrix(1484,j)<<") ";
        }
    cout<<endl;
//*/

//for (SparsityPattern::iterator col=jacobian_sparsity_pattern.begin(3196); col!=jacobian_sparsity_pattern.end(3196); ++col)
//        {
//        unsigned int j = col->column();
//        cout<<" "<<3196<<" "<<j<<" "<<jacobian_matrix(3196,j)<<" "<<alpha*jacobian_dot_matrix(3196,j)<<" * "<<src(j)<<endl;
//        }
  
  jacobian_dot_matrix.vmult(dst,src);

  dst*=alpha;
  jacobian_matrix.vmult_add(dst,src);

  for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
      {
      if ( (!constraints.is_constrained(i)) &&
           ((comp_dom.flags[i] & water) || (comp_dom.flags[i] & pressure)) &&
           (!(comp_dom.flags[i] & transom_on_water)))
         {//cout<<i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()<<" ***  "<<bem_dphi_dn(i)<<endl;
         dst(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()) = bem_dphi_dn(i) - dphi_dn(i);
         }
      else if ( (!constraints.is_constrained(i)) &&
                (!(comp_dom.flags[i] & transom_on_water)) )
         {
         dst(i+comp_dom.vector_dh.n_dofs()) = bem_phi(i) - phi(i);
         }
      else if ( (comp_dom.flags[i] & transom_on_water) )
         {
         dst(i+comp_dom.vector_dh.n_dofs()) = bem_phi(i) - phi(i);
         }
      }



  dae_linear_step_residual = dae_nonlin_residual;
  dae_linear_step_residual.add(dst);

  cout<<"Linear step residual: "<<dae_linear_step_residual.l2_norm()<<endl;

//  for (unsigned int i=0; i<dst.size(); ++i)
//      cout<<i<<" "<<dst(i)<<" "<<src_yy(i)<<" "<<src_yp(i)<<" "<<src(i)<<endl;


  static unsigned int jac_evals = 0;
     jac_evals++;
     std::cout << "Jacobian matrix-vector product evaluations: " << jac_evals << std::endl;


  return 0;
}   
   

				     /** This function computes either DAE residual
                                         or corrensponding Jacobian matrix vector product
                                         with vector src. Boolean is_jacobian determines if
                                         this function is used for residual (false) or for
                                         Jacobian (true). In residual case, src and alpha
                                         values assigned are disregarded.*/
template <int dim>
int FreeSurface<dim>::residual_and_jacobian(const double t,
		                            Vector<double> &dst,  
		                            const Vector<double> &src_yy,
		                            const Vector<double> &src_yp,
		                            const Vector<double> &src,
		                            const double alpha,
                                            const bool is_jacobian)
{
    jacobian_matrix = 0;
    jacobian_dot_matrix = 0;

    static double old_time = -1000;
    dst = 0;


    //if (t != old_time)
       //{
       //comp_dom.old_map_points = comp_dom.map_points;
       //std::string filename1 = ( "new_ts_mesh.vtu" );
       //output_results(filename1, t, src_yy, src_yp);
       //}
//*/




  const VectorView<double> nodes_positions(comp_dom.vector_dh.n_dofs(),src_yy.begin());
  const VectorView<double> nodes_velocities(comp_dom.vector_dh.n_dofs(),src_yp.begin());
  const VectorView<double> phi(comp_dom.dh.n_dofs(),src_yy.begin()+comp_dom.vector_dh.n_dofs());
  const VectorView<double> phi_time_derivs(comp_dom.dh.n_dofs(),src_yp.begin()+comp_dom.vector_dh.n_dofs());
  const VectorView<double> dphi_dn(comp_dom.dh.n_dofs(),src_yy.begin()+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs());
  const VectorView<double> dphi_dn_time_derivs(comp_dom.dh.n_dofs(),src_yp.begin()+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs());

  Point<3> hull_lin_vel;
  Point<3> hull_lin_displ;
  Point<3> hull_lin_vel_dot;
  Point<3> hull_lin_displ_dot;
  Point<3> hull_ang_vel;
  Point<3> hull_quat_vector;
  Point<3> hull_ang_vel_dot;
  Point<3> hull_quat_vector_dot;
  for (unsigned int d=0; d<3;++d)
      {
      hull_lin_vel(d) = src_yy(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+d);
      hull_lin_displ(d) = src_yy(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+3+d);
      hull_lin_vel_dot(d) = src_yp(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+d);
      hull_lin_displ_dot(d) = src_yp(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+3+d);
      hull_ang_vel(d) = src_yy(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+6+d);
      hull_quat_vector(d) = src_yy(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+9+d);
      hull_ang_vel_dot(d) = src_yp(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+6+d);
      hull_quat_vector_dot(d) = src_yp(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+9+d);
      }
  double hull_quat_scalar = src_yy(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+12);
  double hull_quat_scalar_dot = src_yp(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+12);

  //for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
  //       cout<<3*i<<" "<<nodes_velocities(3*i)<<" "<<nodes_velocities(3*i+1)<<" "<<nodes_velocities(3*i+2)<<endl;

    bool new_time_step = false;
    if (t != old_time) 
    {
    old_time = t;
    new_time_step = true;      
    }
 


  // TO BE REMOVED
  DphiDt_sys_rhs.reinit (comp_dom.dh.n_dofs()); 
  DphiDt_sys_rhs_2.reinit (comp_dom.dh.n_dofs());
  DphiDt_sys_rhs_3.reinit (comp_dom.dh.n_dofs());
  DphiDt_sys_rhs_4.reinit (comp_dom.dh.n_dofs());
  DphiDt_sys_matrix.reinit(DphiDt_sparsity_pattern);
  DphiDt_sys_matrix_2.reinit(DphiDt_sparsity_pattern);

  //if(old_time < t) 
   // {
    //std::string filename1 = ( "mesh.vtu" );
    //output_results(filename1, t, src_yy, src_yp);
  //  }
  
  
  AssertDimension(dst.size(), src_yy.size());
  wind.set_time(t);
  inflow_norm_potential_grad.set_time(t);
  Vector<double> instantWindValue(dim);
  Point<dim> zero(0,0,0);
  wind.vector_value(zero,instantWindValue);
  std::cout<<std::endl<<"Simulation time= "<<t<<"   Vinf= ";
  instantWindValue.print(std::cout,4,false,true);
  std::cout << "Ndofs ODE= " << dst.size();
  std::cout << "  Ndofs BEM= " << comp_dom.dh.n_dofs();
  std::cout<<std::endl;  
  wind.set_time(t);

  cout<<"Hull Rigid Displacement: "<<hull_lin_displ<<endl;
  cout<<"Restart Hull Rigid Displacement: "<<restart_hull_displacement<<endl;
  cout<<"Hull Rigid Velocity: "<<hull_lin_vel<<endl;
  cout<<"Hull Angular Velocity: "<<hull_ang_vel<<endl;

  // these two for loops set the jacobian for all the constrained (hanging node dofs)
  // first one is for position dofs  
  for (unsigned int i=0; i<comp_dom.vector_dh.n_dofs(); ++i)
      { 
      if (vector_constraints.is_constrained(i))
         {
         jacobian_matrix.add(i,i,1.0);
         std::vector< std::pair< unsigned int, double > > entries = *vector_constraints.get_constraint_entries(i);
         for (unsigned int k=0; k<entries.size(); ++k)
             {//cout<<"* "<<i<<endl;
             //cout<<i<<"  "<<entries[k].first<<" "<<-entries[k].second<<endl;
             jacobian_matrix.add(i,entries[k].first,-entries[k].second);
             }
         }
      }

    // second one is for potential and potential normal derivative dofs
    for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
      {//cout<<i<<" "<<comp_dom.surface_nodes(i)<<"  ("<<comp_dom.support_points[i]<<")"<<endl;
      if (constraints.is_constrained(i))
         {
         jacobian_matrix.add(i+comp_dom.vector_dh.n_dofs(),i+comp_dom.vector_dh.n_dofs(),-1.0);
         jacobian_matrix.add(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs(),i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs(),-1.0);
         std::vector< std::pair< unsigned int, double > > entries = *constraints.get_constraint_entries(i);
         for (unsigned int k=0; k<entries.size(); ++k)
             {//cout<<i<<"  "<<entries[k].first<<" "<<-entries[k].second<<endl;
             jacobian_matrix.add(i+comp_dom.vector_dh.n_dofs(),entries[k].first+comp_dom.vector_dh.n_dofs(),entries[k].second);
             jacobian_matrix.add(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs(),
                                 entries[k].first+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs(),entries[k].second);
             }
         }
      }

//////////////////////////////////////////////////////////////////////
  // here we take care of the geometric part of the variables
//////////////////////////////////////////////////////////////////////


  // these lines take care of the hull (and consequential hull internal grid nodes) displacement
  // these lines will only act when the time step is new
  
  //if (new_time_step)
  //   {
     //double rigid_vertical_displacement = 0.01*sin(2.0*dealii::numbers::PI*t);
     //double rigid_angular_displacement = 0.01*sin(1.0*dealii::numbers::PI*t);
     //cout<<"Displacement: "<<rigid_vertical_displacement<<"   "<<sin(2.0*(dealii::numbers::PI)*t)<<endl;
     // this is the displacement of the hull geometry (all projectors and smoothers will adapt
     // to the change)

     
     gp_Trsf reference_to_current_transformation = comp_dom.boat_model.set_current_position(hull_lin_displ,
                                                                                            hull_quat_scalar,
                                                                                            hull_quat_vector);

     // this moves the rigid_motion_map_points vector (target positions) for the internal nodes of the hull 
     // mesh (all nodes except for the free surface ones)
     for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
         {
         if ( ((comp_dom.flags[i] & boat) &&
              !(comp_dom.flags[i] & near_water) ) ||
              (comp_dom.flags[i] & transom_on_water) )
            {
            //if (fabs(t-0.2) <1e-5)
            //cout<<"BEFORE: "<<comp_dom.ref_points[3*i](0)+comp_dom.old_map_points(3*i)<<" "
            //                <<comp_dom.ref_points[3*i](1)+comp_dom.old_map_points(3*i+1)<<" "
            //                <<comp_dom.ref_points[3*i](2)+comp_dom.old_map_points(3*i+2)<<endl;
            gp_Pnt original_boat_mesh_point = Pnt(Point<3>(comp_dom.ref_points[3*i](0)+comp_dom.old_map_points(3*i),
                                                           comp_dom.ref_points[3*i](1)+comp_dom.old_map_points(3*i+1),
                                                           comp_dom.ref_points[3*i](2)+comp_dom.old_map_points(3*i+2)));
            gp_Pnt boat_mesh_point = original_boat_mesh_point;
            // we first take this point (which is in the RESTART hull location) and transform it to be in the
            // REFERENCE configuration 
            boat_mesh_point.Transform(restart_hull_location.Inverted());

            // now we use sacado to compute the residual at this dof, along with its derivatives with respect to 
            // the 7 dofs associated to the rigid linear and angular displacements
            fad_double s_x,s_y,s_z,v_x,v_y,v_z,s;

            s_x = hull_lin_displ(0);
            s_y = hull_lin_displ(1);
            s_z = hull_lin_displ(2);
            v_x = hull_quat_vector(0);
            v_y = hull_quat_vector(1);
            v_z = hull_quat_vector(2);
            s = hull_quat_scalar;

            s_x.diff(0,7);
            s_y.diff(1,7);
            s_z.diff(2,7);
            v_x.diff(3,7);
            v_y.diff(4,7);
            v_z.diff(5,7);
            s.diff(6,7);

            Point<3,fad_double> ref_baricenter_pos(comp_dom.boat_model.reference_hull_baricenter(0),
                                                   comp_dom.boat_model.reference_hull_baricenter(1),
                                                   comp_dom.boat_model.reference_hull_baricenter(2));
   
            Point<3,fad_double> baricenter_pos(s_x+comp_dom.boat_model.reference_hull_baricenter(0),
                                               s_y+comp_dom.boat_model.reference_hull_baricenter(1),
                                               s_z+comp_dom.boat_model.reference_hull_baricenter(2));

            Point<3,fad_double> RotMatRow1(1-2*v_y*v_y-2*v_z*v_z, 2*v_x*v_y-2*s*v_z, 2*v_x*v_z+2*s*v_y);
            Point<3,fad_double> RotMatRow2(2*v_x*v_y+2*s*v_z, 1-2*v_x*v_x-2*v_z*v_z, 2*v_y*v_z-2*s*v_x);
            Point<3,fad_double> RotMatRow3(2*v_x*v_z-2*s*v_y, 2*v_y*v_z+2*s*v_x, 1-2*v_y*v_y-2*v_x*v_x);

            Point<3,fad_double> ref_point_pos(boat_mesh_point.X(),boat_mesh_point.Y(),boat_mesh_point.Z());
            Point<3,fad_double> rigid_lin_displ(s_x,s_y,s_z);
            Point<3,fad_double> target_point_pos(RotMatRow1*(ref_point_pos+(fad_double(-1.0))*ref_baricenter_pos),
                                                 RotMatRow2*(ref_point_pos+(fad_double(-1.0))*ref_baricenter_pos),
                                                 RotMatRow3*(ref_point_pos+(fad_double(-1.0))*ref_baricenter_pos));
            target_point_pos += baricenter_pos;
            // now we have the point on the REFERENCE hull, and transform it to go onto the CURRENT hull
            //boat_mesh_point.Transform(reference_to_current_transformation);
            // the rigid motion map points is the difference between the reference point position and the
            // current (rigidly displaced) node position
            
            comp_dom.rigid_motion_map_points(3*i) = target_point_pos(0).val()-ref_point_pos(0).val();
            comp_dom.rigid_motion_map_points(3*i+1) = target_point_pos(1).val()-ref_point_pos(1).val();
            comp_dom.rigid_motion_map_points(3*i+2) = target_point_pos(2).val()-ref_point_pos(2).val();

            //if (fabs(t-0.2) <1e-5)
            //cout<<"AFTER: "<<Pnt(boat_mesh_point)<<" vs "<<Pnt(original_boat_mesh_point)<<endl;
            //cout<<"RMMP: "<<comp_dom.rigid_motion_map_points(3*i)<<" "
            //              <<comp_dom.rigid_motion_map_points(3*i+1)<<" "
            //              <<comp_dom.rigid_motion_map_points(3*i+2)<<endl;
            //cout<<"NN: "<<comp_dom.rigid_motion_map_points(3*i)+nodes_positions(3*i)+comp_dom.ref_points[3*i](0)<<" "
            //            <<comp_dom.rigid_motion_map_points(3*i+1)+nodes_positions(3*i+1)+comp_dom.ref_points[3*i](1)<<" "
            //            <<comp_dom.rigid_motion_map_points(3*i+2)+nodes_positions(3*i+2)+comp_dom.ref_points[3*i](2)<<endl;
            }
         }
  
  //   }

  



  Vector<double> full_map_points(comp_dom.rigid_motion_map_points);
  full_map_points.add(nodes_positions);
  comp_dom.update_mapping(full_map_points);
  comp_dom.update_support_points();

  
  //we work on a local COPY of map_points
  working_map_points = comp_dom.map_points;
  nodes_pos_res = 0;

  //we enforce constraint on the new geometry assigned by the DAE solver
  vector_constraints.distribute(working_map_points);

  // as for now x and y coordinates of nodes are not moved (no surface smoothing)
  // so on x and y coordinates (except for water nodes) we only need to put a one
  // on the matrix diagonal
  for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
      {
      if ( !(comp_dom.flags[i] & near_boat) &&
           (comp_dom.flags[i] & water) &&
           (comp_dom.flags[i] & edge) &&
           !(comp_dom.flags[i] & transom_on_water) &&
           (!constraints.is_constrained(i)))
         {
         working_map_points(3*i) = comp_dom.old_map_points(3*i);
         working_map_points(3*i+1) = comp_dom.old_map_points(3*i+1);
         jacobian_matrix.add(3*i,3*i,1.0);
         jacobian_matrix.add(3*i+1,3*i+1,1.0);
         if (comp_dom.flags[i] & near_inflow)
            {
            working_map_points(3*i+2) = comp_dom.old_map_points(3*i+2);
            jacobian_matrix.add(3*i+2,3*i+2,1.0);
            }
         //cout<<"& "<<3*i<<" "<<3*i+1<<endl;
         }
      else if ( !(comp_dom.flags[i] & near_water) &&
                (comp_dom.flags[i] & boat) &&
               (!constraints.is_constrained(i)))
         {
         gp_Pnt original_boat_mesh_point = Pnt(Point<3>(comp_dom.ref_points[3*i](0)+comp_dom.old_map_points(3*i),
                                                        comp_dom.ref_points[3*i](1)+comp_dom.old_map_points(3*i+1),
                                                        comp_dom.ref_points[3*i](2)+comp_dom.old_map_points(3*i+2)));
         gp_Pnt boat_mesh_point = original_boat_mesh_point;
         // we first take this point (which is in the RESTART hull location) and transform it to be in the
         // REFERENCE configuration 
         boat_mesh_point.Transform(restart_hull_location.Inverted());
         Point<3,fad_double> ref_point_pos(boat_mesh_point.X(),boat_mesh_point.Y(),boat_mesh_point.Z());
         
         working_map_points(3*i) = ref_point_pos(0).val()-comp_dom.ref_points[3*i](0)+comp_dom.rigid_motion_map_points(3*i);
         working_map_points(3*i+1) = ref_point_pos(1).val()-comp_dom.ref_points[3*i](1)+comp_dom.rigid_motion_map_points(3*i+1);
         working_map_points(3*i+2) = ref_point_pos(2).val()-comp_dom.ref_points[3*i](2)+comp_dom.rigid_motion_map_points(3*i+2);
         jacobian_matrix.add(3*i,3*i,1.0);
         jacobian_matrix.add(3*i+1,3*i+1,1.0);
         jacobian_matrix.add(3*i+2,3*i+2,1.0);
         }
      else if ( !(comp_dom.flags[i] & near_water) &&
                !(comp_dom.flags[i] & water) &&
                !(comp_dom.flags[i] & boat) &&
               (!constraints.is_constrained(i)))
         {
         working_map_points(3*i) = comp_dom.old_map_points(3*i);
         working_map_points(3*i+1) = comp_dom.old_map_points(3*i+1);
         working_map_points(3*i+2) = comp_dom.old_map_points(3*i+2);
         jacobian_matrix.add(3*i,3*i,1.0);
         jacobian_matrix.add(3*i+1,3*i+1,1.0);
         jacobian_matrix.add(3*i+2,3*i+2,1.0);
         //cout<<"&& "<<3*i<<" "<<3*i+1<<" "<<3*i+2<<endl;         
         }
      else if ((comp_dom.flags[i] & transom_on_water) )
         {
         gp_Pnt original_boat_mesh_point = Pnt(Point<3>(comp_dom.ref_points[3*i](0)+comp_dom.old_map_points(3*i),
                                                        comp_dom.ref_points[3*i](1)+comp_dom.old_map_points(3*i+1),
                                                        comp_dom.ref_points[3*i](2)+comp_dom.old_map_points(3*i+2)));
         gp_Pnt boat_mesh_point = original_boat_mesh_point;
         // we first take this point (which is in the RESTART hull location) and transform it to be in the
         // REFERENCE configuration 
         boat_mesh_point.Transform(restart_hull_location.Inverted());
         Point<3,fad_double> ref_point_pos(boat_mesh_point.X(),boat_mesh_point.Y(),boat_mesh_point.Z());
         
         working_map_points(3*i) = ref_point_pos(0).val()-comp_dom.ref_points[3*i](0)+comp_dom.rigid_motion_map_points(3*i);
         working_map_points(3*i+1) = ref_point_pos(1).val()-comp_dom.ref_points[3*i](1)+comp_dom.rigid_motion_map_points(3*i+1);
         working_map_points(3*i+2) = ref_point_pos(2).val()-comp_dom.ref_points[3*i](2)+comp_dom.rigid_motion_map_points(3*i+2);
         jacobian_matrix.add(3*i,3*i,1.0);
         jacobian_matrix.add(3*i+1,3*i+1,1.0);
         jacobian_matrix.add(3*i+2,3*i+2,1.0);
         }          
      }

  // blending factor is needed to avoid that right after restart
  // the free surface mesh smoothing causes infinite horizontal nodes velocity:
  // in fact, at restars water line nodes are moved along water lines for surface
  // smoothing, while surrounding free surface nodes are still. as soon as
  // time integration is restarted the smoothing forces them to move of a finite
  // displacement over a zero time span. using this blending factor could be avoided if
  // surface smoothing was done at restarts too. but this is somehow troublesome as
  // the interpolation if potential on the new (smoothed) surface mesh is troublesome 
  double blend_factor = 0.0;
  if (restart_flag)
     {
     }
  else
     {

     if (t - last_remesh_time < 0.5*remeshing_period)
        blend_factor = sin(3.141592654*(t-last_remesh_time)/remeshing_period);
     else
        blend_factor = 1.0;
     //std::cout<<"t "<<t<<"  last_remesh_time "<<last_remesh_time<<"  remeshing_period/5 "<<remeshing_period/5<<std::endl;
     std::cout<<"blend_factor = "<<blend_factor<<std::endl;
     }

    //this takes care of the right water line nodes projection (without smoothing)
     if (!comp_dom.no_boat) 
     for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
         { 
         fad_double s_x,s_y,s_z,v_x,v_y,v_z,s;

         s_x = hull_lin_displ(0);
         s_y = hull_lin_displ(1);
         s_z = hull_lin_displ(2);
         v_x = hull_quat_vector(0);
         v_y = hull_quat_vector(1);
         v_z = hull_quat_vector(2);
         s = hull_quat_scalar;

         s_x.diff(0,7);
         s_y.diff(1,7);
         s_z.diff(2,7);
         v_x.diff(3,7);
         v_y.diff(4,7);
         v_z.diff(5,7);
         s.diff(6,7);

         Point<3,fad_double> ref_baricenter_pos(comp_dom.boat_model.reference_hull_baricenter(0),
                                                comp_dom.boat_model.reference_hull_baricenter(1),
                                                comp_dom.boat_model.reference_hull_baricenter(2));


         Point<3,fad_double> baricenter_pos(s_x+comp_dom.boat_model.reference_hull_baricenter(0),
                                            s_y+comp_dom.boat_model.reference_hull_baricenter(1),
                                            s_z+comp_dom.boat_model.reference_hull_baricenter(2));

         Point<3,fad_double> RotMatRow1(1-2*v_y*v_y-2*v_z*v_z, 2*v_x*v_y-2*s*v_z, 2*v_x*v_z+2*s*v_y);
         Point<3,fad_double> RotMatRow2(2*v_x*v_y+2*s*v_z, 1-2*v_x*v_x-2*v_z*v_z, 2*v_y*v_z-2*s*v_x);
         Point<3,fad_double> RotMatRow3(2*v_x*v_z-2*s*v_y, 2*v_y*v_z+2*s*v_x, 1-2*v_y*v_y-2*v_x*v_x);

         if ( (comp_dom.flags[i] & water) &&
              (comp_dom.flags[i] & near_boat) &&
              (comp_dom.flags[i] & right_side) &&
              !(comp_dom.flags[i] & transom_on_water) &&
              (comp_dom.moving_point_ids[3] != i) &&
              (comp_dom.moving_point_ids[4] != i) &&
              (comp_dom.moving_point_ids[5] != i) &&
              (comp_dom.moving_point_ids[6] != i) ) // to avoid the bow and stern node
            {//cout<<"**** "<<i<<endl;
            Point<3> proj_node;
            Point<3> iges_normal;
            double iges_curvature;
            Point<3> direction(comp_dom.iges_normals[i](0),comp_dom.iges_normals[i](1),0.0);
            //if (fabs(comp_dom.old_iges_normals[i](0)) > 0.001)
            //   cout<<3*i+1<<"   "<<comp_dom.support_points[i]<<"   ("<<comp_dom.old_iges_normals[i]<<")"<<endl;
            if (fabs(comp_dom.old_iges_normals[i](0))<sqrt(3)/3*fabs(comp_dom.old_iges_normals[i](1)))
               direction(0) = 0.0;
            else
               direction(1) = 0.0;
            //if (fabs(comp_dom.iges_normals[i](0))<0.001)
            //   cout<<3*i<<"  dir:    ("<<direction<<")"<<endl;
            comp_dom.boat_model.boat_water_line_right->assigned_axis_projection_and_diff_forms(proj_node,
                                                                                               comp_dom.iges_normals[i],
                                                                                               comp_dom.iges_mean_curvatures[i],
                                                                                               comp_dom.support_points[i],
                                                                                               direction);  // hor normal dir projection
            //comp_dom.boat_model.boat_water_line_right->axis_projection_and_diff_forms(proj_node,
            //                                                                          comp_dom.iges_normals[i],
            //                                                                          comp_dom.iges_mean_curvatures[i],
            //                                                                          comp_dom.support_points[i]);  // y axis projection
            if (fabs(comp_dom.old_iges_normals[i](0))<sqrt(3)/3*fabs(comp_dom.old_iges_normals[i](1)))
               {
               working_map_points(3*i) = comp_dom.old_map_points(3*i); // x of the node must not change
               working_map_points(3*i+1) = proj_node(1) - comp_dom.ref_points[3*i](1);
               }
            else
               {
               working_map_points(3*i) = proj_node(0) - comp_dom.ref_points[3*i](0);
               working_map_points(3*i+1) = comp_dom.old_map_points(3*i+1); // y of the node must not change
               }            
            working_map_points(3*i+2) = proj_node(2) - comp_dom.ref_points[3*i](2);

            // here we compute the derivatives: we are basically computing the same residual on the tangent plane, which
            // will have the same derivatives as the surface
            gp_Pnt ref_point(proj_node(0),proj_node(1),proj_node(2));
            ref_point.Transform(reference_to_current_transformation.Inverted());
            gp_Vec curr_normal(comp_dom.iges_normals[i](0),comp_dom.iges_normals[i](1),comp_dom.iges_normals[i](2));
            gp_Dir ref_normal_dir(curr_normal);
            ref_normal_dir.Transform(reference_to_current_transformation.Inverted());
            
            Point<3,fad_double> ref_point_pos(ref_point.X(),ref_point.Y(),ref_point.Z());
            Point<3,fad_double> rigid_lin_displ(s_x,s_y,s_z);
            Point<3,fad_double> target_point_pos(RotMatRow1*(ref_point_pos+(fad_double(-1.0))*ref_baricenter_pos),
                                                 RotMatRow2*(ref_point_pos+(fad_double(-1.0))*ref_baricenter_pos),
                                                 RotMatRow3*(ref_point_pos+(fad_double(-1.0))*ref_baricenter_pos));
            target_point_pos += baricenter_pos;
            Point<3,fad_double> ref_normal(ref_normal_dir.X(),ref_normal_dir.Y(),ref_normal_dir.Z());
            Point<3,fad_double> target_point_normal(RotMatRow1*ref_normal,
                                                    RotMatRow2*ref_normal,
                                                    RotMatRow3*ref_normal);
            Point<3,fad_double> guessed_point(comp_dom.support_points[i](0),comp_dom.support_points[i](1),comp_dom.support_points[i](2));


            if (fabs(comp_dom.old_iges_normals[i](0))<sqrt(3)/3*fabs(comp_dom.old_iges_normals[i](1)))
               {
               fad_double residual = guessed_point(1)+
                                     target_point_normal(0)/target_point_normal(1)*(guessed_point(0)-target_point_pos(0))+
                                     target_point_normal(2)/target_point_normal(1)*(guessed_point(2)-target_point_pos(2))-
                                     target_point_pos(1);
               jacobian_matrix.add(3*i,3*i,1.0);
               jacobian_matrix.add(3*i+1,3*i+1,1.0);
               jacobian_matrix.add(3*i+1,3*i+2,comp_dom.iges_normals[i](2)/comp_dom.iges_normals[i](1));
               jacobian_matrix.add(3*i+1,3*i,comp_dom.iges_normals[i](0)/comp_dom.iges_normals[i](1));
               for (unsigned int d=0; d<3;++d)
                   jacobian_matrix.add(3*i+1,comp_dom.vector_dh.n_dofs()+
                                             comp_dom.dh.n_dofs()+
                                             comp_dom.dh.n_dofs()+
                                             3+d,residual.fastAccessDx(d));
               for (unsigned int d=0; d<4;++d)
                   jacobian_matrix.add(3*i+1,comp_dom.vector_dh.n_dofs()+
                                             comp_dom.dh.n_dofs()+
                                             comp_dom.dh.n_dofs()+
                                             9+d,residual.fastAccessDx(3+d));
               }
            else
               {
               fad_double residual = guessed_point(0)+
                                     target_point_normal(1)/target_point_normal(0)*(guessed_point(1)-target_point_pos(1))+
                                     target_point_normal(2)/target_point_normal(0)*(guessed_point(2)-target_point_pos(2))-
                                     target_point_pos(1);
               jacobian_matrix.add(3*i,3*i,1.0); 
               jacobian_matrix.add(3*i,3*i+2,comp_dom.iges_normals[i](2)/comp_dom.iges_normals[i](0));
               jacobian_matrix.add(3*i,3*i+1,comp_dom.iges_normals[i](1)/comp_dom.iges_normals[i](0));
               jacobian_matrix.add(3*i+1,3*i+1,1.0);
               for (unsigned int d=0; d<3;++d)
                   jacobian_matrix.add(3*i,comp_dom.vector_dh.n_dofs()+
                                           comp_dom.dh.n_dofs()+
                                           comp_dom.dh.n_dofs()+
                                           3+d,residual.fastAccessDx(d));
               for (unsigned int d=0; d<4;++d)
                   jacobian_matrix.add(3*i,comp_dom.vector_dh.n_dofs()+
                                           comp_dom.dh.n_dofs()+
                                           comp_dom.dh.n_dofs()+
                                           9+d,residual.fastAccessDx(3+d));
               }

                        // we're doing this thing on the water side, but the iges_normal and iges_mean curvature belong to the boat side
            std::set<unsigned int> duplicates = comp_dom.double_nodes_set[i];
            for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
                {
	        comp_dom.iges_normals[*pos] = comp_dom.iges_normals[i];
                comp_dom.iges_mean_curvatures[*pos] = comp_dom.iges_mean_curvatures[i];
                }
            }              
         }

     //this takes care of the left water line nodes projection (without smoothing)
     if (!comp_dom.no_boat) 
     for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
         { 
         fad_double s_x,s_y,s_z,v_x,v_y,v_z,s;

         s_x = hull_lin_displ(0);
         s_y = hull_lin_displ(1);
         s_z = hull_lin_displ(2);
         v_x = hull_quat_vector(0);
         v_y = hull_quat_vector(1);
         v_z = hull_quat_vector(2);
         s = hull_quat_scalar;

         s_x.diff(0,7);
         s_y.diff(1,7);
         s_z.diff(2,7);
         v_x.diff(3,7);
         v_y.diff(4,7);
         v_z.diff(5,7);
         s.diff(6,7);

         Point<3,fad_double> ref_baricenter_pos(comp_dom.boat_model.reference_hull_baricenter(0),
                                                comp_dom.boat_model.reference_hull_baricenter(1),
                                                comp_dom.boat_model.reference_hull_baricenter(2));

         Point<3,fad_double> baricenter_pos(s_x+comp_dom.boat_model.reference_hull_baricenter(0),
                                            s_y+comp_dom.boat_model.reference_hull_baricenter(1),
                                            s_z+comp_dom.boat_model.reference_hull_baricenter(2));

         Point<3,fad_double> RotMatRow1(1-2*v_y*v_y-2*v_z*v_z, 2*v_x*v_y-2*s*v_z, 2*v_x*v_z+2*s*v_y);
         Point<3,fad_double> RotMatRow2(2*v_x*v_y+2*s*v_z, 1-2*v_x*v_x-2*v_z*v_z, 2*v_y*v_z-2*s*v_x);
         Point<3,fad_double> RotMatRow3(2*v_x*v_z-2*s*v_y, 2*v_y*v_z+2*s*v_x, 1-2*v_y*v_y-2*v_x*v_x);

         if ( (comp_dom.flags[i] & water) &&
              (comp_dom.flags[i] & near_boat) &&
              (comp_dom.flags[i] & left_side) &&
              !(comp_dom.flags[i] & transom_on_water) &&
              (comp_dom.moving_point_ids[3] != i) &&
              (comp_dom.moving_point_ids[4] != i) &&
              (comp_dom.moving_point_ids[5] != i) &&
              (comp_dom.moving_point_ids[6] != i) ) // to avoid the bow and stern node
            {//cout<<"**** "<<i<<endl;
            Point<3> proj_node;
            Point<3> iges_normal;
            double iges_curvature;
            Point<3> direction(comp_dom.iges_normals[i](0),comp_dom.iges_normals[i](1),0.0);
            if (fabs(comp_dom.old_iges_normals[i](0))<sqrt(3)/3*fabs(comp_dom.old_iges_normals[i](1)))
               direction(0) = 0.0;
            else
               direction(1) = 0.0;
            comp_dom.boat_model.boat_water_line_left->assigned_axis_projection_and_diff_forms(proj_node,
                                                                                              comp_dom.iges_normals[i],
                                                                                              comp_dom.iges_mean_curvatures[i],
                                                                                              comp_dom.support_points[i],
                                                                                              direction);  // hor normal dir projection
            //comp_dom.boat_model.boat_water_line_left->axis_projection_and_diff_forms(proj_node,
            //                                                                         comp_dom.iges_normals[i],
            //                                                                         comp_dom.iges_mean_curvatures[i],
            //                                                                         comp_dom.support_points[i]);  // y axis projection
            if (fabs(comp_dom.old_iges_normals[i](0))<sqrt(3)/3*fabs(comp_dom.old_iges_normals[i](1)))
               {
               working_map_points(3*i) = comp_dom.old_map_points(3*i); // x of the node must not change
               working_map_points(3*i+1) = proj_node(1) - comp_dom.ref_points[3*i](1);
               }
            else
               {
               working_map_points(3*i) = proj_node(0) - comp_dom.ref_points[3*i](0);
               working_map_points(3*i+1) = comp_dom.old_map_points(3*i+1); // y of the node must not change
               }            
            working_map_points(3*i+2) = proj_node(2) - comp_dom.ref_points[3*i](2);

            // here we compute the derivatives: we are basically computing the same residual on the tangent plane, which
            // will have the same derivatives as the surface
            gp_Pnt ref_point(proj_node(0),proj_node(1),proj_node(2));
            ref_point.Transform(reference_to_current_transformation.Inverted());
            gp_Vec curr_normal(comp_dom.iges_normals[i](0),comp_dom.iges_normals[i](1),comp_dom.iges_normals[i](2));
            gp_Dir ref_normal_dir(curr_normal);
            ref_normal_dir.Transform(reference_to_current_transformation.Inverted());
            
            Point<3,fad_double> ref_point_pos(ref_point.X(),ref_point.Y(),ref_point.Z());
            Point<3,fad_double> rigid_lin_displ(s_x,s_y,s_z);
            Point<3,fad_double> target_point_pos(RotMatRow1*(ref_point_pos+(fad_double(-1.0))*ref_baricenter_pos),
                                                 RotMatRow2*(ref_point_pos+(fad_double(-1.0))*ref_baricenter_pos),
                                                 RotMatRow3*(ref_point_pos+(fad_double(-1.0))*ref_baricenter_pos));
            target_point_pos += baricenter_pos;
            Point<3,fad_double> ref_normal(ref_normal_dir.X(),ref_normal_dir.Y(),ref_normal_dir.Z());
            Point<3,fad_double> target_point_normal(RotMatRow1*ref_normal,
                                                    RotMatRow2*ref_normal,
                                                    RotMatRow3*ref_normal);
            Point<3,fad_double> guessed_point(comp_dom.support_points[i](0),comp_dom.support_points[i](1),comp_dom.support_points[i](2));

            if (fabs(comp_dom.old_iges_normals[i](0))<sqrt(3)/3*fabs(comp_dom.old_iges_normals[i](1)))
               {
               fad_double residual = guessed_point(1)+
                                     target_point_normal(0)/target_point_normal(1)*(guessed_point(0)-target_point_pos(0))+
                                     target_point_normal(2)/target_point_normal(1)*(guessed_point(2)-target_point_pos(2))-
                                     target_point_pos(1);
               
               jacobian_matrix.add(3*i,3*i,1.0);
               jacobian_matrix.add(3*i+1,3*i+1,1.0);
               jacobian_matrix.add(3*i+1,3*i+2,comp_dom.iges_normals[i](2)/comp_dom.iges_normals[i](1));
               jacobian_matrix.add(3*i+1,3*i,comp_dom.iges_normals[i](0)/comp_dom.iges_normals[i](1));
               for (unsigned int d=0; d<3;++d)
                   jacobian_matrix.add(3*i+1,comp_dom.vector_dh.n_dofs()+
                                             comp_dom.dh.n_dofs()+
                                             comp_dom.dh.n_dofs()+
                                             3+d,residual.fastAccessDx(d));
               for (unsigned int d=0; d<4;++d)
                   jacobian_matrix.add(3*i+1,comp_dom.vector_dh.n_dofs()+
                                             comp_dom.dh.n_dofs()+
                                             comp_dom.dh.n_dofs()+
                                             9+d,residual.fastAccessDx(3+d));
               }
            else
               {
               fad_double residual = guessed_point(0)+
                                     target_point_normal(1)/target_point_normal(0)*(guessed_point(1)-target_point_pos(1))+
                                     target_point_normal(2)/target_point_normal(0)*(guessed_point(2)-target_point_pos(2))-
                                     target_point_pos(1);
               jacobian_matrix.add(3*i,3*i,1.0); 
               jacobian_matrix.add(3*i,3*i+2,comp_dom.iges_normals[i](2)/comp_dom.iges_normals[i](0));
               jacobian_matrix.add(3*i,3*i+1,comp_dom.iges_normals[i](1)/comp_dom.iges_normals[i](0));
               jacobian_matrix.add(3*i+1,3*i+1,1.0);
               for (unsigned int d=0; d<3;++d)
                   jacobian_matrix.add(3*i,comp_dom.vector_dh.n_dofs()+
                                           comp_dom.dh.n_dofs()+
                                           comp_dom.dh.n_dofs()+
                                           3+d,residual.fastAccessDx(d));
               for (unsigned int d=0; d<4;++d)
                   jacobian_matrix.add(3*i,comp_dom.vector_dh.n_dofs()+
                                           comp_dom.dh.n_dofs()+
                                           comp_dom.dh.n_dofs()+
                                           9+d,residual.fastAccessDx(3+d));
               }
            //cout<<i<<"   "<<temp_src(3*i+1)<<"   ("<<comp_dom.iges_normals[i]<<")"<<endl;

            // we're doing this thing on the water side, but the iges_normal and iges_mean curvature belong to the boat side
            std::set<unsigned int> duplicates = comp_dom.double_nodes_set[i];
            for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
                {
	        comp_dom.iges_normals[*pos] = comp_dom.iges_normals[i];
                comp_dom.iges_mean_curvatures[*pos] = comp_dom.iges_mean_curvatures[i];
                }
            }              
         }

     //this takes care of the bow and stern nodes
     if (!comp_dom.no_boat) 
        for (unsigned int k=3; k<7; ++k)
            { 
            fad_double s_x,s_y,s_z,v_x,v_y,v_z,s;

            s_x = hull_lin_displ(0);
            s_y = hull_lin_displ(1);
            s_z = hull_lin_displ(2);
            v_x = hull_quat_vector(0);
            v_y = hull_quat_vector(1);
            v_z = hull_quat_vector(2);
            s = hull_quat_scalar;

            s_x.diff(0,7);
            s_y.diff(1,7);
            s_z.diff(2,7);
            v_x.diff(3,7);
            v_y.diff(4,7);
            v_z.diff(5,7);
            s.diff(6,7);

            Point<3,fad_double> ref_baricenter_pos(comp_dom.boat_model.reference_hull_baricenter(0),
                                                   comp_dom.boat_model.reference_hull_baricenter(1),
                                                   comp_dom.boat_model.reference_hull_baricenter(2));

            Point<3,fad_double> baricenter_pos(s_x+comp_dom.boat_model.reference_hull_baricenter(0),
                                               s_y+comp_dom.boat_model.reference_hull_baricenter(1),
                                               s_z+comp_dom.boat_model.reference_hull_baricenter(2));

            Point<3,fad_double> RotMatRow1(1-2*v_y*v_y-2*v_z*v_z, 2*v_x*v_y-2*s*v_z, 2*v_x*v_z+2*s*v_y);
            Point<3,fad_double> RotMatRow2(2*v_x*v_y+2*s*v_z, 1-2*v_x*v_x-2*v_z*v_z, 2*v_y*v_z-2*s*v_x);
            Point<3,fad_double> RotMatRow3(2*v_x*v_z-2*s*v_y, 2*v_y*v_z+2*s*v_x, 1-2*v_y*v_y-2*v_x*v_x);

            unsigned int i = comp_dom.moving_point_ids[k];
            //cout<<"Moving point id: "<<i<<endl;
            {
            Point <3> dP0 = comp_dom.support_points[i];
            Point <3> dP;
         				   //this is the horizontal plane
            Handle(Geom_Plane) horPlane = new Geom_Plane(0.,0.,1.,-dP0(2));
            Handle(Geom_Curve) curve;
            TopLoc_Location L = comp_dom.boat_model.current_loc;
            TopLoc_Location L_inv = L.Inverted();
            //gp_Pnt test_point(0.0,0.0,0.0);
            //test_point.Transform(L_inv.Transformation());
            //cout<<"####### "<<Pnt(test_point)<<endl;
            horPlane->Transform(L_inv.Transformation());
            if (comp_dom.boat_model.is_transom)
               {
               if (k==3 || k==4)
                  curve = comp_dom.boat_model.equiv_keel_bspline;
               else if (k == 6)
                  curve = comp_dom.boat_model.left_transom_bspline;
               else
                  curve = comp_dom.boat_model.right_transom_bspline;
               }
            else
               {
               curve = comp_dom.boat_model.equiv_keel_bspline;
               }

            TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(curve);
            edge.Location(L);
            BRepAdaptor_Curve AC(edge);
            gp_Pnt P;
            gp_Vec V1;
            GeomAPI_IntCS Intersector(curve, horPlane);
            int npoints = Intersector.NbPoints();

            AssertThrow((npoints != 0), ExcMessage("Keel or transom curve is not intersecting with horizontal plane!"));
            double minDistance=1e7;
            double t,u,v;
            for (int j=0; j<npoints;++j)
                {
                gp_Pnt int_point = Intersector.Point(j+1);
                int_point.Transform(L.Transformation());
                Point<3> inters = Pnt(int_point);
                Intersector.Parameters(j+1,u,v,t);
                if (dP0.distance(inters) < minDistance)
                   {
                   minDistance = dP0.distance(inters);
                   dP = inters;
                   AC.D1(t,P,V1);
                   }
                }



            //cout<<"Check plane-curve intersection:"<<endl;
            //cout<<"Origin: "<<dP0<<"   Proj: "<<dP<<"  dist: "<<minDistance<<endl;
            //cout<<Pnt(P)<<endl;
            /*
            // here temporarily for kcs hull tests
            if (minDistance > 0.5*comp_dom.boat_model.boatWetLength)
               {
               Standard_Real First = curve->FirstParameter();
               Standard_Real Last = curve->LastParameter();
               gp_Pnt PIn(0.0,0.0,0.0);
               gp_Pnt PFin(0.0,0.0,0.0);
               gp_Vec VIn;
               gp_Vec VFin;
               curve->D1(First,PIn,VIn);
               curve->D1(Last,PFin,VFin);
               cout<<"New part one: "<<Pnt(PIn)<<" | "<<Pnt(PFin)<<endl;
               if (dP0.distance(Pnt(PIn)) < dP0.distance(Pnt(PFin)))
                  {
                  double delta_z = dP0(2) - PIn.Z();
                  dP = Point<3>(PIn.X()+delta_z*VIn.X()/VIn.Z(),PIn.Y()+delta_z*VIn.Y()/VIn.Z(),dP0(2));
                  V1 = VIn;
                  }
               else
                  {
                  double delta_z = dP0(2) - PFin.Z();
                  dP = Point<3>(PFin.X()+delta_z*VFin.X()/VFin.Z(),PIn.Y()+delta_z*VFin.Y()/VFin.Z(),dP0(2));
                  V1 = VFin;
                  }
               cout<<"New part two: "<<dP<<" | "<<V1.X()<<" "<<V1.Y()<<" "<<V1.Z()<<" | "<<dP0<<endl;
               }
            */
            //cout<<k<<"("<<i<<") ---> ("<<dP0<<") vs ("<<dP<<")"<<endl;
            working_map_points(3*i) = dP(0)-comp_dom.ref_points[3*i](0);
            working_map_points(3*i+1) = dP(1)-comp_dom.ref_points[3*i](1);
            working_map_points(3*i+2) = dP(2)-comp_dom.ref_points[3*i](2);
            comp_dom.edges_tangents[3*i] = V1.X();
            comp_dom.edges_tangents[3*i+1] = V1.Y();
            comp_dom.edges_tangents[3*i+2] = V1.Z();

            // here we compute the derivatives: we are basically computing the same residual on the tangent plane, which
            // will have the same derivatives as the surface
            gp_Pnt ref_point(dP(0),dP(1),dP(2));
            ref_point.Transform(reference_to_current_transformation.Inverted());
            gp_Vec curr_tangent(V1.X(),V1.Y(),V1.Z());
            gp_Dir ref_tangent_dir(curr_tangent);
            ref_tangent_dir.Transform(reference_to_current_transformation.Inverted());
            
            Point<3,fad_double> ref_point_pos(ref_point.X(),ref_point.Y(),ref_point.Z());
            Point<3,fad_double> rigid_lin_displ(s_x,s_y,s_z);
            Point<3,fad_double> target_point_pos(RotMatRow1*(ref_point_pos+(fad_double(-1.0))*ref_baricenter_pos),
                                                 RotMatRow2*(ref_point_pos+(fad_double(-1.0))*ref_baricenter_pos),
                                                 RotMatRow3*(ref_point_pos+(fad_double(-1.0))*ref_baricenter_pos));
            target_point_pos += baricenter_pos;
            
            Point<3,fad_double> ref_tangent(ref_tangent_dir.X(),ref_tangent_dir.Y(),ref_tangent_dir.Z());
            Point<3,fad_double> target_point_tangent(RotMatRow1*ref_tangent,
                                                     RotMatRow2*ref_tangent,
                                                     RotMatRow3*ref_tangent);
            
            Point<3,fad_double> guessed_point(comp_dom.support_points[i](0),comp_dom.support_points[i](1),comp_dom.support_points[i](2));

            fad_double residual_x = guessed_point(0) -
                                    target_point_tangent(0)/target_point_tangent(2)*(guessed_point(2)-target_point_pos(2)) - 
                                    target_point_pos(0);
            fad_double residual_y = guessed_point(1) -
                                    target_point_tangent(1)/target_point_tangent(2)*(guessed_point(2)-target_point_pos(2)) - 
                                    target_point_pos(1);

            
            //cout<<"oooyooo "<<comp_dom.map_points(3*i+1)-working_map_points(3*i+1)<<" vs "<<residual_y.val()<<" --> "
            //    <<comp_dom.map_points(3*i+1)-working_map_points(3*i+1)-residual_y.val()<<endl;

            jacobian_matrix.set(3*i,3*i,1.0);
            
            jacobian_matrix.set(3*i,3*i+2,-comp_dom.edges_tangents(3*i)/comp_dom.edges_tangents(3*i+2));
            
            for (unsigned int d=0; d<3;++d)
                jacobian_matrix.set(3*i,comp_dom.vector_dh.n_dofs()+
                                        comp_dom.dh.n_dofs()+
                                        comp_dom.dh.n_dofs()+
                                        3+d,residual_x.fastAccessDx(d));
            for (unsigned int d=0; d<4;++d)
                jacobian_matrix.set(3*i,comp_dom.vector_dh.n_dofs()+
                                        comp_dom.dh.n_dofs()+
                                        comp_dom.dh.n_dofs()+
                                        9+d,residual_x.fastAccessDx(3+d));
            jacobian_matrix.set(3*i+1,3*i+1,1.0);
            jacobian_matrix.set(3*i+1,3*i+2,-comp_dom.edges_tangents(3*i+1)/comp_dom.edges_tangents(3*i+2));
            for (unsigned int d=0; d<3;++d)
                jacobian_matrix.set(3*i+1,comp_dom.vector_dh.n_dofs()+
                                        comp_dom.dh.n_dofs()+
                                        comp_dom.dh.n_dofs()+
                                        3+d,residual_y.fastAccessDx(d));
            for (unsigned int d=0; d<4;++d)
                jacobian_matrix.set(3*i+1,comp_dom.vector_dh.n_dofs()+
                                        comp_dom.dh.n_dofs()+
                                        comp_dom.dh.n_dofs()+
                                        9+d,residual_y.fastAccessDx(3+d));
            //cout<<i<<" (point) "<<comp_dom.support_points[i]<<endl;
            //cout<<i<<" (edges_tangents) "<<comp_dom.edges_tangents(3*i)<<","<<comp_dom.edges_tangents(3*i+1)<<","<<comp_dom.edges_tangents(3*i+2)<<endl;   
            }              
            }

     // this cycle hooks the boat and far field double nodes
     // to their water twins that have been moved
     for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
         {
         if ( (comp_dom.flags[i] & water) &&
              (comp_dom.flags[i] & edge) )
            {
            fad_double s_x,s_y,s_z,v_x,v_y,v_z,s;

            s_x = hull_lin_displ(0);
            s_y = hull_lin_displ(1);
            s_z = hull_lin_displ(2);
            v_x = hull_quat_vector(0);
            v_y = hull_quat_vector(1);
            v_z = hull_quat_vector(2);
            s = hull_quat_scalar;

            s_x.diff(0,10);
            s_y.diff(1,10);
            s_z.diff(2,10);
            v_x.diff(3,10);
            v_y.diff(4,10);
            v_z.diff(5,10);
            s.diff(6,10);
            
            fad_double map_x,map_y,map_z;

            map_x = nodes_positions(3*i);
            map_y = nodes_positions(3*i+1);
            map_z = nodes_positions(3*i+2);
 
            map_x.diff(7,10);
            map_y.diff(8,10);
            map_z.diff(9,10);

            Point<3,fad_double> ref_baricenter_pos(comp_dom.boat_model.reference_hull_baricenter(0),
                                                   comp_dom.boat_model.reference_hull_baricenter(1),
                                                   comp_dom.boat_model.reference_hull_baricenter(2));

            Point<3,fad_double> baricenter_pos(s_x+comp_dom.boat_model.reference_hull_baricenter(0),
                                               s_y+comp_dom.boat_model.reference_hull_baricenter(1),
                                               s_z+comp_dom.boat_model.reference_hull_baricenter(2));

            Point<3,fad_double> RotMatRow1(1-2*v_y*v_y-2*v_z*v_z, 2*v_x*v_y-2*s*v_z, 2*v_x*v_z+2*s*v_y);
            Point<3,fad_double> RotMatRow2(2*v_x*v_y+2*s*v_z, 1-2*v_x*v_x-2*v_z*v_z, 2*v_y*v_z-2*s*v_x);
            Point<3,fad_double> RotMatRow3(2*v_x*v_z-2*s*v_y, 2*v_y*v_z+2*s*v_x, 1-2*v_y*v_y-2*v_x*v_x);

            gp_Pnt original_boat_mesh_point = Pnt(Point<3>(comp_dom.ref_points[3*i](0)+comp_dom.old_map_points(3*i),
                                                           comp_dom.ref_points[3*i](1)+comp_dom.old_map_points(3*i+1),
                                                           comp_dom.ref_points[3*i](2)+comp_dom.old_map_points(3*i+2)));
            gp_Pnt boat_mesh_point = original_boat_mesh_point;
            // we first take this point (which is in the RESTART hull location) and transform it to be in the
            // REFERENCE configuration 
            boat_mesh_point.Transform(restart_hull_location.Inverted());

            Point<3,fad_double> pp_orig_ref(boat_mesh_point.X(),boat_mesh_point.Y(),boat_mesh_point.Z());

            Point<3,fad_double> map_p_transf_ref(RotMatRow1*(pp_orig_ref+(fad_double(-1.0))*ref_baricenter_pos)+s_x-boat_mesh_point.X()+map_x,
                                                 RotMatRow2*(pp_orig_ref+(fad_double(-1.0))*ref_baricenter_pos)+s_y-boat_mesh_point.Y()+map_y,
                                                 RotMatRow3*(pp_orig_ref+(fad_double(-1.0))*ref_baricenter_pos)+s_z-boat_mesh_point.Z()+map_z);


            std::set<unsigned int> duplicates = comp_dom.double_nodes_set[i];
            duplicates.erase(i); 
            for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
                {//cout<<*pos<<"("<<i<<") "<<comp_dom.map_points(i)<<" vs "<<comp_dom.map_points(*pos)<<"   diff "<<comp_dom.map_points(i)-comp_dom.map_points(*pos)<<endl;
                for (unsigned int k=0; k<3; k++)
                    {
                    working_map_points(3*(*pos)+k) = comp_dom.map_points(3*i+k);
                
                    jacobian_matrix.add(3*(*pos)+k,3*(*pos)+k,1.0);
                if (comp_dom.flags[i] & transom_on_water)
                       {
                     
                       for (unsigned int d=0; d<3; d++)
                           jacobian_matrix.add(3*(*pos)+k,3*i+d,-map_p_transf_ref(k).fastAccessDx(7+d));
                       
                       for (unsigned int d=0; d<3; d++)
                           jacobian_matrix.add(3*(*pos)+k,comp_dom.vector_dh.n_dofs()+
                                                          comp_dom.dh.n_dofs()+
                                                          comp_dom.dh.n_dofs()+3+d,-map_p_transf_ref(k).fastAccessDx(d));
                        
                       for (unsigned int d=0; d<4; d++)
                           jacobian_matrix.add(3*(*pos)+k,comp_dom.vector_dh.n_dofs()+
                                                          comp_dom.dh.n_dofs()+
                                                          comp_dom.dh.n_dofs()+9+d,-map_p_transf_ref(k).fastAccessDx(d+3));
                       }
                    else
                       {
                       jacobian_matrix.add(3*(*pos)+k,3*i+k,-1.0);
                       }
                    }
                }
            }
         }


     nodes_pos_res = working_map_points;
     nodes_pos_res*=-1;
     nodes_pos_res.add(comp_dom.map_points);

     //cout<<"?? "<<nodes_pos_res(1674)<<" ("<<comp_dom.map_points(1674)<<" vs "<<working_map_points(1674)<<")"<<endl;



//////////////////////////////////////////////////////////////////////
  // here we take care of the bem part of the variables
////////////////////////////////////////////////////////////////////// 
     //cout<<"BEFORE "<<endl;
     //for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
     //    if (constraints.is_constrained(i))
     //       cout<<i<<" "<<src_yy(i+comp_dom.dh.n_dofs()+comp_dom.vector_dh.n_dofs())<<endl;

     bem_phi = (const Vector<double> &)phi;
     constraints.distribute(bem_phi);

     bem_dphi_dn = (const Vector<double> &)dphi_dn;
     constraints.distribute(bem_dphi_dn); 

     //cout<<"AFTER "<<endl;
     //for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
     //    if (constraints.is_constrained(i))
     //       cout<<i<<" "<<src_yy(i+comp_dom.dh.n_dofs()+comp_dom.vector_dh.n_dofs())<<endl;
     

     Vector<double> bem_bc(comp_dom.dh.n_dofs());
     for (unsigned int i=0; i<comp_dom.dh.n_dofs(); i++)
         {
         if ((comp_dom.flags[i] & water) ||
             (comp_dom.flags[i] & pressure) )
            bem_bc(i) = phi(i);
         else
            bem_bc(i) = dphi_dn(i);
         }

      // trying a fix for transom stern nodes
      for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
          {
          if ( comp_dom.flags[i] & transom_on_water )
             {
	     comp_dom.surface_nodes(i) = 0;
	     comp_dom.other_nodes(i) = 1;
             std::set<unsigned int> duplicates = comp_dom.double_nodes_set[i];
             duplicates.erase(i); 
             bem_bc(i) = 0;
             jacobian_matrix.add(i+comp_dom.dh.n_dofs()+comp_dom.vector_dh.n_dofs(),
                                 i+comp_dom.dh.n_dofs()+comp_dom.vector_dh.n_dofs(),
                                 -1.0);
             for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
                 {
                 bem_bc(i) += bem_dphi_dn(*pos)/duplicates.size();
                 jacobian_matrix.add(i+comp_dom.dh.n_dofs()+comp_dom.vector_dh.n_dofs(),
                                     *pos+comp_dom.dh.n_dofs()+comp_dom.vector_dh.n_dofs(),
                                     1.0/duplicates.size());
                 }
             bem_dphi_dn(i) = bem_bc(i);
             }
          }   
     // trying a fix for water/pressure nodes (pressure side)
      for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
          {
          if ( (comp_dom.flags[i] & pressure) &&  (comp_dom.flags[i] & near_water) )
             {
             std::set<unsigned int> duplicates = comp_dom.double_nodes_set[i];
             duplicates.erase(i);
             //cout<<i<<" "<<bem_bc(i)<<" "<<phi(i)<<" "<<dphi_dn(i)<<endl;
             unsigned int count=0;
             bem_bc(i) = 0;
             jacobian_matrix.add(i+comp_dom.vector_dh.n_dofs(),
                                 i+comp_dom.vector_dh.n_dofs(),
                                 -1.0);
             for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
                 if (comp_dom.flags[*pos] & water)
                    count++;
                 
             for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
                 {
                 if (comp_dom.flags[*pos] & water)
                    {
                    bem_bc(i) += phi(*pos)/count;
                    //cout<<*pos<<" ("<<i<<") "<<bem_bc(i)<<" "<<phi(*pos)<<" "<<dphi_dn(*pos)<<"  ? "<<count<<endl;
                    jacobian_matrix.add(i+comp_dom.vector_dh.n_dofs(),
                                        *pos+comp_dom.vector_dh.n_dofs(),
                                        1.0/count);
                    }
                 }
             bem_phi(i) = bem_bc(i);
             }
          }     


     if ( (sync_bem_with_geometry == true) || (new_time_step) )
        {
        //bem.assemble_system();
        //bem.residual(bem_residual, phi, dphi_dn);
        bem.solve(bem_phi, bem_dphi_dn, bem_bc);
        }
     else
        {
        //bem.residual(bem_residual, phi, dphi_dn);
        bem.solve_system(bem_phi, bem_dphi_dn, bem_bc);
        }


       // this is to enforce constraints in a more strict way
       // the vector given back by bem has constraints
       // imposed up to the bem GMRES tolerance
       for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
           {
           if ( constraints.is_constrained(i))
              {
              bem_phi(i) = 0;
              bem_dphi_dn(i) = 0;
              std::vector< std::pair< unsigned int, double > > entries = *constraints.get_constraint_entries(i);
              for (unsigned int k=0; k<entries.size(); ++k)
                  {
                  bem_phi(i) += entries[k].second*phi(entries[k].first);
                  bem_dphi_dn(i) += entries[k].second*dphi_dn(entries[k].first);
                  }
              }
           }

             
//////////////////////////////////////////////////////////////////////
  // here we compute eta_dry (just for a check)
  Vector<double> wind_valuee(dim);
  wind.vector_value(Point<3>(0.0,0.0,0.0),wind_valuee);
  Point<dim> Vinff;
  for (unsigned int i = 0; i < dim; i++)
      Vinff(i) = wind_valuee(i);
  double eta_dryy = 1.0;
  if (!comp_dom.no_boat && comp_dom.boat_model.is_transom && sqrt(Vinff*Vinff) > 0.0)
     {
     double transom_beamm = comp_dom.boat_model.CurrentPointRightTransom(1) - comp_dom.boat_model.CurrentPointLeftTransom(1);
     double transom_draftt = -comp_dom.boat_model.CurrentPointCenterTransom(2);
     double transom_aspect_ratioo = transom_beamm/transom_draftt;
     //cout<<restart_transom_right_tangent<<endl;
     //cout<<"b: "<<transom_beam.val()<<" d:"<<transom_draft.val()<<" "<<"AR: "<<transom_aspect_ratio.val()<<endl;
     double FrTr = sqrt(Vinff*Vinff)/sqrt(9.81*transom_draftt);
     double ReTr = sqrt(9.81*pow(transom_draftt,3.0))/1.307e-6;
     eta_dryy = 0.05*pow(FrTr,2.834)*pow(transom_aspect_ratioo,0.1352)*pow(ReTr,0.01338);
     if (eta_dryy > 1.0)
        eta_dryy = 1.0;

     cout<<"Current eta_dry: "<<eta_dryy<<endl;
     }

//////////////////////////////////////////////////////////////////////
  // here we take care of the free surface boundary condition (differential)
  // part of the variables, and the boat neumann condition variables
////////////////////////////////////////////////////////////////////// 



// building reference cell
  Triangulation<2> ref_triangulation;

  std::vector<Point<2> > ref_vertices;
  std::vector<CellData<2> > ref_cells;
  SubCellData ref_subcelldata;

  ref_vertices.resize(4);
  ref_vertices[0](0)= 0.0; ref_vertices[0](1)= 0.0; //ref_vertices[0](2)=0.0;
  ref_vertices[1](0)= 1.0; ref_vertices[1](1)= 0.0; //ref_vertices[1](2)=0.0;
  ref_vertices[2](0)= 0.0; ref_vertices[2](1)= 1.0; //ref_vertices[2](2)=0.0;
  ref_vertices[3](0)= 1.0; ref_vertices[3](1)= 1.0; //ref_vertices[3](2)=0.0;

  ref_cells.resize(1);

  ref_cells[0].vertices[0]=0; ref_cells[0].vertices[1]=1; ref_cells[0].vertices[2]=3; ref_cells[0].vertices[3]=2;
  ref_cells[0].material_id = 1;

  GridTools::delete_unused_vertices (ref_vertices, ref_cells, ref_subcelldata);
  GridReordering<2>::reorder_cells (ref_cells);

  ref_triangulation.create_triangulation_compatibility(ref_vertices, ref_cells, ref_subcelldata );

  FE_Q<2> fe(1);
 

  DoFHandler<2> ref_dh(ref_triangulation);
  ref_dh.distribute_dofs(fe);

  FEValues<2> ref_fe_v(StaticMappingQ1<2>::mapping, fe, *comp_dom.quadrature,
   		         update_values | update_gradients |
		         update_quadrature_points |
		         update_JxW_values);

  const unsigned int n_q_points = ref_fe_v.n_quadrature_points;
  const unsigned int  dofs_per_cell   = fe.dofs_per_cell;

  //typedef typename DoFHandler<dim-1,dim>::active_cell_iterator cell_it;
  //typedef typename Triangulation<dim-1,dim>::active_cell_iterator tria_it;  
  DoFHandler<2>::active_cell_iterator ref_cell = ref_dh.begin_active();
  ref_fe_v.reinit(ref_cell);  

///////////////////////////////////

// test: let's try assemble actual mass matrix
  double g = 9.81;
  DphiDt_sys_matrix = 0;
  DphiDt_sys_solution = 0;
  DphiDt_sys_solution_2 = 0;
  DphiDt_sys_solution_3 = 0;
  DphiDt_sys_rhs = 0;
  DphiDt_sys_rhs_2 = 0;
  DphiDt_sys_rhs_3 = 0;
  DphiDt_sys_rhs_4 = 0;


//let's build the residual on the free surface cells (differential components)
  std::vector<double> eta_res(comp_dom.dh.n_dofs(),0.0);
  std::vector<double> phi_res(comp_dom.dh.n_dofs(),0.0);
  std::vector<double> dphi_dn_res(comp_dom.dh.n_dofs(),0.0);
  std::vector<double> x_smoothing_res(comp_dom.dh.n_dofs(),0.0);
  std::vector<double> y_smoothing_res(comp_dom.dh.n_dofs(),0.0);
// we'll need Vinf
  Vector<double> wind_value(dim);
  wind.vector_value(Point<3>(0.0,0.0,0.0),wind_value);
  Point<dim> Vinf;
  for (unsigned int i = 0; i < dim; i++)
      Vinf(i) = wind_value(i);

  double amplitude = 0.02;
  double rho = 1025.1;
  double ref_height;
  double Fn;
  if (comp_dom.no_boat)
     {
     Fn = Vinf(0)/sqrt(max_z_coor_value*g);
     ref_height = amplitude;
     }
  else
     {
     ref_height = fabs(comp_dom.boat_model.PointMidBot(2));
     Fn = Vinf(0)/sqrt(comp_dom.boat_model.boatWetLength*g);
     }

  //this is to assign correct values to the (possible) exact pressure to be used in pressure boundary condition
  Vector <double> ex_press(comp_dom.dh.n_dofs());
  for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
      {
      if (comp_dom.flags[i] & pressure && !(comp_dom.flags[i] & near_water) )
         {
         double k=0.62994; double omega=2.4835; double h=5.5;  double a=0.00;
         double time_factor = 1.0;
         double time_factor_deriv = 0.0;
         double ramp_length = 20.0;
         if (t<ramp_length)
            {
            time_factor = 0.5*sin(3.141592654*(t)/ramp_length-3.141592654/2)+0.5;
            time_factor_deriv = 0.5*3.141592654/ramp_length*cos(3.141592654*(t)/ramp_length-3.141592654/2);
            }                  
         double dphi_dt = omega*omega*a/k*cosh(k*(comp_dom.support_points[i](2)+h))/sinh(k*h)*cos(k*comp_dom.support_points[i](0)+omega*t)*time_factor +
                          omega*a/k*cosh(k*(comp_dom.support_points[i](2)+h))/sinh(k*h)*sin(k*comp_dom.support_points[i](0)+omega*t)*time_factor_deriv;
         Point<3> grad_phi(omega*a*cosh(k*(comp_dom.support_points[i](2)+h))/sinh(k*h)*cos(k*comp_dom.support_points[i](0)+omega*t)*time_factor,
                           0.0*time_factor,
                           omega*a*sinh(k*(comp_dom.support_points[i](2)+h))/sinh(k*h)*sin(k*comp_dom.support_points[i](0)+omega*t)*time_factor);



         ex_press(i) = ( -dphi_dt - 0.5*(grad_phi*grad_phi) - 
                         grad_phi(0)*instantWindValue(0)-grad_phi(1)*instantWindValue(1)-grad_phi(2)*instantWindValue(2) -
                         comp_dom.support_points[i](2)*g ) * rho;
         }
      }


  FullMatrix<double>   local_DphiDt_matrix (dofs_per_cell, dofs_per_cell);
  FullMatrix<double>   local_DphiDt_matrix_2 (dofs_per_cell, dofs_per_cell);
  Vector<double>       local_DphiDt_rhs (dofs_per_cell);
  Vector<double>       local_DphiDt_rhs_2 (dofs_per_cell);
  Vector<double>       local_DphiDt_rhs_3 (dofs_per_cell);
  Vector<double>       local_DphiDt_rhs_4 (dofs_per_cell);

  cell_it
  cell = comp_dom.dh.begin_active(),
  endc = comp_dom.dh.end();

  std::vector<fad_double> ref_coors(3*dofs_per_cell);
  std::vector<fad_double> coors(3*dofs_per_cell);
  std::vector<fad_double> phis(dofs_per_cell);
  std::vector<fad_double> dphi_dns(dofs_per_cell);
  std::vector<fad_double> ref_coors_dot(3*dofs_per_cell);
  std::vector<fad_double> coors_dot(3*dofs_per_cell);
  std::vector<fad_double> phis_dot(dofs_per_cell);
  std::vector<fad_double> dphi_dns_dot(dofs_per_cell);
  std::vector<fad_double> x_displs(dofs_per_cell);
  std::vector<fad_double> y_displs(dofs_per_cell);

  std::vector<fad_double> hull_baricenter_displ(3);
  std::vector<fad_double> hull_baricenter_vel(3);
  std::vector<fad_double> hull_quaternion(4);
  std::vector<fad_double> hull_omega(3); 
 
  std::vector<fad_double> loc_eta_res(dofs_per_cell);
  std::vector<fad_double> loc_phi_res(dofs_per_cell);
  std::vector<fad_double> loc_dphi_dn_res(dofs_per_cell);
  std::vector<fad_double> loc_x_smooth_res(dofs_per_cell);
  std::vector<fad_double> loc_y_smooth_res(dofs_per_cell);
  Point<3,fad_double> loc_pressure_force;
  Point<3,fad_double> loc_pressure_moment;
  Point<3> pressure_force(0.0,0.0,0.0);
  Point<3> pressure_moment(0.0,0.0,0.0);

  std::vector< std::vector<fad_double> > loc_stiffness_matrix(dofs_per_cell);
  std::vector< std::vector<fad_double> > loc_mass_matrix(dofs_per_cell);
  std::vector< std::vector<fad_double> > loc_supg_mass_matrix(dofs_per_cell);

  for (unsigned int i=0; i<dofs_per_cell; ++i)
      {
      loc_stiffness_matrix[i].resize(dofs_per_cell);
      loc_mass_matrix[i].resize(dofs_per_cell);
      loc_supg_mass_matrix[i].resize(dofs_per_cell);
      }

  std::vector<Point<dim> > initial_support_points(comp_dom.dh.n_dofs());
  for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
      {
      initial_support_points[i] = Point<3>(comp_dom.ref_points[3*i](0)+comp_dom.initial_map_points(3*i),
                                           comp_dom.ref_points[3*i](1)+comp_dom.initial_map_points(3*i+1),
                                           comp_dom.ref_points[3*i](2)+comp_dom.initial_map_points(3*i+2));
      } 


  std::vector<unsigned int> local_dof_indices(dofs_per_cell);

   
  for (; cell!=endc; ++cell)
      {

      //if (cell->material_id() == comp_dom.free_sur_ID1 ||
      //    cell->material_id() == comp_dom.free_sur_ID2 ||
      //    cell->material_id() == comp_dom.free_sur_ID3 ||
      //    cell->material_id() == comp_dom.wall_sur_ID1 ||
      //    cell->material_id() == comp_dom.wall_sur_ID2 ||
      //    cell->material_id() == comp_dom.wall_sur_ID3 )
         {
         //if (fabs(t-0.005) < 1e-4)
         //std::cout<<std::endl;
         //std::cout<<"Cell: "<<cell<<std::endl;

         local_DphiDt_matrix = 0;
         local_DphiDt_matrix_2 = 0;
         local_DphiDt_rhs = 0;
         local_DphiDt_rhs_2 = 0;
         local_DphiDt_rhs_3 = 0;
         local_DphiDt_rhs_4 = 0;

         loc_pressure_force = Point<3,fad_double>(0.0,0.0,0.0);
         loc_pressure_moment = Point<3,fad_double>(0.0,0.0,0.0);

         cell->get_dof_indices(local_dof_indices);
         for (unsigned int i=0; i<dofs_per_cell; ++i)
             {
             loc_eta_res[i] = 0;
             loc_phi_res[i] = 0;
             loc_dphi_dn_res[i] = 0;
             loc_x_smooth_res[i] = 0;
             loc_y_smooth_res[i] = 0;
             for (unsigned int j=0; j<dofs_per_cell; ++j)
                 {
                 loc_mass_matrix[i][j] = 0;
                 loc_supg_mass_matrix[i][j] = 0;
                 loc_stiffness_matrix[i][j] = 0;
                 }
             //cout<<local_dof_indices[i]<<" ";
             }
            /*
             for (unsigned int i=0; i<dofs_per_cell; ++i)
                {
                if (local_dof_indices[i]== 44)
                   {
                   cout<<"44!!: "<<cell<<endl;
                   for (unsigned int k=0; k<dofs_per_cell; ++k)
                       {
                       for (unsigned int j=0; j<3; ++j)
                           {
                           cout<<3*local_dof_indices[k]+j<<" ";
                           }
                       cout<<local_dof_indices[k]+comp_dom.vector_dh.n_dofs()<<" ";
                       cout<<local_dof_indices[k]+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()<<" ";
                       }
                   cout<<endl;
                   }
                }
            */

         for (unsigned int j=0; j<3; ++j)
             {
             hull_baricenter_vel[j] = hull_lin_vel(j);
             hull_baricenter_displ[j] = hull_lin_displ(j);
             hull_omega[j] = hull_ang_vel(j);
             hull_quaternion[j] = hull_quat_vector(j);
             hull_baricenter_vel[j].diff(10*dofs_per_cell+j,10*dofs_per_cell+13);
             hull_baricenter_displ[j].diff(10*dofs_per_cell+3+j,10*dofs_per_cell+13);
             hull_omega[j].diff(10*dofs_per_cell+6+j,10*dofs_per_cell+13);
             hull_quaternion[j].diff(10*dofs_per_cell+9+j,10*dofs_per_cell+13);
             }
         hull_quaternion[3] = hull_quat_scalar;
         hull_quaternion[3].diff(10*dofs_per_cell+12,10*dofs_per_cell+13);

         Point<3,fad_double> ref_baricenter_pos(comp_dom.boat_model.reference_hull_baricenter(0),
                                                comp_dom.boat_model.reference_hull_baricenter(1),
                                                comp_dom.boat_model.reference_hull_baricenter(2));

         Point<3,fad_double> baricenter_pos(hull_baricenter_displ[0]+comp_dom.boat_model.reference_hull_baricenter(0),
                                            hull_baricenter_displ[1]+comp_dom.boat_model.reference_hull_baricenter(1),
                                            hull_baricenter_displ[2]+comp_dom.boat_model.reference_hull_baricenter(2));

         Point<3,fad_double> RotMatRow1(1-2*hull_quaternion[1]*hull_quaternion[1]-2*hull_quaternion[2]*hull_quaternion[2],
                                        2*hull_quaternion[0]*hull_quaternion[1]-2*hull_quaternion[3]*hull_quaternion[2],
                                        2*hull_quaternion[0]*hull_quaternion[2]+2*hull_quaternion[3]*hull_quaternion[1]);
         Point<3,fad_double> RotMatRow2(2*hull_quaternion[0]*hull_quaternion[1]+2*hull_quaternion[3]*hull_quaternion[2],
                                        1-2*hull_quaternion[0]*hull_quaternion[0]-2*hull_quaternion[2]*hull_quaternion[2],
                                        2*hull_quaternion[1]*hull_quaternion[2]-2*hull_quaternion[3]*hull_quaternion[0]);
         Point<3,fad_double> RotMatRow3(2*hull_quaternion[0]*hull_quaternion[2]-2*hull_quaternion[3]*hull_quaternion[1],
                                        2*hull_quaternion[1]*hull_quaternion[2]+2*hull_quaternion[3]*hull_quaternion[0],
                                        1-2*hull_quaternion[1]*hull_quaternion[1]-2*hull_quaternion[0]*hull_quaternion[0]);


         for (unsigned int i=0; i<dofs_per_cell; ++i)
             if ( ((comp_dom.flags[local_dof_indices[i]] & boat) &&
                  !(comp_dom.flags[local_dof_indices[i]] & near_water) ) ||
                  (comp_dom.flags[local_dof_indices[i]] & transom_on_water) )
                {
                Point<3> reference_point_position(comp_dom.ref_points[3*local_dof_indices[i]](0)+
                                                  comp_dom.map_points(3*local_dof_indices[i])-
                                                  comp_dom.rigid_motion_map_points(3*local_dof_indices[i]),
                                                  comp_dom.ref_points[3*local_dof_indices[i]](1)+
                                                  comp_dom.map_points(3*local_dof_indices[i]+1)-
                                                  comp_dom.rigid_motion_map_points(3*local_dof_indices[i]+1),
                                                  comp_dom.ref_points[3*local_dof_indices[i]](2)+
                                                  comp_dom.map_points(3*local_dof_indices[i]+2)-
                                                  comp_dom.rigid_motion_map_points(3*local_dof_indices[i]+2));
                
                for (unsigned int j=0; j<3; ++j)
                    {
                    ref_coors[3*i+j] = reference_point_position(j);
                    (ref_coors[3*i+j]).diff(3*i+j,10*dofs_per_cell+13);
                    
                    //cout<<"TESTTTTTT:::: "<<coors[3*i+j]-comp_dom.support_points[local_dof_indices[i]](j)<<endl;
                    ref_coors_dot[3*i+j] = nodes_velocities(3*local_dof_indices[i]+j);
                    (ref_coors_dot[3*i+j]).diff(3*i+j+5*dofs_per_cell,10*dofs_per_cell+13);
                    }
                 Point<3,fad_double> sacado_orig_point(ref_coors[3*i],ref_coors[3*i+1],ref_coors[3*i+2]);
                 Point<3,fad_double> sacado_curr_point(RotMatRow1*(sacado_orig_point+(fad_double(-1.0))*ref_baricenter_pos)+
                                                       baricenter_pos(0),
                                                       RotMatRow2*(sacado_orig_point+(fad_double(-1.0))*ref_baricenter_pos)+
                                                       baricenter_pos(1),
                                                       RotMatRow3*(sacado_orig_point+(fad_double(-1.0))*ref_baricenter_pos)+
                                                       baricenter_pos(2));
                 //cout<<"Orig: "<<sacado_orig_point<<endl;//(0)<<" "<<sacado_orig_point(1).val()<<" "<<sacado_orig_point(2).val()<<endl;
                 //cout<<"Point "<<local_dof_indices[i]<<endl;
                 //cout<<"Curr sacado: "<<sacado_curr_point(0).val()<<" "<<sacado_curr_point(1).val()<<" "<<sacado_curr_point(2).val()<<endl;
                 //cout<<"Curr true: "<<comp_dom.ref_points[3*local_dof_indices[i]](0)+comp_dom.map_points(3*local_dof_indices[i])+comp_dom.rigid_motion_map_points(3*local_dof_indices[i])<<" "
                 //                   <<comp_dom.ref_points[3*local_dof_indices[i]](1)+comp_dom.map_points(3*local_dof_indices[i]+1)+comp_dom.rigid_motion_map_points(3*local_dof_indices[i]+1)<<" "
                //                    <<comp_dom.ref_points[3*local_dof_indices[i]](2)+comp_dom.map_points(3*local_dof_indices[i]+2)+comp_dom.rigid_motion_map_points(3*local_dof_indices[i]+2)<<endl;
     //cout<<"Mismatch: "<<comp_dom.ref_points[3*local_dof_indices[i]](0)+comp_dom.map_points(3*local_dof_indices[i])-sacado_curr_point(0).val()<<" "
     //             <<comp_dom.ref_points[3*local_dof_indices[i]](1)+comp_dom.map_points(3*local_dof_indices[i]+1)-sacado_curr_point(1).val()<<" "
     //             <<comp_dom.ref_points[3*local_dof_indices[i]](2)+comp_dom.map_points(3*local_dof_indices[i]+2)-sacado_curr_point(2).val()<<endl;

                 gp_Pnt original_target_boat_mesh_point = Pnt(Point<3>(comp_dom.ref_points[3*local_dof_indices[i]](0)+comp_dom.old_map_points(3*i),
                                                           comp_dom.ref_points[3*local_dof_indices[i]](1)+comp_dom.old_map_points(3*i+1),
                                                           comp_dom.ref_points[3*local_dof_indices[i]](2)+comp_dom.old_map_points(3*i+2)));
                 gp_Pnt target_boat_mesh_point = original_target_boat_mesh_point;
                 // we first take this point (which is in the RESTART hull location) and transform it to be in the
                 // REFERENCE configuration 
                 target_boat_mesh_point.Transform(restart_hull_location.Inverted());
                 Point<3,fad_double> target_ref_point(target_boat_mesh_point.X(),target_boat_mesh_point.Y(),target_boat_mesh_point.Z());
                 Point<3,fad_double> sacado_target_point(RotMatRow1*(target_ref_point+(fad_double(-1.0))*ref_baricenter_pos)+
                                                         baricenter_pos(0),
                                                         RotMatRow2*(target_ref_point+(fad_double(-1.0))*ref_baricenter_pos)+
                                                         baricenter_pos(1),
                                                         RotMatRow3*(target_ref_point+(fad_double(-1.0))*ref_baricenter_pos)+
                                                         baricenter_pos(2));

                 //questa velocità RIGIDA va preparata utilizzando i TARGET POINTS, altrimenti contiene due volte la derivata rispetto
                 //a map_points.rigid_motions_map_points
                 Point<3,fad_double> sacado_non_rig_vel(ref_coors_dot[3*i],ref_coors_dot[3*i+1],ref_coors_dot[3*i+2]);
                 Point<3,fad_double> sacado_rig_vel(hull_omega[1]*(sacado_target_point[2]-baricenter_pos(2))-
                                                    hull_omega[2]*(sacado_target_point[1]-baricenter_pos(1))+hull_baricenter_vel[0],
                                                    hull_omega[2]*(sacado_target_point[0]-baricenter_pos(0))-
                                                    hull_omega[0]*(sacado_target_point[2]-baricenter_pos(2))+hull_baricenter_vel[1],
                                                    hull_omega[0]*(sacado_target_point[1]-baricenter_pos(1))-
                                                    hull_omega[1]*(sacado_target_point[0]-baricenter_pos(0))+hull_baricenter_vel[2]);
                 
                 //if (comp_dom.ref_points[3*i].distance(Point<3>(-2.388065,0.0,-0.355460)) < 0.005)
                 //   {
                 //   cout<<"Non Rigid Vel: "<<sacado_non_rig_vel(0).val()<<" "<<sacado_non_rig_vel(1).val()<<" "<<sacado_non_rig_vel(2).val()<<endl;
                 //   cout<<"Rigid Vel: "<<sacado_rig_vel(0).val()<<" "<<sacado_rig_vel(1).val()<<" "<<sacado_rig_vel(2).val()<<endl;
                 //   }
                 for (unsigned int j=0; j<3; ++j)
                    {
                    coors[3*i+j] = sacado_curr_point(j);
                    coors_dot[3*i+j] = sacado_non_rig_vel(j)+sacado_rig_vel(j);
                    }
                }
             else
                for (unsigned int j=0; j<3; ++j)
                    {
                    coors[3*i+j] = comp_dom.support_points[local_dof_indices[i]](j);
                    (coors[3*i+j]).diff(3*i+j,10*dofs_per_cell+13);
                    coors_dot[3*i+j] = nodes_velocities(3*local_dof_indices[i]+j);
                    (coors_dot[3*i+j]).diff(3*i+j+5*dofs_per_cell,10*dofs_per_cell+13);
                    }
         for (unsigned int i=0; i<dofs_per_cell; ++i)
             {
             phis[i] = phi(local_dof_indices[i]);
             phis[i].diff(i+3*dofs_per_cell,10*dofs_per_cell+13);
             phis_dot[i] = phi_time_derivs(local_dof_indices[i]);
             phis_dot[i].diff(i+8*dofs_per_cell,10*dofs_per_cell+13);
             dphi_dns[i] = dphi_dn(local_dof_indices[i]);
             dphi_dns[i].diff(i+4*dofs_per_cell,10*dofs_per_cell+13);
             dphi_dns_dot[i] = dphi_dn_time_derivs(local_dof_indices[i]);
             dphi_dns_dot[i].diff(i+9*dofs_per_cell,10*dofs_per_cell+13);
             //std::cout<<i<<"--> "<<local_dof_indices[i]<<"--------->"<<bem_phi(local_dof_indices[i])<<"  "<<bem_dphi_dn(local_dof_indices[i])<<endl;
             }

         //just in case of additional pressure to be added past wet transom stern
         fad_double eta_dry = 1.0;
         if (!comp_dom.no_boat && comp_dom.boat_model.is_transom && sqrt(Vinf*Vinf) > 0.0)
            {
            gp_Pnt ref_left_transom_pnt = Pnt(restart_transom_left_point);
            gp_Pnt ref_right_transom_pnt = Pnt(restart_transom_right_point);
            gp_Pnt ref_center_transom_pnt = Pnt(restart_transom_center_point);
            // we first take these points (which is in the RESTART hull location) and transform it to be in the
                // REFERENCE configuration 
            ref_left_transom_pnt.Transform(restart_hull_location.Inverted());
            ref_right_transom_pnt.Transform(restart_hull_location.Inverted());
            ref_center_transom_pnt.Transform(restart_hull_location.Inverted());
            Point<3,fad_double> ref_left_transom_point(ref_left_transom_pnt.X(),ref_left_transom_pnt.Y(),ref_left_transom_pnt.Z());
            Point<3,fad_double> ref_right_transom_point(ref_right_transom_pnt.X(),ref_right_transom_pnt.Y(),ref_right_transom_pnt.Z());
            Point<3,fad_double> ref_center_transom_point(ref_center_transom_pnt.X(),ref_center_transom_pnt.Y(),ref_center_transom_pnt.Z());
            Point<3,fad_double> curr_left_transom_point(baricenter_pos(0)+
                                                        RotMatRow1*(ref_left_transom_point+(fad_double(-1.0))*ref_baricenter_pos),
                                                        baricenter_pos(1)+
                                                        RotMatRow2*(ref_left_transom_point+(fad_double(-1.0))*ref_baricenter_pos),
                                                        baricenter_pos(2)+
                                                        RotMatRow3*(ref_left_transom_point+(fad_double(-1.0))*ref_baricenter_pos));
            Point<3,fad_double> curr_right_transom_point(baricenter_pos(0)+
                                                         RotMatRow1*(ref_right_transom_point+(fad_double(-1.0))*ref_baricenter_pos),
                                                         baricenter_pos(1)+
                                                         RotMatRow2*(ref_right_transom_point+(fad_double(-1.0))*ref_baricenter_pos),
                                                         baricenter_pos(2)+
                                                         RotMatRow3*(ref_right_transom_point+(fad_double(-1.0))*ref_baricenter_pos));
            Point<3,fad_double> curr_center_transom_point(baricenter_pos(0)+
                                                          RotMatRow1*(ref_center_transom_point+(fad_double(-1.0))*ref_baricenter_pos),
                                                          baricenter_pos(1)+
                                                          RotMatRow2*(ref_center_transom_point+(fad_double(-1.0))*ref_baricenter_pos),
                                                          baricenter_pos(2)+
                                                          RotMatRow3*(ref_center_transom_point+(fad_double(-1.0))*ref_baricenter_pos));
            fad_double transom_beam = restart_transom_right_point(1) - restart_transom_left_point(1) +
                    restart_transom_right_tangent(1)/restart_transom_right_tangent(2)*(curr_right_transom_point(2)-restart_transom_right_point(2))-
                    restart_transom_left_tangent(1)/restart_transom_left_tangent(2)*(curr_left_transom_point(2)-restart_transom_left_point(2));
            fad_double transom_draft = -curr_center_transom_point(2);
            fad_double transom_aspect_ratio = transom_beam/transom_draft;
            //cout<<restart_transom_right_tangent<<endl;
            //cout<<"b: "<<transom_beam.val()<<" d:"<<transom_draft.val()<<" "<<"AR: "<<transom_aspect_ratio.val()<<endl;
            fad_double FrT = sqrt(Vinf*Vinf)/sqrt(9.81*transom_draft);
            fad_double ReT = sqrt(9.81*pow(transom_draft,3.0))/1.307e-6;
            eta_dry = 0.05*pow(FrT,2.834)*pow(transom_aspect_ratio,0.1352)*pow(ReT,0.01338);
            if (eta_dry.val() > 1)
               eta_dry = 1.0;
            if (cell == comp_dom.dh.begin_active())
               cout<<"Sacado eta_dry: "<<eta_dry.val()<<endl;
            }

         // computation of displacements
         for (unsigned int i=0; i<dofs_per_cell; ++i)
             {
             x_displs[i] = coors[3*i] - comp_dom.ref_points[3*local_dof_indices[i]](0);
             y_displs[i] = coors[3*i+1] - comp_dom.ref_points[3*local_dof_indices[i]](1);
             //cout<<i<<" "<<coors[3*i].val()<<" "<<comp_dom.ref_points[3*local_dof_indices[i]](0)<<" "<<comp_dom.support_points[local_dof_indices[i]](0)<<endl;
             //cout<<i<<" "<<coors[3*i+1].val()<<" "<<comp_dom.ref_points[3*local_dof_indices[i]](1)<<" "<<comp_dom.support_points[local_dof_indices[i]](1)<<endl;
             }
         // computation of cell center
         Point<3,fad_double> center(0.0,0.0,0.0);
         for (unsigned int i=0; i<dofs_per_cell; ++i)
             {
             center += (Point<3,fad_double>(coors[3*i],coors[3*i+1],coors[3*i+2]))/dofs_per_cell;
             }
         // computation of cell diameter
         fad_double cell_diameter = 0.0;
         for (unsigned int i=0; i<dofs_per_cell; ++i)
             {
             fad_double dof_center_distance = sqrt( (center(0)-coors[3*i])*(center(0)-coors[3*i])+
                                                    (center(1)-coors[3*i+1])*(center(1)-coors[3*i+1])+
                                                    (center(2)-coors[3*i+2])*(center(2)-coors[3*i+2]) );
             cell_diameter += 2.0*(dof_center_distance)/dofs_per_cell;
             }
         
         Point<3,fad_double> gg(fad_double(0.0),fad_double(0.0),fad_double(g));

         std::vector<fad_double> eta_dot_rhs_fun(n_q_points);
         std::vector<fad_double> phi_dot_rhs_fun(n_q_points);
         std::vector< Point<3,fad_double> > fluid_vel(n_q_points);
         std::vector<fad_double> q_JxW(n_q_points);
   
   
         for (unsigned int q=0; q<n_q_points; ++q)
             {
             Point<3,fad_double> q_point(fad_double(0.0),fad_double(0.0),fad_double(0.0));
             Point<3,fad_double> u_deriv_pos(fad_double(0.0),fad_double(0.0),fad_double(0.0));
             Point<3,fad_double> v_deriv_pos(fad_double(0.0),fad_double(0.0),fad_double(0.0));
             fad_double u_deriv_phi = fad_double(0.0);
             fad_double v_deriv_phi = fad_double(0.0);
             fad_double q_dphi_dn = fad_double(0.0);
             fad_double q_phi_dot = fad_double(0.0);
             fad_double q_x_dot = fad_double(0.0);
             fad_double q_y_dot = fad_double(0.0);
             fad_double q_z_dot = fad_double(0.0);
             fad_double q_eta = fad_double(0.0);
             fad_double q_ex_press = fad_double(0.0);
             Point<3,fad_double> q_init(fad_double(0.0),fad_double(0.0),fad_double(0.0));
             cout.precision(10);
             for (unsigned int i=0; i<dofs_per_cell; ++i)
                 {
                 unsigned int index = local_dof_indices[i];
                 q_point += fad_double(ref_fe_v.shape_value(i,q))*Point<3,fad_double>(coors[3*i],coors[3*i+1],coors[3*i+2]);
                 u_deriv_pos += fad_double(ref_fe_v.shape_grad(i,q)[0])*Point<3,fad_double>(coors[3*i],coors[3*i+1],coors[3*i+2]);
                 v_deriv_pos += fad_double(ref_fe_v.shape_grad(i,q)[1])*Point<3,fad_double>(coors[3*i],coors[3*i+1],coors[3*i+2]);
                 u_deriv_phi += fad_double(ref_fe_v.shape_grad(i,q)[0])*phis[i];
                 v_deriv_phi += fad_double(ref_fe_v.shape_grad(i,q)[1])*phis[i];
                 q_dphi_dn += fad_double(ref_fe_v.shape_value(i,q))*dphi_dns[i];
                 q_phi_dot += fad_double(ref_fe_v.shape_value(i,q))*phis_dot[i];
                 q_x_dot += fad_double(ref_fe_v.shape_value(i,q))*coors_dot[3*i];
                 q_y_dot += fad_double(ref_fe_v.shape_value(i,q))*coors_dot[3*i+1];
                 q_z_dot += fad_double(ref_fe_v.shape_value(i,q))*coors_dot[3*i+2];
                 q_eta +=  fad_double(ref_fe_v.shape_value(i,q))*coors[3*i+2];
                 q_ex_press += fad_double(ref_fe_v.shape_value(i,q))*ex_press[local_dof_indices[i]];
                 q_init(0) += ref_fe_v.shape_value(i,q)*fad_double(initial_support_points[local_dof_indices[i]](0));
                 q_init(1) += ref_fe_v.shape_value(i,q)*fad_double(initial_support_points[local_dof_indices[i]](1));
                 q_init(2) += ref_fe_v.shape_value(i,q)*fad_double(initial_support_points[local_dof_indices[i]](2));
                 //q_init += ref_fe_v.shape_value(i,q)*comp_dom.ref_points[3*local_dof_indices[i]];
                 //std::cout<<i<<"-------> "<<u_deriv_pos<<" "<<v_deriv_pos<<" "<<u_deriv_phi<<" "<<v_deriv_phi<<endl;
                 //std::cout<<i<<"-------> "<<coors[3*i]<<" "<<coors[3*i+1]<<" "<<coors[3*i+2]<<" "<<endl;
                 //cout<<"Earlier: "<<"i "<<i<<"   q "<<q<<"  "<<ref_fe_v.shape_value(i,q)*gg(2).val()<<endl;
                 }
             Point<3> cen = cell->center();
             if (fabs(cen(2)) > 1e-5 )
                {
                q_init(0) += (hull_baricenter_displ[0]-fad_double(restart_hull_displacement(0)));
                q_init(1) += (hull_baricenter_displ[1]-fad_double(restart_hull_displacement(1)));
                q_init(2) += (hull_baricenter_displ[2]-fad_double(restart_hull_displacement(2)));
                }
             fad_double transom_added_pressure = 0.0;

             Point<3,fad_double> q_normal(u_deriv_pos(1)*v_deriv_pos(2)-u_deriv_pos(2)*v_deriv_pos(1),
                                          u_deriv_pos(2)*v_deriv_pos(0)-u_deriv_pos(0)*v_deriv_pos(2),
                                          u_deriv_pos(0)*v_deriv_pos(1)-u_deriv_pos(1)*v_deriv_pos(0));
             //std::cout<<"q_normal="<<q_normal<<std::endl;
             //std::cout<<"q_y_dot="<<q_y_dot<<std::endl;
             fad_double q_jac_det = q_normal.norm();
             q_normal/=q_jac_det;
             fad_double a = 1.0/((u_deriv_pos*u_deriv_pos)*(v_deriv_pos*v_deriv_pos)-(u_deriv_pos*v_deriv_pos)*(u_deriv_pos*v_deriv_pos));
             fad_double d11 = a*(u_deriv_pos(0)*v_deriv_pos*v_deriv_pos-v_deriv_pos(0)*u_deriv_pos*v_deriv_pos);
             fad_double d21 = a*(u_deriv_pos(1)*v_deriv_pos*v_deriv_pos-v_deriv_pos(1)*u_deriv_pos*v_deriv_pos);
             fad_double d31 = a*(u_deriv_pos(2)*v_deriv_pos*v_deriv_pos-v_deriv_pos(2)*u_deriv_pos*v_deriv_pos);
             fad_double d12 = a*(v_deriv_pos(0)*u_deriv_pos*u_deriv_pos-u_deriv_pos(0)*u_deriv_pos*v_deriv_pos);
             fad_double d22 = a*(v_deriv_pos(1)*u_deriv_pos*u_deriv_pos-u_deriv_pos(1)*u_deriv_pos*v_deriv_pos);
             fad_double d32 = a*(v_deriv_pos(2)*u_deriv_pos*u_deriv_pos-u_deriv_pos(2)*u_deriv_pos*v_deriv_pos);
             Point<3,fad_double> phi_surf_grad(d11*u_deriv_phi+d12*v_deriv_phi,
                                               d21*u_deriv_phi+d22*v_deriv_phi,
                                               d31*u_deriv_phi+d32*v_deriv_phi);
             Point<3,fad_double> phi_surf_grad_corrected(phi_surf_grad(0) - phi_surf_grad(2)*q_normal(0)/q_normal(2),
                                                         phi_surf_grad(1) - phi_surf_grad(2)*q_normal(1)/q_normal(2),
                                                         0.0);
             Point<3,fad_double> phi_grad = phi_surf_grad + q_normal*q_dphi_dn;

             //std::cout<<"q_point="<<q_point<<"   q_normal="<<q_normal<<"   q_dphi_dn="<<q_dphi_dn<<std::endl;
             //cout<<q<<" phi_grad("<<phi_grad<<")  phi_surf_grad("<<phi_surf_grad<<")"<<endl;
             Point<3,fad_double> eta_grad(-q_normal(0)/q_normal(2),-q_normal(1)/q_normal(2),0.0);
             Point<3,fad_double> q_nodes_vel(q_x_dot,q_y_dot,q_z_dot);
             fluid_vel[q] = Point<3,fad_double>(fad_double(Vinf(0)),fad_double(Vinf(1)),fad_double(Vinf(2))) + phi_grad;
             fad_double fluid_vel_norm = fluid_vel[q].norm();
             fad_double horiz_fluid_vel_norm = sqrt(pow(fluid_vel[q](0),2.0)+pow(fluid_vel[q](1),2.0));
             fad_double eta_grad_norm = eta_grad.norm();
             if (fluid_vel_norm < 1e-3)
                fluid_vel_norm = -8.0e+05*pow(fluid_vel_norm,3.0) + 1.7e+03*pow(fluid_vel_norm,2.0) + 0.0001;
             if (horiz_fluid_vel_norm < 1e-3)
                horiz_fluid_vel_norm = -8.0e+05*pow(horiz_fluid_vel_norm,3.0) + 1.7e+03*pow(horiz_fluid_vel_norm,2.0) + 0.0001;
             if (eta_grad_norm < 1e-3)
                eta_grad_norm = -8.0e+05*pow(eta_grad_norm,3.0) + 1.7e+03*pow(eta_grad_norm,2.0) + 0.0001;
             Point<3,fad_double> horiz_vel_unit_vect(fluid_vel[q](0)/horiz_fluid_vel_norm,fluid_vel[q](1)/horiz_fluid_vel_norm,0.0);
             //fad_double cell_diameter;
             //for (unsigned int i=0; i<dofs_per_cell; ++i)
             //    {
             //    cell_diameter += pow(fluid_vel[q]*(Point<3,fad_double>(coors[3*i],coors[3*i+1],coors[3*i+2])-center),2.0)/dofs_per_cell;
             //    }
             //cell_diameter = sqrt(cell_diameter)*2;
             
             fad_double breaking_wave_added_pressure = 0.0;
             if ( (cell->material_id() == comp_dom.free_sur_ID1 ||
                   cell->material_id() == comp_dom.free_sur_ID2 ||
                   cell->material_id() == comp_dom.free_sur_ID3  ))
                if (sqrt(q_eta*q_eta+eta_grad*eta_grad*pow(Fn,4.0))-0.1725*Fn*Fn > 0)
                //if (eta_grad.norm() > 3.0)
                   {
                   //breaking_wave_added_pressure = pow(eta_grad*horiz_vel_unit_vect-7.0,2.0)*rho*g*cell_diameter;
                   //cout<<"Breaking wave damping ON at "<<q_point(0).val()<<","<<q_point(1).val()<<","<<q_point(2).val();
                   //cout<<"1)  ("<<eta_grad(0).val()<<","<<eta_grad(1).val()<<","<<eta_grad(2).val()<<") * ("<<fluid_vel[q](0).val()<<","<<fluid_vel[q](1).val()<<","<<fluid_vel[q](2).val()<<")"<<endl;
                   //cout<<"2) "<<q_eta.val()<<" "<<" "<<pow(Fn,4.0)<<" "<<0.1725*Fn*Fn<<endl;
                   breaking_wave_added_pressure = pow(sqrt(q_eta*q_eta+eta_grad*eta_grad*pow(Fn,4.0))-0.1725*Fn*Fn,2.0);
                   //if (breaking_wave_added_pressure < 0)
                   //   {
                   //   cout<<"1: "<<breaking_wave_added_pressure.val()<<endl;
                   //   }
                   fad_double factor;
                   if (q_eta+0.117*Fn*Fn < 0)
                      factor = 0;
                   else
                     if (q_eta+0.117*Fn*Fn < 0.01)
                        factor = pow(-1.0000e+04*pow(q_eta+0.117*Fn*Fn,3.0)+2.0000e+02*pow(q_eta+0.117*Fn*Fn,2.0),2.0);
                     else
                        factor = 2.0*(q_eta+0.117*Fn*Fn);
                   if (eta_grad*fluid_vel[q] > 0 )
                      breaking_wave_added_pressure *= eta_grad*eta_grad*factor;
                   else
                      breaking_wave_added_pressure *= eta_grad*eta_grad*factor*0.2; 
                   //if (breaking_wave_added_pressure < 0)
                   //   {
                   //   cout<<"2: "<<breaking_wave_added_pressure.val()<<"  factor: "<<factor<<endl;
                   //   }
                   //cout<<"3) "<<breaking_wave_added_pressure.val()<<" "<<endl;
                   breaking_wave_added_pressure *= rho*g/ref_height/12.0*(fluid_vel[q]*eta_grad)/fluid_vel_norm/eta_grad_norm;
                   //if (breaking_wave_added_pressure < 0)
                   //   {
                   //   cout<<"3: "<<breaking_wave_added_pressure.val()<<endl;
                   //   }
                   //fad_double test = fluid_vel[q]*eta_grad;
                   //fad_double test2 = 1/fluid_vel_norm/eta_grad_norm;
                   //cout<<"4) "<<breaking_wave_added_pressure.val()<<" "<<test.val()<<" "<<test2.val()<<" "<<rho*g/ref_height*2.0<<endl;
                   //cout<<"4) "<<breaking_wave_added_pressure.val()<<" "<<rho<<" "<<g<<" "<<ref_height<<endl;
                   //cout<<"Breaking wave damping ON at "<<q_point(0).val()<<","<<q_point(1).val()<<","<<q_point(2).val()<<"  --->  "<<breaking_wave_added_pressure.val()<<endl;
                   //cout<<"0) "<<breaking_wave_added_pressure.val()<<" "<<endl;
                   }
/*
               // this if is needed to compute the pressure patch behind the transom stern
               // pressure is estimated via Doctors regression formulas
               if (comp_dom.boat_model.is_transom)
                  if ( (cell->material_id() == comp_dom.free_sur_ID1 ||
                        cell->material_id() == comp_dom.free_sur_ID2 ||
                        cell->material_id() == comp_dom.free_sur_ID3  ))
                     if ( (q_point(1) < comp_dom.boat_model.PointRightTransom(1)) &&
                          (q_point(1) > comp_dom.boat_model.PointLeftTransom(1)) )
                        {
                        if (q_point(1) < 0)
                           curve = comp_dom.boat_model.left_transom_bspline;
                        else
                           curve = comp_dom.boat_model.right_transom_bspline;
           				   //this is the horizontal plane
                        Handle(Geom_Plane) zxPlane = new Geom_Plane(0.,1.,0.,q_point(1));
                        Handle(Geom_Curve) curve;
                        GeomAPI_IntCS Intersector(curve, horPlane);
                        gp_Pnt P;
                        gp_Vec V1;
                        int npoints = Intersector.NbPoints();
                        AssertThrow((npoints != 0), ExcMessage("Keel or transom curve is not intersecting with horizontal plane!"));

            double minDistance=1e7;
            double t,u,v;
            for (int j=0; j<npoints;++j)
                {
                Point<3> inters = Pnt(Intersector.Point(j+1));
                Intersector.Parameters(j+1,u,v,t);
                if (dP0.distance(inters) < minDistance)
                   {
                   minDistance = dP0.distance(inters);
                   dP = inters;
                   curve->D1(t,P,V1);
                   }
                }
                        }
*/
/*
                    if (q_point(0).val()<102 &&
                        q_point(0).val()>36.9 &&
                        abs(q_point(1).val())< 5.6 &&
                        abs(q_point(2).val())< 0.5 &&
                        fabs(t-0.005) < 1e-4 )             
                {
                std::cout<<"q_point=("<<q_point(0).val()<<","<<q_point(1).val()<<","<<q_point(2).val()<<")  q_dphi_dn="<<q_dphi_dn.val()<<endl;
                std::cout<<"phi_grad=("<<phi_grad(0).val()<<","<<phi_grad(1).val()<<","<<phi_grad(2).val()<<")"<<endl;
                std::cout<<"fluid_vel=("<<fluid_vel[q](0).val()<<","<<fluid_vel[q](1).val()<<","<<fluid_vel[q](2).val()<<")"<<endl;
                std::cout<<"q_nodes_vel=("<<q_nodes_vel(0).val()<<","<<q_nodes_vel(1).val()<<","<<q_nodes_vel(2).val()<<")"<<endl;
                std::cout<<"phi_surf_grad_corrected=("<<phi_surf_grad_corrected(0).val()<<","<<phi_surf_grad_corrected(1).val()<<","<<phi_surf_grad_corrected(2).val()<<")"<<endl;
                cout<<"Elevations:  "<<q_eta.val()<<"   VS   "<<q_init(2)<<endl;
                cout<<q<<" erhs("<<eta_dot_rhs_fun[q].val()<<")  prhs("<<phi_dot_rhs_fun[q].val()<<")"<<endl;
                //std::cout<<"local_DphiDt_rhs_2(i)="<<local_DphiDt_rhs_2(i)<<"local_DphiDt_rhs_2(i) ="<<local_DphiDt_rhs_2(i)<<")"<<endl;
                }
*/          
             transom_added_pressure = (1-eta_dry)*fad_double(g*q_init(2));
             //if (q_eta.val() > 1e-5)
             //cout<<"Later: "<<q_point(2)*gg(2)<<endl;
             q_JxW[q] = q_jac_det*ref_fe_v.JxW(q);

             //cout<<cell<<" "<<q<<" "<<phi_dot_rhs_fun[q]<<endl;
             //cout<<cell<<" "<<q<<" "<<q_JxW[q]<<endl;
             //cout.precision(0);
             //if ( (q_point(1).val() < comp_dom.boat_model.PointRightTransom(1)) &&
               //   (q_point(1).val() >= 0.0) &&
               //   (q_point(0).val() > comp_dom.boat_model.PointCenterTransom(0)-fabs(comp_dom.boat_model.PointCenterTransom(2)) ) &&
               //   (q_point(0).val() < comp_dom.boat_model.PointCenterTransom(0)+5*fabs(comp_dom.boat_model.PointCenterTransom(2)) )  )
              //  {
                //cout<<(int)cell->material_id()<<endl;
                //cout<<q<<"   "<<q_point(0).val()<<","<<q_point(1).val()<<","<<q_point(2).val()<<endl;
              //  cout<<transom_added_pressure.val()-g*q_eta.val()<<"    ("<<transom_added_pressure.val()<<" vs "<<g*q_eta.val()<<")"<<endl;
                //cout<<q<<" fvel("<<fluid_vel[q]<<")  fvel_norm="<<fluid_vel_norm<<"   q_JxW="<<q_JxW[q]<<endl;
                //cout<<q<<" erhs("<<eta_dot_rhs_fun[q]<<")  prhs("<<phi_dot_rhs_fun[q]<<")"<<endl;
                //cout<<q<<" phi_grad("<<phi_surf_grad_corrected(0).val()<<","<<phi_surf_grad_corrected(1).val()<<","<<phi_surf_grad_corrected(2).val()<<")"<<endl;//  phi_surf_grad("<<phi_surf_grad<<")"<<endl;
                //cout<<q<<"   "<<phi_dot_rhs_fun[q].val()<<endl;//" "<<phi_dot_rhs_fun[q].val()<<endl;
             //   }
             if (cell->material_id() == comp_dom.free_sur_ID1 ||
                 cell->material_id() == comp_dom.free_sur_ID2 ||
                 cell->material_id() == comp_dom.free_sur_ID3 )
                { //cout<<q<<"   "<<phi_dot_rhs_fun[q].val()<<endl;
                fad_double wave_damping_pressure = 0.0;
                //if (comp_dom.no_boat && q_point(0).val() > 0.0)
                fad_double Lx_boat = comp_dom.boat_model.boatWetLength;
                if (q_point(0).val() > Lx_boat*2.0)
                   wave_damping_pressure = -fad_double(1.0)*pow(q_point(0).val()-Lx_boat*2.0,2.0)/pow(Lx_boat*4.0,2.0)*q_dphi_dn;

                eta_dot_rhs_fun[q] = phi_grad*Point<3,fad_double>(fad_double(0.0),fad_double(0.0),fad_double(1.0)) +
                                     eta_grad*(q_nodes_vel-fluid_vel[q]);
                phi_dot_rhs_fun[q] = phi_grad*phi_grad/2 - q_point*gg + phi_surf_grad_corrected*(q_nodes_vel-fluid_vel[q])-
                                     breaking_wave_added_pressure +
                                     + (1-eta_dry)*fad_double(g*q_init(2)) + wave_damping_pressure;
                //if ( (q_point(0).val() < 3.10) && (q_point(0).val() > 3.03) &&
                //     (q_point(1).val() < 0.12) && (q_point(1).val() > 0) &&
                //     (q_point(2).val() > -1.0) )
                //     cout<<phi_dot_rhs_fun[q].val()<<endl;
                for (unsigned int i=0;i<dofs_per_cell;++i)
                    {
                    Point<3,fad_double> N_i_surf_grad(d11*ref_fe_v.shape_grad(i,q)[0]+d12*ref_fe_v.shape_grad(i,q)[1],
                                                     d21*ref_fe_v.shape_grad(i,q)[0]+d22*ref_fe_v.shape_grad(i,q)[1],
                                                     d31*ref_fe_v.shape_grad(i,q)[0]+d32*ref_fe_v.shape_grad(i,q)[1]);
                    fad_double N_i_supg = fad_double(ref_fe_v.shape_value(i,q)) +
                                          N_i_surf_grad*fluid_vel[q]/fluid_vel_norm*cell_diameter/sqrt(2);
                    loc_eta_res[i] -= eta_dot_rhs_fun[q]*N_i_supg*q_JxW[q];
                    loc_phi_res[i] -= phi_dot_rhs_fun[q]*N_i_supg*q_JxW[q];
                    loc_x_smooth_res[i] -= blend_factor*0*fad_double(ref_fe_v.shape_value(i,q))*q_JxW[q];
                    loc_y_smooth_res[i] -= blend_factor*0*fad_double(ref_fe_v.shape_value(i,q))*q_JxW[q];
                    local_DphiDt_rhs(i) += (phi_dot_rhs_fun[q]*N_i_supg*q_JxW[q]).val();
                    local_DphiDt_rhs_2(i) += (eta_dot_rhs_fun[q]*N_i_supg*q_JxW[q]).val();
                    //local_DphiDt_rhs_4(i) += (breaking_wave_added_pressure*fad_double(ref_fe_v.shape_value(i,q))*q_JxW[q]).val();
                    //local_DphiDt_rhs_4(i) += (transom_added_pressure*fad_double(ref_fe_v.shape_value(i,q))*q_JxW[q]).val();
                    local_DphiDt_rhs_4(i) += (wave_damping_pressure*fad_double(ref_fe_v.shape_value(i,q))*q_JxW[q]).val();
                    //cout<<q<<"  "<<i<<"   "<<phi_grad(2)<<"    "<<eta_grad<<"    "<<q_nodes_vel-fluid_vel[q]<<endl;
                    //cout<<q<<"  "<<i<<"   "<<N_i_surf_grad<<"    "<<fluid_vel[q]/fluid_vel_norm<<"   "<<cell_diameter/sqrt(2)<<endl;
                    //cout<<q<<"  "<<i<<"   "<<N_i_supg.val()<<"   "<<phi_dot_rhs_fun[q].val()<<"   "<<q_JxW[q].val()<<endl;
                    
                    //cout<<q<<"  "<<i<<" "<<q_JxW[q]<<endl;
                    for (unsigned int j=0;j<dofs_per_cell;++j)
                        {
                        //loc_eta_res[i] += fad_double(ref_fe_v.shape_value(j,q))*coors_dot[3*j+2]*N_i_supg*q_JxW[q];
                        //loc_phi_res[i] += fad_double(ref_fe_v.shape_value(j,q))*phis_dot[j]*N_i_supg*q_JxW[q];
                        //local_DphiDt_matrix.add(i,j,ref_fe_v.shape_value(j,q)*(N_i_supg*q_JxW[q]).val());
                        loc_supg_mass_matrix[i][j] += fad_double(ref_fe_v.shape_value(j,q))*N_i_supg*q_JxW[q];
                        Point<3,fad_double> N_j_surf_grad(d11*ref_fe_v.shape_grad(j,q)[0]+d12*ref_fe_v.shape_grad(j,q)[1],
                                                          d21*ref_fe_v.shape_grad(j,q)[0]+d22*ref_fe_v.shape_grad(j,q)[1],
                                                          d31*ref_fe_v.shape_grad(j,q)[0]+d32*ref_fe_v.shape_grad(j,q)[1]);
                        loc_stiffness_matrix[i][j] += N_i_surf_grad*N_j_surf_grad*q_JxW[q];
                        }
                    //if (fmax(abs(loc_eta_res[i].val()),abs(loc_phi_res[i].val()))>1e-6)
                    //   cout<<q<<"  "<<i<<"   "<<loc_eta_res[i].val()<<"("<<coors_dot[3*i+2].val()<<")  "<<loc_phi_res[i].val()<<"("<<phis_dot[i].val()<<")  "<<endl;   
                    }
                
                }

             if (cell->material_id() == comp_dom.pressure_sur_ID )
                { //cout<<q<<"   "<<phi_dot_rhs_fun[q].val()<<endl;
                //fad_double pressure = -q_point*gg*fad_double(rho);
                phi_dot_rhs_fun[q] = -q_ex_press/fad_double(rho) + phi_grad*phi_grad/2 - q_point*gg +
                                     phi_grad*(q_nodes_vel-fluid_vel[q]);
                for (unsigned int i=0;i<dofs_per_cell;++i)
                    {
                    Point<3,fad_double> N_i_surf_grad(d11*ref_fe_v.shape_grad(i,q)[0]+d12*ref_fe_v.shape_grad(i,q)[1],
                                                     d21*ref_fe_v.shape_grad(i,q)[0]+d22*ref_fe_v.shape_grad(i,q)[1],
                                                     d31*ref_fe_v.shape_grad(i,q)[0]+d32*ref_fe_v.shape_grad(i,q)[1]);
                    fad_double N_i_supg = fad_double(ref_fe_v.shape_value(i,q)) +
                                          N_i_surf_grad*fluid_vel[q]/fluid_vel_norm*cell_diameter/sqrt(2);
                    loc_phi_res[i] -= phi_dot_rhs_fun[q]*N_i_supg*q_JxW[q];
                    loc_x_smooth_res[i] -= blend_factor*0*fad_double(ref_fe_v.shape_value(i,q))*q_JxW[q];
                    loc_y_smooth_res[i] -= blend_factor*0*fad_double(ref_fe_v.shape_value(i,q))*q_JxW[q];
                    local_DphiDt_rhs(i) += (phi_dot_rhs_fun[q]*N_i_supg*q_JxW[q]).val();
                    //local_DphiDt_rhs_4(i) += (breaking_wave_added_pressure*fad_double(ref_fe_v.shape_value(i,q))*q_JxW[q]).val();
                    //local_DphiDt_rhs_4(i) += (transom_added_pressure*fad_double(ref_fe_v.shape_value(i,q))*q_JxW[q]).val();
                    local_DphiDt_rhs_4(i) += (0*fad_double(ref_fe_v.shape_value(i,q))*q_JxW[q]).val();
                    //cout<<q<<"  "<<i<<"   ("<<phi_surf_grad_corrected(0).val()<<","<<phi_surf_grad_corrected(1).val()<<","<<phi_surf_grad_corrected(2).val()<<")    ("<<q_nodes_vel(0).val()-fluid_vel[q](0).val()<<","<<q_nodes_vel(1).val()-fluid_vel[q](1).val()<<","<<q_nodes_vel(2).val()-fluid_vel[q](2).val()<<")"<<endl;
                    //cout<<q<<"  "<<i<<"   "<<N_i_surf_grad<<"    "<<fluid_vel[q]/fluid_vel_norm<<"   "<<cell_diameter/sqrt(2)<<endl;
                    //cout<<q<<"  "<<i<<"   "<<N_i_supg.val()<<"   "<<phi_dot_rhs_fun[q].val()<<"   "<<q_JxW[q].val()<<"   "<<loc_phi_res[i].val()<<endl;
                    
                    //cout<<q<<"  "<<i<<" "<<q_JxW[q]<<endl;
                    for (unsigned int j=0;j<dofs_per_cell;++j)
                        {
                        //loc_eta_res[i] += fad_double(ref_fe_v.shape_value(j,q))*coors_dot[3*j+2]*N_i_supg*q_JxW[q];
                        //loc_phi_res[i] += fad_double(ref_fe_v.shape_value(j,q))*phis_dot[j]*N_i_supg*q_JxW[q];
                        //local_DphiDt_matrix.add(i,j,ref_fe_v.shape_value(j,q)*(N_i_supg*q_JxW[q]).val());
                        loc_supg_mass_matrix[i][j] += fad_double(ref_fe_v.shape_value(j,q))*N_i_supg*q_JxW[q];
                        Point<3,fad_double> N_j_surf_grad(d11*ref_fe_v.shape_grad(j,q)[0]+d12*ref_fe_v.shape_grad(j,q)[1],
                                                          d21*ref_fe_v.shape_grad(j,q)[0]+d22*ref_fe_v.shape_grad(j,q)[1],
                                                          d31*ref_fe_v.shape_grad(j,q)[0]+d32*ref_fe_v.shape_grad(j,q)[1]);
                        loc_stiffness_matrix[i][j] += N_i_surf_grad*N_j_surf_grad*q_JxW[q];
                        }
                    //if (fmax(abs(loc_eta_res[i].val()),abs(loc_phi_res[i].val()))>1e-6)
                    //   cout<<q<<"  "<<i<<"   "<<loc_eta_res[i].val()<<"("<<coors_dot[3*i+2].val()<<")  "<<loc_phi_res[i].val()<<"("<<phis_dot[i].val()<<")  "<<endl;   
                    }
                
                }

             if (cell->material_id() == comp_dom.wall_sur_ID1 ||
                 cell->material_id() == comp_dom.wall_sur_ID2 ||
                 cell->material_id() == comp_dom.wall_sur_ID3 )
                {
                Point<3,fad_double> VV_inf(fad_double(Vinf(0)),fad_double(Vinf(1)),fad_double(Vinf(2)));
                fad_double dphi_dt = q_phi_dot - q_nodes_vel*phi_grad;
                //if (q==10)
                //   {
                //   cout<<"VEL: "<<q_nodes_vel(0).val()<<" "<<q_nodes_vel(1).val()<<" "<<q_nodes_vel(2).val()<<endl;
                //   cout<<dphi_dt.val()<<" "<<q_phi_dot.val()<<" "<<(q_nodes_vel*phi_grad).val()<<endl;
                //   }
                fad_double local_pressure = -(fad_double(rho)*(dphi_dt+q_point*gg+phi_grad*(phi_grad/2.0+VV_inf))*q_JxW[q]);
                loc_pressure_force += local_pressure*q_normal;

                Point<3,fad_double> q_mom((q_point(1)-baricenter_pos(1))*q_normal(2)-(q_point(2)-baricenter_pos(2))*q_normal(1),
                                          (q_point(2)-baricenter_pos(2))*q_normal(0)-(q_point(0)-baricenter_pos(0))*q_normal(2),
                                          (q_point(0)-baricenter_pos(0))*q_normal(1)-(q_point(1)-baricenter_pos(1))*q_normal(0));
                q_mom = q_mom*local_pressure;
                loc_pressure_moment+= q_mom;
                for (unsigned int i=0;i<dofs_per_cell;++i)
                    {
                    loc_dphi_dn_res[i] -= (q_normal*(q_nodes_vel-VV_inf))*
                                          fad_double(ref_fe_v.shape_value(i,q))*q_JxW[q];
                    //cout<<q<<"  "<<i<<"   "<<-(q_normal*Point<3,fad_double>(fad_double(Vinf(0)),fad_double(Vinf(1)),fad_double(Vinf(2)))).val()<<"    "<<cell->center()<<"    "<<endl;
                    //cout<<q<<"  "<<i<<"   "<<N_i_surf_grad<<"    "<<fluid_vel[q]/fluid_vel_norm<<"   "<<cell_diameter/sqrt(2)<<endl;
                    //cout<<q<<"  "<<i<<"   "<<N_i_supg.val()<<"   "<<phi_dot_res_fun[q].val()<<"   "<<q_JxW[q].val()<<endl;
                    local_DphiDt_rhs_3(i) += (-(q_normal*(q_nodes_vel-VV_inf))*
                                           fad_double(ref_fe_v.shape_value(i,q))*q_JxW[q]).val();
                    //cout<<"**** "<<loc_dphi_dn_res[i].val()<<" "<<local_DphiDt_rhs_3(i)<<" "<<loc_dphi_dn_res[i].val()+local_DphiDt_rhs_3(i)<<endl;
                    for (unsigned int j=0;j<dofs_per_cell;++j)
                        {
                        //loc_dphi_dn_res[i] += fad_double(ref_fe_v.shape_value(i,q))*fad_double(ref_fe_v.shape_value(j,q))*dphi_dns[j]*q_JxW[q];
                        loc_mass_matrix[i][j] += fad_double(ref_fe_v.shape_value(j,q))*fad_double(ref_fe_v.shape_value(i,q))*q_JxW[q]; 
                        }
                    //if (abs(loc_dphi_dn_res[i].val())>1e-7)
                    //   cout<<q<<"  "<<i<<"   "<<loc_dphi_dn_res[i].val()<<endl;   
                    }
                }

             if (cell->material_id() == comp_dom.inflow_sur_ID1 ||
                 cell->material_id() == comp_dom.inflow_sur_ID2  )
                {
                for (unsigned int i=0;i<dofs_per_cell;++i)
                    {
                    //fad_double k=0.62994; fad_double omega=2.4835; fad_double h=5.5;  fad_double a=0.1;
                    //fad_double inflow_norm_pot_grad = omega*a*cosh(k*(q_point(2)+h))/sinh(k*h)*cos(k*q_point(0)-omega*t);
                    //fad_double inflow_norm_pot_grad = k*a*cos(k*q_point(0)-omega*t);
                    //if (t<40.0)
                    //   inflow_norm_pot_grad *= 0.5*sin(3.141592654*(t)/40.0-3.141592654/2)+0.5;
                    //cout<<q<<"  "<<i<<"   "<<inflow_norm_pot_grad.val()<<"  x: "<<q_point(0).val()<<"  z: "<<q_point(2).val()<<endl; 
                    //loc_dphi_dn_res[i] -= inflow_norm_pot_grad*fad_double(ref_fe_v.shape_value(i,q))*q_JxW[q];
                    loc_dphi_dn_res[i] -= -(inflow_norm_potential_grad.value(Point<3>(q_point(0).val(),q_point(1).val(),q_point(2).val())))*
                                           fad_double(ref_fe_v.shape_value(i,q))*q_JxW[q];
                    //cout<<q<<"  "<<i<<"   "<<-(inflow_norm_potential_grad.value(Point<3>(q_point(0).val(),q_point(1).val(),q_point(2).val())))<<"    "<<cell->center()<<"    "<<endl;
                    //cout<<q<<"  "<<i<<"   "<<N_i_surf_grad<<"    "<<fluid_vel[q]/fluid_vel_norm<<"   "<<cell_diameter/sqrt(2)<<endl;
                    //cout<<q<<"  "<<i<<"   "<<N_i_supg.val()<<"   "<<phi_dot_res_fun[q].val()<<"   "<<q_JxW[q].val()<<endl;
                    local_DphiDt_rhs_3(i) += inflow_norm_potential_grad.value(Point<3>(q_point(0).val(),q_point(1).val(),q_point(2).val()))*
                                           (fad_double(ref_fe_v.shape_value(i,q))*q_JxW[q]).val();
                    //cout<<"**** "<<loc_dphi_dn_res[i].val()<<" "<<local_DphiDt_rhs_3(i)<<" "<<loc_dphi_dn_res[i].val()+local_DphiDt_rhs_3(i)<<endl;
                    for (unsigned int j=0;j<dofs_per_cell;++j)
                        {
                        //loc_dphi_dn_res[i] += fad_double(ref_fe_v.shape_value(i,q))*fad_double(ref_fe_v.shape_value(j,q))*dphi_dns[j]*q_JxW[q];
                        loc_mass_matrix[i][j] += fad_double(ref_fe_v.shape_value(j,q))*fad_double(ref_fe_v.shape_value(i,q))*q_JxW[q]; 
                        }
                    //if (abs(loc_dphi_dn_res[i].val())>1e-7)
                    //   cout<<q<<"  "<<i<<"   "<<loc_dphi_dn_res[i].val()<<endl;   
                    }
                }

             if (cell->material_id() != comp_dom.wall_sur_ID1 &&
                 cell->material_id() != comp_dom.wall_sur_ID2 &&
                 cell->material_id() != comp_dom.wall_sur_ID3 &&
                 cell->material_id() != comp_dom.free_sur_ID1 &&
                 cell->material_id() != comp_dom.free_sur_ID2 &&
                 cell->material_id() != comp_dom.free_sur_ID3 &&
                 cell->material_id() != comp_dom.inflow_sur_ID1 &&
                 cell->material_id() != comp_dom.inflow_sur_ID2 &&
                 cell->material_id() != comp_dom.pressure_sur_ID)
                {
                for (unsigned int i=0;i<dofs_per_cell;++i)
                    {
                    loc_dphi_dn_res[i] -= 0;
                    //cout<<q<<"  "<<i<<"   "<<-(q_normal*Point<3,fad_double>(fad_double(Vinf(0)),fad_double(Vinf(1)),fad_double(Vinf(2)))).val()<<"    "<<cell->center()<<"    "<<endl;
                    //cout<<q<<"  "<<i<<"   "<<N_i_surf_grad<<"    "<<fluid_vel[q]/fluid_vel_norm<<"   "<<cell_diameter/sqrt(2)<<endl;
                    //cout<<q<<"  "<<i<<"   "<<N_i_supg.val()<<"   "<<phi_dot_res_fun[q].val()<<"   "<<q_JxW[q].val()<<endl;
                    local_DphiDt_rhs_3(i) += 0;
                    //cout<<"**** "<<loc_dphi_dn_res[i].val()<<" "<<local_DphiDt_rhs_3(i)<<" "<<loc_dphi_dn_res[i].val()+local_DphiDt_rhs_3(i)<<endl;
                    for (unsigned int j=0;j<dofs_per_cell;++j)
                        {
                        //loc_dphi_dn_res[i] += fad_double(ref_fe_v.shape_value(i,q))*fad_double(ref_fe_v.shape_value(j,q))*dphi_dns[j]*q_JxW[q];
                        loc_mass_matrix[i][j] += fad_double(ref_fe_v.shape_value(j,q))*fad_double(ref_fe_v.shape_value(i,q))*q_JxW[q]; 
                        }
                    //if (abs(loc_dphi_dn_res[i].val())>1e-7)
                    //   cout<<q<<"  "<<i<<"   "<<loc_dphi_dn_res[i].val()<<endl;   
                    }
                }

             for (unsigned int i=0;i<dofs_per_cell;++i)
                 {
                 Point<3,fad_double> N_i_surf_grad(d11*ref_fe_v.shape_grad(i,q)[0]+d12*ref_fe_v.shape_grad(i,q)[1],
                                                   d21*ref_fe_v.shape_grad(i,q)[0]+d22*ref_fe_v.shape_grad(i,q)[1],
                                                   d31*ref_fe_v.shape_grad(i,q)[0]+d32*ref_fe_v.shape_grad(i,q)[1]);
                 fad_double N_i_supg = fad_double(ref_fe_v.shape_value(i,q)) +
                                       N_i_surf_grad*fluid_vel[q]/fluid_vel_norm*cell_diameter/sqrt(2);
                 for (unsigned int j=0;j<dofs_per_cell;++j)
                     {
                     local_DphiDt_matrix.add(i,j,ref_fe_v.shape_value(j,q)*(N_i_supg*q_JxW[q]).val());
                     local_DphiDt_matrix_2.add(i,j,ref_fe_v.shape_value(j,q)*ref_fe_v.shape_value(i,q)*(q_JxW[q]).val());
                     //loc_supg_mass_matrix[i][j] += fad_double(ref_fe_v.shape_value(j,q))*N_i_supg*q_JxW[q];
                     //loc_mass_matrix[i][j] += fad_double(ref_fe_v.shape_value(j,q))*fad_double(ref_fe_v.shape_value(i,q))*q_JxW[q]; 
                     }
                 } 
  	     }

         for (unsigned int i=0;i<dofs_per_cell;++i)
                 {
                 DphiDt_sys_rhs(local_dof_indices[i]) += local_DphiDt_rhs(i);
                 DphiDt_sys_rhs_2(local_dof_indices[i]) += local_DphiDt_rhs_2(i);
                 DphiDt_sys_rhs_3(local_dof_indices[i]) += local_DphiDt_rhs_3(i);
                 DphiDt_sys_rhs_4(local_dof_indices[i]) += local_DphiDt_rhs_4(i);
                 for (unsigned int j=0;j<dofs_per_cell;++j)
                     {
                     DphiDt_sys_matrix.add(local_dof_indices[i],local_dof_indices[j],local_DphiDt_matrix(i,j));
                     DphiDt_sys_matrix_2.add(local_dof_indices[i],local_dof_indices[j],local_DphiDt_matrix_2(i,j));
                     }
                 } 

         if (cell->material_id() == comp_dom.free_sur_ID1 ||
             cell->material_id() == comp_dom.free_sur_ID2 ||
             cell->material_id() == comp_dom.free_sur_ID3 )
             {
             for (unsigned int i=0;i<dofs_per_cell;++i)
                 {
                 //if ( (cell->center()(0) < 3.20) && (cell->center()(0) > 3.03) &&
                 //    (cell->center()(1) < 0.12) && (cell->center()(1) > 0) &&
                 //    (cell->center()(2) > -1.0) )
                 //    cout<<local_dof_indices[i]<<" *****   "<<phis_dot[i].val()<<" "<<phis[i].val()<<" "<<dphi_dns[i].val()<<endl;
                 for (unsigned int j=0;j<dofs_per_cell;++j)
                     {
                     //if ( (cell->center()(0) < 3.20) && (cell->center()(0) > 3.03) &&
                     //     (cell->center()(1) < 0.12) && (cell->center()(1) > 0) &&
                     //     (cell->center()(2) > -1.0) )
                     //     cout<<" o*****   "<<loc_supg_mass_matrix[i][j].val()<<endl;
                     loc_eta_res[i] += loc_supg_mass_matrix[i][j]*coors_dot[3*j+2];
                     loc_phi_res[i] += loc_supg_mass_matrix[i][j]*phis_dot[j];
                     //if (fabs(t<0.005)<1e-4) {cout<<cell<<" "<<i<<" "<<phis_dot[j]<<endl;}
                     loc_x_smooth_res[i] += loc_stiffness_matrix[i][j]*(x_displs[j]-(1-blend_factor)*comp_dom.old_map_points(3*local_dof_indices[j]));
                     loc_y_smooth_res[i] += loc_stiffness_matrix[i][j]*(y_displs[j]-(1-blend_factor)*comp_dom.old_map_points(3*local_dof_indices[j]+1));
                     }
                 if ( !constraints.is_constrained(local_dof_indices[i]) &&
                      !(comp_dom.flags[local_dof_indices[i]] & transom_on_water))
                    {
                    unsigned int ii = local_dof_indices[i];
                    eta_res[ii] += loc_eta_res[i].val();
                    phi_res[ii] += loc_phi_res[i].val();
                    //cout<<"* "<<cell<<" "<<local_dof_indices[i]<<" "<<loc_phi_res[i]<<endl;
                    for (unsigned int j=0;j<dofs_per_cell;++j)
                        {
                        unsigned int jj = local_dof_indices[j];
                        for (unsigned int k=0; k<dim; ++k)
                            jacobian_matrix.add(3*ii+2,
                                                3*jj+k,
                                                loc_eta_res[i].fastAccessDx(3*j+k));
                        jacobian_matrix.add(3*ii+2,
                                            jj+comp_dom.vector_dh.n_dofs(),
                                            loc_eta_res[i].fastAccessDx(3*dofs_per_cell+j));
                        jacobian_matrix.add(3*ii+2,
                                            jj+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            loc_eta_res[i].fastAccessDx(4*dofs_per_cell+j));
                       
                        for (unsigned int k=0; k<dim; ++k)
                            jacobian_matrix.add(ii+comp_dom.vector_dh.n_dofs(),
                                                3*jj+k,
                                                loc_phi_res[i].fastAccessDx(3*j+k));
                        jacobian_matrix.add(ii+comp_dom.vector_dh.n_dofs(),
                                            jj+comp_dom.vector_dh.n_dofs(),
                                            loc_phi_res[i].fastAccessDx(3*dofs_per_cell+j));
                        jacobian_matrix.add(ii+comp_dom.vector_dh.n_dofs(),
                                            jj+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            loc_phi_res[i].fastAccessDx(4*dofs_per_cell+j));
                                                      
                        for (unsigned int k=0; k<dim; ++k)
                            jacobian_dot_matrix.add(3*ii+2,
                                                    3*jj+k,
                                                    loc_eta_res[i].fastAccessDx(5*dofs_per_cell+3*j+k));
                        jacobian_dot_matrix.add(3*ii+2,
                                                jj+comp_dom.vector_dh.n_dofs(),
                                                loc_eta_res[i].fastAccessDx(8*dofs_per_cell+j));
                        jacobian_dot_matrix.add(3*ii+2,
                                                jj+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs(),
                                                loc_eta_res[i].fastAccessDx(9*dofs_per_cell+j));
                        for (unsigned int k=0; k<dim; ++k)
                            jacobian_dot_matrix.add(ii+comp_dom.vector_dh.n_dofs(),
                                                    3*jj+k,
                                                    loc_phi_res[i].fastAccessDx(5*dofs_per_cell+3*j+k));
                        jacobian_dot_matrix.add(ii+comp_dom.vector_dh.n_dofs(),
                                                jj+comp_dom.vector_dh.n_dofs(),
                                                loc_phi_res[i].fastAccessDx(8*dofs_per_cell+j));
                        jacobian_dot_matrix.add(ii+comp_dom.vector_dh.n_dofs(),
                                                jj+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs(),
                                                loc_phi_res[i].fastAccessDx(9*dofs_per_cell+j));
                        }

                    for (unsigned int k=0; k<dim; ++k)
                        jacobian_matrix.add(3*ii+2,
                                            k+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            loc_eta_res[i].fastAccessDx(k+10*dofs_per_cell));
                    for (unsigned int k=0; k<dim; ++k)
                        jacobian_matrix.add(3*ii+2,
                                            k+3+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            loc_eta_res[i].fastAccessDx(k+3+10*dofs_per_cell));
                    for (unsigned int k=0; k<dim; ++k)
                        jacobian_matrix.add(3*ii+2,
                                            k+6+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            loc_eta_res[i].fastAccessDx(k+6+10*dofs_per_cell));
                    for (unsigned int k=0; k<4; ++k)
                        jacobian_matrix.add(3*ii+2,
                                            k+9+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            loc_eta_res[i].fastAccessDx(k+9+10*dofs_per_cell));
                    for (unsigned int k=0; k<dim; ++k)
                        jacobian_matrix.add(ii+comp_dom.vector_dh.n_dofs(),
                                            k+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            loc_phi_res[i].fastAccessDx(k+10*dofs_per_cell));
                    for (unsigned int k=0; k<dim; ++k)
                        jacobian_matrix.add(ii+comp_dom.vector_dh.n_dofs(),
                                            k+3+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            loc_phi_res[i].fastAccessDx(k+3+10*dofs_per_cell));
                    for (unsigned int k=0; k<dim; ++k)
                        jacobian_matrix.add(ii+comp_dom.vector_dh.n_dofs(),
                                            k+6+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            loc_phi_res[i].fastAccessDx(k+6+10*dofs_per_cell));
                    for (unsigned int k=0; k<4; ++k)
                        jacobian_matrix.add(ii+comp_dom.vector_dh.n_dofs(),
                                            k+9+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            loc_phi_res[i].fastAccessDx(k+9+10*dofs_per_cell));

                    }
                 if ( !(constraints.is_constrained(local_dof_indices[i])) &&
                      !(comp_dom.flags[local_dof_indices[i]] & edge) )
                    {
                    unsigned int ii = local_dof_indices[i];
                    //if (local_dof_indices[i] == 601)
                    //   cout<<"££MATCH££ "<<local_dof_indices[i]<<": "<<loc_x_smooth_res[i]<<endl;
                    x_smoothing_res[ii] += loc_x_smooth_res[i].val();
                    y_smoothing_res[ii] += loc_y_smooth_res[i].val();
                    for (unsigned int j=0;j<dofs_per_cell;++j)
                        {
                        //if (local_dof_indices[i] == 601)
                        //   {
                        //   cout<<"("<<local_dof_indices[j]<<")"<<endl;
                        //   for (unsigned int k=0; k<dim; ++k)
                        //       cout<<"(-"<<coors[3*j+k]<<"-)"<<endl;
                        //   }
                        unsigned int jj = local_dof_indices[j];
                        for (unsigned int k=0; k<dim; ++k)
                            jacobian_matrix.add(3*ii,
                                                3*jj+k,
                                                loc_x_smooth_res[i].fastAccessDx(3*j+k));
                        jacobian_matrix.add(3*ii,
                                            jj+comp_dom.vector_dh.n_dofs(),
                                            loc_x_smooth_res[i].fastAccessDx(3*dofs_per_cell+j));
                        jacobian_matrix.add(3*ii,
                                            jj+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            loc_x_smooth_res[i].fastAccessDx(4*dofs_per_cell+j));


                        for (unsigned int k=0; k<dim; ++k)
                            jacobian_matrix.add(3*ii+1,
                                                3*jj+k,
                                                loc_y_smooth_res[i].fastAccessDx(3*j+k));
                        jacobian_matrix.add(3*ii+1,
                                            jj+comp_dom.vector_dh.n_dofs(),
                                            loc_y_smooth_res[i].fastAccessDx(3*dofs_per_cell+j));
                        jacobian_matrix.add(3*ii+1,
                                            jj+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            loc_y_smooth_res[i].fastAccessDx(4*dofs_per_cell+j));
                        }
                     for (unsigned int k=0; k<dim; ++k)
                         jacobian_matrix.add(3*ii,
                                             k+3+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                             loc_x_smooth_res[i].fastAccessDx(k+3+10*dofs_per_cell));
                     for (unsigned int k=0; k<4; ++k)
                         jacobian_matrix.add(3*ii,
                                             k+9+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                             loc_x_smooth_res[i].fastAccessDx(k+9+10*dofs_per_cell));
                     for (unsigned int k=0; k<dim; ++k)
                         jacobian_matrix.add(3*ii+1,
                                             k+3+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                             loc_y_smooth_res[i].fastAccessDx(k+3+10*dofs_per_cell));
                     for (unsigned int k=0; k<dim; ++k)
                         jacobian_matrix.add(3*ii+1,
                                             k+9+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                             loc_y_smooth_res[i].fastAccessDx(k+9+10*dofs_per_cell));
                    }
                 }
             }
         if (cell->material_id() == comp_dom.pressure_sur_ID)
             {
             for (unsigned int i=0;i<dofs_per_cell;++i)
                 {
                 for (unsigned int j=0;j<dofs_per_cell;++j)
                     {
                     loc_phi_res[i] += loc_supg_mass_matrix[i][j]*phis_dot[j];
                     //if (fabs(t<0.005)<1e-4) {cout<<cell<<" "<<i<<" "<<phis_dot[j]<<endl;}
                     }

                 if ( !constraints.is_constrained(local_dof_indices[i]) && 
                      !(comp_dom.flags[local_dof_indices[i]] & near_water))
                    {
                    unsigned int ii = local_dof_indices[i];
                    phi_res[ii] += loc_phi_res[i].val();
                    //cout<<"* "<<cell<<" "<<local_dof_indices[i]<<" "<<loc_phi_res[i]<<endl;
                    for (unsigned int j=0;j<dofs_per_cell;++j)
                        {
                        unsigned int jj = local_dof_indices[j];
                        //if (local_dof_indices[i]==44)
                        //   {
                        //   cout<<"44!!!: "<<cell<<endl;
                        //   for (unsigned int k=0; k<dim; ++k)
                        //       cout<<3*jj+k<<" ";
                        //   cout<<jj+comp_dom.vector_dh.n_dofs()<<" ";
                        //   cout<<jj+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()<<" ";
                        //   cout<<endl;
                        //   } 
                        //for (unsigned int k=0; k<dim; ++k)
                        //cout<<"ooo "<<loc_phi_res[i].fastAccessDx(3*j+k)<<endl;
                        for (unsigned int k=0; k<dim; ++k)
                            jacobian_matrix.add(ii+comp_dom.vector_dh.n_dofs(),
                                                3*jj+k,
                                                loc_phi_res[i].fastAccessDx(3*j+k));
                        jacobian_matrix.add(ii+comp_dom.vector_dh.n_dofs(),
                                            jj+comp_dom.vector_dh.n_dofs(),
                                            loc_phi_res[i].fastAccessDx(3*dofs_per_cell+j));
                        jacobian_matrix.add(ii+comp_dom.vector_dh.n_dofs(),
                                            jj+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            loc_phi_res[i].fastAccessDx(4*dofs_per_cell+j));                                                       
                        for (unsigned int k=0; k<dim; ++k)
                            jacobian_dot_matrix.add(ii+comp_dom.vector_dh.n_dofs(),
                                                    3*jj+k,
                                                    loc_phi_res[i].fastAccessDx(5*dofs_per_cell+3*j+k));
                        jacobian_dot_matrix.add(ii+comp_dom.vector_dh.n_dofs(),
                                                jj+comp_dom.vector_dh.n_dofs(),
                                                loc_phi_res[i].fastAccessDx(8*dofs_per_cell+j));
                        jacobian_dot_matrix.add(ii+comp_dom.vector_dh.n_dofs(),
                                                jj+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs(),
                                                loc_phi_res[i].fastAccessDx(9*dofs_per_cell+j));
                        }                     
                    }

                 }
             }
         if (cell->material_id() != comp_dom.free_sur_ID1 &&
             cell->material_id() != comp_dom.free_sur_ID2 &&
             cell->material_id() != comp_dom.free_sur_ID3 )
             {
             
             if (cell->material_id() == comp_dom.wall_sur_ID1 ||
                 cell->material_id() == comp_dom.wall_sur_ID2 ||
                 cell->material_id() == comp_dom.wall_sur_ID3 )
                {
                for (unsigned int d=0;d<3;++d)
                    {
                    pressure_force(d)+=loc_pressure_force(d).val();
                    pressure_moment(d)+=loc_pressure_moment(d).val();
                    }
                for (unsigned int j=0;j<dofs_per_cell;++j)
                    {
                    unsigned int jj = local_dof_indices[j];
                    
                    if (is_hull_x_translation_imposed != true)
                        {
                        for (unsigned int d=0;d<3;++d)
                            jacobian_matrix.add(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                                3*jj+d,
                                                -loc_pressure_force(0).fastAccessDx(3*j+d));
                        jacobian_matrix.add(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            comp_dom.vector_dh.n_dofs()+jj,
                                            -loc_pressure_force(0).fastAccessDx(j+3*dofs_per_cell));
                        jacobian_matrix.add(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+jj,
                                            -loc_pressure_force(0).fastAccessDx(j+4*dofs_per_cell));
                        for (unsigned int d=0;d<3;++d)
                            jacobian_dot_matrix.add(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                                    3*jj+d,
                                                    -loc_pressure_force(0).fastAccessDx(3*j+d+5*dofs_per_cell));
                        jacobian_dot_matrix.add(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                                comp_dom.vector_dh.n_dofs()+jj,
                                                -loc_pressure_force(0).fastAccessDx(j+8*dofs_per_cell));
                        jacobian_dot_matrix.add(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                                comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+jj,
                                                -loc_pressure_force(0).fastAccessDx(j+9*dofs_per_cell));
                        }
                    if (is_hull_y_translation_imposed != true)
                        {
                        for (unsigned int d=0;d<3;++d)
                            jacobian_matrix.add(1+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                                3*jj+d,
                                                -loc_pressure_force(1).fastAccessDx(3*j+d));
                        jacobian_matrix.add(1+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            comp_dom.vector_dh.n_dofs()+jj,
                                            -loc_pressure_force(1).fastAccessDx(j+3*dofs_per_cell));
                        jacobian_matrix.add(1+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+jj,
                                            -loc_pressure_force(1).fastAccessDx(j+4*dofs_per_cell));
                        for (unsigned int d=0;d<3;++d)
                            jacobian_dot_matrix.add(1+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                                    3*jj+d,
                                                    -loc_pressure_force(1).fastAccessDx(3*j+d+5*dofs_per_cell));
                        jacobian_dot_matrix.add(1+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                                comp_dom.vector_dh.n_dofs()+jj,
                                                -loc_pressure_force(1).fastAccessDx(j+8*dofs_per_cell));
                        jacobian_dot_matrix.add(1+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                                comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+jj,
                                                -loc_pressure_force(1).fastAccessDx(j+9*dofs_per_cell));
                        }
                    if (is_hull_z_translation_imposed != true)
                        {
                        for (unsigned int d=0;d<3;++d)
                            jacobian_matrix.add(2+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                                3*jj+d,
                                                -loc_pressure_force(2).fastAccessDx(3*j+d));
                        jacobian_matrix.add(2+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            comp_dom.vector_dh.n_dofs()+jj,
                                            -loc_pressure_force(2).fastAccessDx(j+3*dofs_per_cell));
                        jacobian_matrix.add(2+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+jj,
                                            -loc_pressure_force(2).fastAccessDx(j+4*dofs_per_cell));
                        for (unsigned int d=0;d<3;++d)
                            jacobian_dot_matrix.add(2+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                                    3*jj+d,
                                                    -loc_pressure_force(2).fastAccessDx(3*j+d+5*dofs_per_cell));
                        jacobian_dot_matrix.add(2+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                                comp_dom.vector_dh.n_dofs()+jj,
                                                -loc_pressure_force(2).fastAccessDx(j+8*dofs_per_cell));
                        jacobian_dot_matrix.add(2+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                                comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+jj,
                                                -loc_pressure_force(2).fastAccessDx(j+9*dofs_per_cell));
                        }
                    
                    //if (is_hull_x_translation_imposed != true)
                    /*
                        {
                        for (unsigned int d=0;d<3;++d)
                            jacobian_matrix.add(6+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                                3*jj+d,
                                                -loc_pressure_moment(0).fastAccessDx(3*j+d));
                        jacobian_matrix.add(6+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            comp_dom.vector_dh.n_dofs()+jj,
                                            -loc_pressure_moment(0).fastAccessDx(j+3*dofs_per_cell));
                        jacobian_matrix.add(6+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+jj,
                                            -loc_pressure_moment(0).fastAccessDx(j+4*dofs_per_cell));
                        for (unsigned int d=0;d<3;++d)
                            jacobian_dot_matrix.add(6+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                                    3*jj+d,
                                                    -loc_pressure_moment(0).fastAccessDx(3*j+d+5*dofs_per_cell));
                        jacobian_dot_matrix.add(6+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                                comp_dom.vector_dh.n_dofs()+jj,
                                                -loc_pressure_moment(0).fastAccessDx(j+8*dofs_per_cell));
                        jacobian_dot_matrix.add(6+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                                comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+jj,
                                                -loc_pressure_moment(0).fastAccessDx(j+9*dofs_per_cell));
                        }
                    */
                    //if (is_hull_y_translation_imposed != true)
                        {
                        for (unsigned int d=0;d<3;++d)
                            jacobian_matrix.add(7+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                                3*jj+d,
                                                -loc_pressure_moment(1).fastAccessDx(3*j+d));
                        jacobian_matrix.add(7+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            comp_dom.vector_dh.n_dofs()+jj,
                                            -loc_pressure_moment(1).fastAccessDx(j+3*dofs_per_cell));
                        jacobian_matrix.add(7+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+jj,
                                            -loc_pressure_moment(1).fastAccessDx(j+4*dofs_per_cell));
                        for (unsigned int d=0;d<3;++d)
                            jacobian_dot_matrix.add(7+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                                    3*jj+d,
                                                    -loc_pressure_moment(1).fastAccessDx(3*j+d+5*dofs_per_cell));
                        jacobian_dot_matrix.add(7+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                                comp_dom.vector_dh.n_dofs()+jj,
                                                -loc_pressure_moment(1).fastAccessDx(j+8*dofs_per_cell));
                        jacobian_dot_matrix.add(7+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                                comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+jj,
                                                -loc_pressure_moment(1).fastAccessDx(j+9*dofs_per_cell));
                        }
                    /*
                    //if (is_hull_z_translation_imposed != true)
                        {
                        for (unsigned int d=0;d<3;++d)
                            jacobian_matrix.add(8+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                                3*jj+d,
                                                -loc_pressure_moment(2).fastAccessDx(3*j+d));
                        jacobian_matrix.add(8+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            comp_dom.vector_dh.n_dofs()+jj,
                                            -loc_pressure_moment(2).fastAccessDx(j+3*dofs_per_cell));
                        jacobian_matrix.add(8+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+jj,
                                            -loc_pressure_moment(2).fastAccessDx(j+4*dofs_per_cell));
                        for (unsigned int d=0;d<3;++d)
                            jacobian_dot_matrix.add(8+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                                    3*jj+d,
                                                    -loc_pressure_moment(2).fastAccessDx(3*j+d+5*dofs_per_cell));
                        jacobian_dot_matrix.add(8+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                                comp_dom.vector_dh.n_dofs()+jj,
                                                -loc_pressure_moment(2).fastAccessDx(j+8*dofs_per_cell));
                        jacobian_dot_matrix.add(8+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                                comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+jj,
                                                -loc_pressure_moment(2).fastAccessDx(j+9*dofs_per_cell));
                        }
                    */
                    }
                
                if (is_hull_x_translation_imposed != true)
                    {
                    
                    for (unsigned int d=0;d<3;++d)
                        jacobian_matrix.add(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            d+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            -loc_pressure_force(0).fastAccessDx(d+10*dofs_per_cell));
                    for (unsigned int d=0;d<3;++d)
                        jacobian_matrix.add(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            d+3+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            -loc_pressure_force(0).fastAccessDx(d+3+10*dofs_per_cell));                        
                    for (unsigned int d=0;d<3;++d)
                        jacobian_matrix.add(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            d+6+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            -loc_pressure_force(0).fastAccessDx(d+6+10*dofs_per_cell));
                    for (unsigned int d=0;d<4;++d)
                        jacobian_matrix.add(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            d+9+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            -loc_pressure_force(0).fastAccessDx(d+9+10*dofs_per_cell));
                    }
                if (is_hull_y_translation_imposed != true)
                    {
                    for (unsigned int d=0;d<3;++d)
                        jacobian_matrix.add(1+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            d+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            -loc_pressure_force(1).fastAccessDx(d+10*dofs_per_cell));
                    for (unsigned int d=0;d<3;++d)
                        jacobian_matrix.add(1+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            d+3+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            -loc_pressure_force(1).fastAccessDx(d+3+10*dofs_per_cell));
                    for (unsigned int d=0;d<3;++d)
                        jacobian_matrix.add(1+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            d+6+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            -loc_pressure_force(1).fastAccessDx(d+6+10*dofs_per_cell));
                    for (unsigned int d=0;d<4;++d)
                        jacobian_matrix.add(1+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            d+9+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            -loc_pressure_force(1).fastAccessDx(d+9+10*dofs_per_cell));                    
                    }
                if (is_hull_z_translation_imposed != true)
                    {
                    for (unsigned int d=0;d<3;++d)
                        jacobian_matrix.add(2+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            d+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            -loc_pressure_force(2).fastAccessDx(d+10*dofs_per_cell));
                    for (unsigned int d=0;d<3;++d)
                        jacobian_matrix.add(2+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            d+3+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            -loc_pressure_force(2).fastAccessDx(d+3+10*dofs_per_cell));
                    for (unsigned int d=0;d<3;++d)
                        jacobian_matrix.add(2+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            d+6+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            -loc_pressure_force(2).fastAccessDx(d+6+10*dofs_per_cell));
                    for (unsigned int d=0;d<4;++d)
                        jacobian_matrix.add(2+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            d+9+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            -loc_pressure_force(2).fastAccessDx(d+9+10*dofs_per_cell)); 
                    }
                
                //if (is_hull_y_translation_imposed != true)
                /*
                    {
                    for (unsigned int d=0;d<3;++d)
                        jacobian_matrix.add(6+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            d+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            -loc_pressure_moment(0).fastAccessDx(d+10*dofs_per_cell));
                    for (unsigned int d=0;d<3;++d)
                        jacobian_matrix.add(6+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            d+3+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            -loc_pressure_moment(0).fastAccessDx(d+3+10*dofs_per_cell));                        
                    for (unsigned int d=0;d<3;++d)
                        jacobian_matrix.add(6+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            d+6+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            -loc_pressure_moment(0).fastAccessDx(d+6+10*dofs_per_cell));
                    for (unsigned int d=0;d<4;++d)
                        jacobian_matrix.add(6+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            d+9+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            -loc_pressure_moment(0).fastAccessDx(d+9+10*dofs_per_cell));
                    }
                */
                //if (is_hull_y_translation_imposed != true)
                    {
                    for (unsigned int d=0;d<3;++d)
                        jacobian_matrix.add(7+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            d+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            -loc_pressure_moment(1).fastAccessDx(d+10*dofs_per_cell));
                    for (unsigned int d=0;d<3;++d)
                        jacobian_matrix.add(7+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            d+3+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            -loc_pressure_moment(1).fastAccessDx(d+3+10*dofs_per_cell));                        
                    for (unsigned int d=0;d<3;++d)
                        jacobian_matrix.add(7+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            d+6+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            -loc_pressure_moment(1).fastAccessDx(d+6+10*dofs_per_cell));
                    for (unsigned int d=0;d<4;++d)
                        jacobian_matrix.add(7+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            d+9+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            -loc_pressure_moment(1).fastAccessDx(d+9+10*dofs_per_cell));
                    }
                /*
                //if (is_hull_y_translation_imposed != true)
                    {
                    for (unsigned int d=0;d<3;++d)
                        jacobian_matrix.add(8+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            d+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            -loc_pressure_moment(2).fastAccessDx(d+10*dofs_per_cell));
                    for (unsigned int d=0;d<3;++d)
                        jacobian_matrix.add(8+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            d+3+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            -loc_pressure_moment(2).fastAccessDx(d+3+10*dofs_per_cell));                        
                    for (unsigned int d=0;d<3;++d)
                        jacobian_matrix.add(8+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            d+6+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            -loc_pressure_moment(2).fastAccessDx(d+6+10*dofs_per_cell));
                    for (unsigned int d=0;d<4;++d)
                        jacobian_matrix.add(8+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            d+9+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            -loc_pressure_moment(2).fastAccessDx(d+9+10*dofs_per_cell));
                    }
                */
                }
                
             for (unsigned int i=0;i<dofs_per_cell;++i)
                 {
                 unsigned int ii = local_dof_indices[i];
                 for (unsigned int j=0;j<dofs_per_cell;++j)
                     {
                     loc_dphi_dn_res[i] += loc_mass_matrix[i][j]*dphi_dns[j];
                     }
                 if (!constraints.is_constrained(ii))
                    {
                    dphi_dn_res[ii] += loc_dphi_dn_res[i].val();
                    
                    for (unsigned int j=0;j<dofs_per_cell;++j)
                        {
                        unsigned int jj = local_dof_indices[j];
                        for (unsigned int k=0; k<dim; ++k)
                            jacobian_matrix.add(ii+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs(),
                                                3*jj+k,
                                                loc_dphi_dn_res[i].fastAccessDx(3*j+k));
                        for (unsigned int k=0; k<dim; ++k)
                            jacobian_dot_matrix.add(ii+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs(),
                                                    3*jj+k,
                                                    loc_dphi_dn_res[i].fastAccessDx(5*dofs_per_cell+3*j+k));
                        jacobian_matrix.add(ii+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            jj+comp_dom.vector_dh.n_dofs(),
                                            loc_dphi_dn_res[i].fastAccessDx(3*dofs_per_cell+j));
                        jacobian_matrix.add(ii+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            jj+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            loc_dphi_dn_res[i].fastAccessDx(4*dofs_per_cell+j));
                        }
                    for (unsigned int k=0; k<dim; ++k)
                        jacobian_matrix.add(ii+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            k+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            loc_dphi_dn_res[i].fastAccessDx(k+10*dofs_per_cell));
                    for (unsigned int k=0; k<dim; ++k)
                        jacobian_matrix.add(ii+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            k+3+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            loc_dphi_dn_res[i].fastAccessDx(k+3+10*dofs_per_cell));
                    for (unsigned int k=0; k<dim; ++k)
                        jacobian_matrix.add(ii+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            k+6+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            loc_dphi_dn_res[i].fastAccessDx(k+6+10*dofs_per_cell));
                    for (unsigned int k=0; k<4; ++k)
                        jacobian_matrix.add(ii+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            k+9+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            loc_dphi_dn_res[i].fastAccessDx(k+9+10*dofs_per_cell));
                    }
                 }
             }
         }
      }

//let's build the residual (and jacobian) on the hull displacements and velocities (differential components) 
  Point<3> hull_lin_vel_res(0.0,0.0,0.0);
  Point<3> hull_lin_displ_res(0.0,0.0,0.0);

  Point<3> g_vec(0.0,0.0,-g);

  hull_lin_vel_res+= comp_dom.boat_model.boat_mass*hull_lin_vel_dot;
// this jacobian will be completed by the derivatives of the pressure force with respect to all the dofs 
  for (unsigned int d=0; d<3;++d)
      jacobian_dot_matrix.add(d+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                              d+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                              comp_dom.boat_model.boat_mass);


  hull_lin_displ_res+= hull_lin_displ_dot - hull_lin_vel;
  for (unsigned int d=0; d<3;++d)
      {
      jacobian_dot_matrix.add(d+3+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                              d+3+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                              1.0);    
      jacobian_matrix.add(d+3+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                          d+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                          -1.0);
      }

//once looping on the cells is over, the pressure force and moment are ready: let's
//add them to the residual for the boat acceleration
cout<<"The current pressure force is: "<<pressure_force<<endl;
cout<<"The current pressure moment is: "<<pressure_moment<<endl;
if (is_hull_x_translation_imposed)
   {
   Point<1> time(t);
   //double accel = hull_x_axis_translation.laplacian(time);
   Point<1> delta_t(1e-5);
   double accel = ( hull_x_axis_translation.value(time+delta_t) -
                    2*hull_x_axis_translation.value(time) +
                    hull_x_axis_translation.value(time+(-1.0*delta_t)) ) / delta_t(0) / delta_t(0);
   cout<<"********Test Laplacian x: "<<accel<<endl;
   hull_lin_vel_res(0)-= accel*comp_dom.boat_model.boat_mass;
   }
else
   {
   hull_lin_vel_res(0)-= pressure_force(0);
   }

if (is_hull_y_translation_imposed)
   {
   Point<1> time(t);
   //double accel = hull_y_axis_translation.laplacian(time);
   Point<1> delta_t(1e-5);
   double accel = ( hull_y_axis_translation.value(time+delta_t) -
                    2*hull_y_axis_translation.value(time) +
                    hull_y_axis_translation.value(time+(-1.0*delta_t)) ) / delta_t(0) / delta_t(0);
   cout<<"********Test Laplacian y: "<<accel<<endl;
   hull_lin_vel_res(1)-= accel*comp_dom.boat_model.boat_mass;
   }
else
   {
   hull_lin_vel_res(1)-= pressure_force(1);
   }

if (is_hull_z_translation_imposed)
   {
   Point<1> time(t);
   //double accel = hull_z_axis_translation.laplacian(time);
   Point<1> delta_t(1e-5);
   double accel = ( hull_z_axis_translation.value(time+delta_t) -
                    2*hull_z_axis_translation.value(time) +
                    hull_z_axis_translation.value(time+(-1.0*delta_t)) ) / delta_t(0) / delta_t(0);
   //double accel = -0.01*pow(2.0*dealii::numbers::PI,2.0)*sin(2.0*dealii::numbers::PI*t); 
   cout<<"********Test Laplacian z: "<<accel<<endl;
   hull_lin_vel_res(2)-= accel*comp_dom.boat_model.boat_mass;
   }
else
   {
   hull_lin_vel_res(2)-= pressure_force(2);
   hull_lin_vel_res(2)-= (-g*comp_dom.boat_model.boat_mass);
   cout<<"************ "<<pressure_force(2)<<" vs "<<-g*comp_dom.boat_model.boat_mass<<endl;
   }

  Point<3,fad_double> hull_quaternions_vect_res(0.0,0.0,0.0); 
  fad_double hull_quaternions_scal_res;
  fad_double omega_x,omega_y,omega_z,v_x,v_y,v_z,s;
  fad_double omega_x_dot,omega_y_dot,omega_z_dot,v_x_dot,v_y_dot,v_z_dot,s_dot;

  omega_x_dot = hull_ang_vel_dot(0);
  omega_y_dot = hull_ang_vel_dot(1);
  omega_z_dot = hull_ang_vel_dot(2);
  v_x_dot = hull_quat_vector_dot(0);
  v_y_dot = hull_quat_vector_dot(1);
  v_z_dot = hull_quat_vector_dot(2);
  s_dot = hull_quat_scalar_dot;
  omega_x = hull_ang_vel(0);
  omega_y = hull_ang_vel(1);
  omega_z = hull_ang_vel(2);
  v_x = hull_quat_vector(0);
  v_y = hull_quat_vector(1);
  v_z = hull_quat_vector(2);
  s = hull_quat_scalar;

  omega_x_dot.diff(0,14);
  omega_y_dot.diff(1,14);
  omega_z_dot.diff(2,14);
  v_x_dot.diff(3,14);
  v_y_dot.diff(4,14);
  v_z_dot.diff(5,14);
  s_dot.diff(6,14);
  omega_x.diff(7,14);
  omega_y.diff(8,14);
  omega_z.diff(9,14);
  v_x.diff(10,14);
  v_y.diff(11,14);
  v_z.diff(12,14);
  s.diff(13,14);

  Point<3,fad_double> RotMatRow1(1-2*v_y*v_y-2*v_z*v_z, 2*v_x*v_y-2*s*v_z, 2*v_x*v_z+2*s*v_y);
  Point<3,fad_double> RotMatRow2(2*v_x*v_y+2*s*v_z, 1-2*v_x*v_x-2*v_z*v_z, 2*v_y*v_z-2*s*v_x);
  Point<3,fad_double> RotMatRow3(2*v_x*v_z-2*s*v_y, 2*v_y*v_z+2*s*v_x, 1-2*v_y*v_y-2*v_x*v_x);
  Point<3,fad_double> RelInertiaMatRow1(comp_dom.Ixx,comp_dom.Ixy,comp_dom.Ixz);
  Point<3,fad_double> RelInertiaMatRow2(comp_dom.Ixy,comp_dom.Iyy,comp_dom.Iyz);
  Point<3,fad_double> RelInertiaMatRow3(comp_dom.Ixz,comp_dom.Iyz,comp_dom.Izz);

  Point<3,fad_double> servMatRow1(RotMatRow1*RelInertiaMatRow1, RotMatRow1*RelInertiaMatRow2, RotMatRow1*RelInertiaMatRow3);
  Point<3,fad_double> servMatRow2(RotMatRow2*RelInertiaMatRow1, RotMatRow2*RelInertiaMatRow2, RotMatRow2*RelInertiaMatRow3);
  Point<3,fad_double> servMatRow3(RotMatRow3*RelInertiaMatRow1, RotMatRow3*RelInertiaMatRow2, RotMatRow3*RelInertiaMatRow3);

  Point<3,fad_double> AbsInertiaMatRow1(servMatRow1*RotMatRow1, servMatRow1*RotMatRow2, servMatRow1*RotMatRow3);
  Point<3,fad_double> AbsInertiaMatRow2(servMatRow2*RotMatRow1, servMatRow2*RotMatRow2, servMatRow2*RotMatRow3);
  Point<3,fad_double> AbsInertiaMatRow3(servMatRow3*RotMatRow1, servMatRow3*RotMatRow2, servMatRow3*RotMatRow3);

  Point<3,fad_double> omega(omega_x,omega_y,omega_z);

  Point<3,fad_double> WMatRow1(     0.0, -omega_z,  omega_y);
  Point<3,fad_double> WMatRow2( omega_z,      0.0, -omega_x);
  Point<3,fad_double> WMatRow3(-omega_y,  omega_x,      0.0);

  Point<3,fad_double> sservMatRow1(WMatRow1*AbsInertiaMatRow1, WMatRow1*AbsInertiaMatRow2, WMatRow1*AbsInertiaMatRow3);
  Point<3,fad_double> sservMatRow2(WMatRow2*AbsInertiaMatRow1, WMatRow2*AbsInertiaMatRow2, WMatRow2*AbsInertiaMatRow3);
  Point<3,fad_double> sservMatRow3(WMatRow3*AbsInertiaMatRow1, WMatRow3*AbsInertiaMatRow2, WMatRow3*AbsInertiaMatRow3);    

  Point<3,fad_double> omega_dot(omega_x_dot,omega_y_dot,omega_z_dot);

  Point<3,fad_double> vv(v_x,v_y,v_z);
  vv/=sqrt(v_x*v_x+v_y*v_y+v_z*v_z+s*s);
  
  fad_double ss = s/sqrt(v_x*v_x+v_y*v_y+v_z*v_z+s*s);

  
  Point<3,fad_double> vv_dot(v_x_dot,v_y_dot,v_z_dot);

  Point<3,fad_double> rhs_quat_vect(ss*omega_x+WMatRow1*vv,ss*omega_y+WMatRow2*vv,ss*omega_z+WMatRow1*vv);
/*
  Point<3,fad_double> hull_ang_vel_res(omega_x_dot,omega_y_dot,omega_z_dot);

  for (unsigned int k=0; k<3; ++k)
      {
      jacobian_dot_matrix.add(6+k+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                              6+k+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                              1.0); 
      } 
*/

  Point<3,fad_double> hull_ang_vel_res(AbsInertiaMatRow1*omega_dot, AbsInertiaMatRow2*omega_dot, AbsInertiaMatRow3*omega_dot);
  hull_ang_vel_res = hull_ang_vel_res + Point<3,fad_double>(sservMatRow1*omega,sservMatRow3*omega,sservMatRow3*omega);
//cout<<"OOOOOOOOOOOOO   "<<hull_ang_vel_res<<endl;
  for (unsigned int k=0; k<3; ++k)
      {
      for (unsigned int d=0; d<7; ++d)
          jacobian_dot_matrix.add(6+k+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                  6+d+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                  hull_ang_vel_res(k).fastAccessDx(d));
      for (unsigned int d=0; d<7; ++d)
          jacobian_matrix.add(6+k+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                              6+d+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                              hull_ang_vel_res(k).fastAccessDx(7+d));
      }

//for (unsigned int k=0; k<3; ++k)
//      {
//      //hull_ang_vel_res(k) = hull_ang_vel_res(k) - 1.0*(t<10.0? t/10.0:1.0)*pressure_moment(k);
//      }
//hull_ang_vel_res(1) = hull_ang_vel_res(1) + 25.0*(t<4.0? t/4.0:1.0)*pow(numbers::PI/2.0,2.0)*cos(numbers::PI/2.0*t);
  cout<<"Moment Fraction Considered: "<<(t<10.0? 0.0:(t<15.0? (t-10.0)/5.0:1.0))<<endl;
  hull_ang_vel_res(1) = hull_ang_vel_res(1) - pressure_moment(1);

  hull_quaternions_scal_res = s_dot + omega*vv/2.0;

  for (unsigned int k=0; k<3; ++k)
      {
      hull_quaternions_vect_res(k) = vv_dot(k) - rhs_quat_vect(k)/2.0;  
      }

  for (unsigned int k=0; k<3; ++k)
      {
      for (unsigned int d=0; d<7; ++d)
          jacobian_dot_matrix.add(9+k+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                  6+d+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                  hull_quaternions_vect_res(k).fastAccessDx(d));
      for (unsigned int d=0; d<7; ++d)
          jacobian_matrix.add(9+k+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                              6+d+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                              hull_quaternions_vect_res(k).fastAccessDx(7+d));
      }
     for (unsigned int d=0; d<7; ++d)
          jacobian_dot_matrix.add(12+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                  6+d+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                  hull_quaternions_scal_res.fastAccessDx(d));
      for (unsigned int d=0; d<7; ++d)
          jacobian_matrix.add(12+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                              6+d+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                              hull_quaternions_scal_res.fastAccessDx(7+d));

  Vector<double> hull_rigid_modes_res(13);
  for (unsigned int k=0; k<3; ++k)
      {
      hull_rigid_modes_res(k) = hull_lin_vel_res(k);
      hull_rigid_modes_res(k+3) = hull_lin_displ_res(k);
      hull_rigid_modes_res(k+6) = hull_ang_vel_res(k).val();
      }
   for (unsigned int k=0; k<3; ++k)
      {
      hull_rigid_modes_res(k+9) = hull_quaternions_vect_res(k).val();
      }
   hull_rigid_modes_res(12) = hull_quaternions_scal_res.val();
//for (unsigned int i=0; i<comp_dom.dh.n_dofs();++i)
//    cout<<i<<"--->  "<<eta_res[i]<<endl;;

/*
typedef const unsigned int *SparsityPattern::iterator;
     if (fabs(t-0.231589)<0.00001)
{
for (unsigned int i=0; i<this->n_dofs(); ++i)
    for (SparsityPattern::iterator col=DphiDt_sparsity_pattern.begin(i); col!=DphiDt_sparsity_pattern.end(i); ++col)
        {
        unsigned int j = col->column();
        cout<<" "<<i<<" "<<j<<" "<<DphiDt_sys_matrix(i,j)<<" "<<DphiDt_sys_matrix_2(i,j)<<endl;
        }
}
*/



     SparseMatrix<double> DphiDt_sys_matrix_copy;
     DphiDt_sys_matrix_copy.reinit(DphiDt_sparsity_pattern);
     DphiDt_sys_matrix_copy.copy_from(DphiDt_sys_matrix);

     SparseMatrix<double> DphiDt_sys_matrix_copy_2;
     DphiDt_sys_matrix_copy_2.reinit(DphiDt_sparsity_pattern);
     DphiDt_sys_matrix_copy_2.copy_from(DphiDt_sys_matrix_2);     

     constraints.condense(DphiDt_sys_matrix,DphiDt_sys_rhs);
     constraints.condense(DphiDt_sys_matrix_copy,DphiDt_sys_rhs_2);
     constraints.condense(DphiDt_sys_matrix_2,DphiDt_sys_rhs_3);
     constraints.condense(DphiDt_sys_matrix_copy_2,DphiDt_sys_rhs_4);
 

     SparseDirectUMFPACK DphiDt_direct;
     SparseDirectUMFPACK DphiDt_direct_copy;
     SparseDirectUMFPACK DphiDt_direct_copy_2;
     SparseDirectUMFPACK DphiDt_direct_2;
     
     DphiDt_direct.initialize(DphiDt_sys_matrix);
     DphiDt_direct_copy.initialize(DphiDt_sys_matrix_copy);
     DphiDt_direct_copy_2.initialize(DphiDt_sys_matrix_copy_2);
     DphiDt_direct_2.initialize(DphiDt_sys_matrix_2);

     DphiDt_direct.vmult(DphiDt_sys_solution, DphiDt_sys_rhs); // solving for phi_dot
     constraints.distribute(DphiDt_sys_solution);
     //cout<<"---------"<<endl;
     //cout<<DphiDt_sys_solution<<endl;
     //cout<<"---------"<<endl;
     DphiDt_direct_copy.vmult(DphiDt_sys_solution_2, DphiDt_sys_rhs_2); // solving for eta_dot
     constraints.distribute(DphiDt_sys_solution_2);
     DphiDt_direct_2.vmult(DphiDt_sys_solution_3, DphiDt_sys_rhs_3); // solving for dphi_dn
     constraints.distribute(DphiDt_sys_solution_3);
     DphiDt_direct_copy_2.vmult(break_wave_press, DphiDt_sys_rhs_4); // solving for eta_dot
     constraints.distribute(break_wave_press);     

     Vector<double> RES(DphiDt_sys_solution.size());
     DphiDt_sys_matrix.vmult(RES,DphiDt_sys_solution_2);
     RES*=-1.0;
     RES.add(DphiDt_sys_rhs_2);
     RES*=-1.0;

     //for (unsigned int i=0; i<comp_dom.dh.n_dofs();i++)
     //    {
     //    if (constraints.is_constrained(i) == 0)
     //       cout<<"*eta_dot("<<i<<") "<<DphiDt_sys_solution_2(i)<<"   res("<<i<<") "<<RES(i)<<"  eta_res("<<i<<") "<<eta_res[i]<<endl;
     //    }
//*/

/*
     SparseDirectUMFPACK DphiDt_direct;
     SparseDirectUMFPACK DphiDt_direct_2;
     SparseDirectUMFPACK DphiDt_direct_3;
     
     DphiDt_direct.initialize(DphiDt_sys_matrix);
     DphiDt_direct_2.initialize(DphiDt_sys_matrix);
     DphiDt_direct_3.initialize(DphiDt_sys_matrix_2);

     DphiDt_direct.vmult(DphiDt_sys_solution, DphiDt_sys_rhs);
     constraints.distribute(DphiDt_sys_solution);
     DphiDt_direct_2.vmult(DphiDt_sys_solution_2, DphiDt_sys_rhs_2);
     constraints.distribute(DphiDt_sys_solution_2);
     DphiDt_direct_3.vmult(DphiDt_sys_solution_3, DphiDt_sys_rhs_3);
     constraints.distribute(DphiDt_sys_solution_3);

     Vector<double> test_phi(comp_dom.dh.n_dofs());
     DphiDt_sys_matrix.vmult(test_phi,(const Vector<double>&)phi_time_derivs);

     Vector<double> eta_dot(comp_dom.dh.n_dofs());
     for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
         eta_dot(i) = nodes_velocities(3*i+2);
     Vector<double> test_eta(comp_dom.dh.n_dofs());

     Vector<double> test_dphi_dn(comp_dom.dh.n_dofs());
     DphiDt_sys_matrix_2.vmult(test_dphi_dn,(const Vector<double>&)dphi_dn);     
   
     //DphiDt_sys_matrix.vmult(test_phi,(const Vector<double>&)phi_time_derivs);     
     DphiDt_sys_matrix.vmult(test_eta,eta_dot);  
*/ 
     //for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
     //    if ( (sys_comp(i+comp_dom.dh.n_dofs()+comp_dom.vector_dh.n_dofs()) == 4) )//&& (fabs(test_phi(i)-DphiDt_sys_rhs(i)-phi_res[i])>1e-15))
     //       cout<<i<<" "<<test_dphi_dn(i)-DphiDt_sys_rhs_3(i)<<"     "<<dphi_dn_res[i]<<"   ("<<test_dphi_dn(i)-DphiDt_sys_rhs_3(i)-dphi_dn_res[i]<<")"<<endl;

            //cout<<i<<" "<<test_phi(i)-DphiDt_sys_rhs(i)<<"     "<<phi_res[i]<<"   ("<<test_phi(i)-DphiDt_sys_rhs(i)-phi_res[i]<<")"<<endl;
     //for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
     //    if ( (sys_comp(3*i+2) == 1) )//&& (fabs(test_eta(i)-DphiDt_sys_rhs_2(i)-eta_res[i])>1e-15))
     //       {
     //       cout<<i<<" "<<test_eta(i)<<"     "<<DphiDt_sys_rhs_2(i)<<"   ("<<test_eta(i)-DphiDt_sys_rhs_2(i)<<")"<<endl;
            //cout<<i<<" "<<test_eta(i)-DphiDt_sys_rhs_2(i)<<"     "<<eta_res[i]<<"   ("<<test_eta(i)-DphiDt_sys_rhs_2(i)-eta_res[i]<<")"<<endl;
            //cout<<DphiDt_sys_solution_2(i)<<" vs "<<eta_dot(i)<<"  --->   "<<DphiDt_sys_solution_2(i)-eta_dot(i)<<endl;
     //       }
     double residual_counter_0 = 0;
     double residual_counter_1 = 0;
     double residual_counter_2 = 0;
     double residual_counter_3 = 0;
     double residual_counter_4 = 0;
     double residual_counter_5 = 0;
     double residual_counter_6 = 0;
     double residual_counter_7 = 0;
     double residual_counter_8 = 0;
     double residual_counter_9 = 0;
     double residual_counter_10 = 0;
     double residual_counter_11 = 0;
     double residual_counter_12 = 0;

     for(unsigned int i=0; i<sys_comp.size(); ++i)
        switch((int) sys_comp(i)) 
          {
	  case 0:
					       // Alg. X
	         dst(i) = nodes_pos_res(i);
                 //if ((fabs(dst(i)) > 1e-4) && fabs(t-1.1)<1e-6)
                 //   std::cout<<"i --> "<<i<<" (0) "<<dst(i)<<" "<<comp_dom.vector_support_points[i]<<std::endl;
                 residual_counter_0 += fabs(dst(i)*dst(i));
	         break;
	  case 1:
	         dst(i) = eta_res[int((i-2)/3)];
                 //std::cout<<i<<" --> "<<int((i-2)/3)<<" (1) "<<dst(i)<<" "<<comp_dom.support_points[int((i-2)/3)]<<std::endl;
	         residual_counter_1 += fabs(dst(i)*dst(i));
	         break;
	  case 2:
	         dst(i) = bem_phi(i-comp_dom.vector_dh.n_dofs()) - src_yy(i);
                 //std::cout<<"i --> "<<i<<" (2) "<<dst(i)<<std::endl;
	         residual_counter_2 += fabs(dst(i)*dst(i));
	         break;
	  case 3:
	         dst(i) = phi_res[i-comp_dom.vector_dh.n_dofs()];
                 //if (fabs(phi_res[i-comp_dom.vector_dh.n_dofs()]) > 1e-5)
                    //std::cout<<"i --> "<<i<<" (3) "<<dst(i)<<" "<<comp_dom.support_points[i-comp_dom.vector_dh.n_dofs()]<<std::endl;
	         residual_counter_3 += fabs(dst(i)*dst(i));
	         break;
	  case 4:
	         dst(i) = //bem_dphi_dn(i-comp_dom.vector_dh.n_dofs()-comp_dom.dh.n_dofs()) - src_yy(i);
                          dphi_dn_res[i-comp_dom.vector_dh.n_dofs()-comp_dom.dh.n_dofs()];
                          //DphiDt_sys_solution_3(i-comp_dom.vector_dh.n_dofs()-comp_dom.dh.n_dofs()) - src_yy(i);
                 //std::cout<<"i --> "<<i<<" (4) "<<dst(i)<<std::endl;
	         residual_counter_4 += fabs(dst(i)*dst(i));
	         break;
	  case 5:
	         dst(i) = bem_dphi_dn(i-comp_dom.vector_dh.n_dofs()-comp_dom.dh.n_dofs()) - src_yy(i);
                 //std::cout<<"i --> "<<i<<" (5) "<<dst(i)<<std::endl;
	         residual_counter_5 += fabs(dst(i)*dst(i));
	         break;
   	  case 6:
	         dst(i) = bem_phi(i-comp_dom.vector_dh.n_dofs()) - src_yy(i);
                 //std::cout<<"i --> "<<i<<" (6) "<<dst(i)<<std::endl;
	         residual_counter_6 += fabs(dst(i)*dst(i));
	         break;
          case 7:
	         dst(i) = bem_dphi_dn(i-comp_dom.vector_dh.n_dofs()-comp_dom.dh.n_dofs()) - src_yy(i);
                          //bem_residual(i-comp_dom.vector_dh.n_dofs()-comp_dom.dh.n_dofs());
                 //std::cout<<"i --> "<<i<<" (7) "<<dst(i)<<std::endl;
	         residual_counter_7 += fabs(dst(i)*dst(i));
	         break;
   	  case 8:
	         dst(i) = x_smoothing_res[int((i)/3)];
                 //if (dst(i)>1e-6)
                 //   std::cout<<"i --> "<<i<<" (8) "<<dst(i)<<" *** "<<comp_dom.vector_support_points[i]<<std::endl;
	         residual_counter_8 += fabs(dst(i)*dst(i));
	         break;
          case 9:
	        dst(i) = y_smoothing_res[int((i-1)/3)]; 
                 //std::cout<<"i --> "<<i<<" (9) "<<dst(i)<<std::endl;
	         residual_counter_9 += fabs(dst(i)*dst(i));
	         break;
   	  case 10:
	         dst(i) = bem_phi(i-comp_dom.vector_dh.n_dofs()) - src_yy(i);
                 //std::cout<<"i --> "<<i<<" (10) "<<dst(i)<<" ("<<bem_phi(i-comp_dom.vector_dh.n_dofs())<<" vs "<<src_yy(i)<<")"<<std::endl;
	         residual_counter_10 += fabs(dst(i)*dst(i));
                 break;
   	  case 11:
                 //cout<<i<<" ++++++ "<<i-comp_dom.vector_dh.n_dofs()-comp_dom.dh.n_dofs()-comp_dom.dh.n_dofs()<<endl;
	         dst(i) = hull_rigid_modes_res(i-comp_dom.vector_dh.n_dofs()-comp_dom.dh.n_dofs()-comp_dom.dh.n_dofs());
                 //std::cout<<"i --> "<<i<<" (11) "<<dst(i)<<" ("<<bem_phi(i-comp_dom.vector_dh.n_dofs())<<" vs "<<src_yy(i)<<")"<<std::endl;
	         residual_counter_11 += fabs(dst(i)*dst(i));
                 break;
	  default:
	         Assert(false, ExcInternalError());
	         break;
          }
     std::cout << "Current total residual: " << dst.l2_norm() << std::endl;
     std::cout << "Current generic algebraic coords residual: " << sqrt(residual_counter_0) << std::endl;
     std::cout << "Current differential coords residual: " << sqrt(residual_counter_1) << std::endl;
     std::cout << "Current algebraic potential residual: " << sqrt(residual_counter_2) << std::endl;
     std::cout << "Current differential potential residual: " << sqrt(residual_counter_3) << std::endl;
     std::cout << "Current algebraic potential normal gradient residual: " << sqrt(residual_counter_4) << std::endl;
     std::cout << "Current algebraic potential normal gradient hanging nodes residual: " << sqrt(residual_counter_5) << std::endl;
     std::cout << "Current algebraic potential hanging nodes residual: " << sqrt(residual_counter_6) << std::endl;
     std::cout << "Current algebraic potential normal gradient bem residual: " << sqrt(residual_counter_7) << std::endl;
     std::cout << "Current algebraic free surface coords x residual: " << sqrt(residual_counter_8) << std::endl;
     std::cout << "Current algebraic free surface coords y residual: " << sqrt(residual_counter_9) << std::endl;
     std::cout << "Current algebraic pressure/free surf phi residual: " << sqrt(residual_counter_10) << std::endl;
     std::cout << "Current differential rigid mode vel residual: " << sqrt(residual_counter_11) << std::endl;
     std::cout << "Current differential rigid mode displ residual: " << sqrt(residual_counter_12) << std::endl;    
     static unsigned int res_evals = 0;
     res_evals++;
     std::cout << "Residual evaluations: " << res_evals << std::endl;



  return 0;

}   

				     /** Jacobian preconditioner
					 vector product for Newton. */
template <int dim>
int FreeSurface<dim>::jacobian_prec_prod(Vector<double> &dst,  
			                 const Vector<double> &src_yy,
			                 const Vector<double> &src)
{
Vector<double> src_copy(src);
jacobian_preconditioner_matrix.vmult(dst,src_copy);

  return 0;
}


				     /** Jacobian inverse preconditioner
					 vector product for Newton. */
template <int dim>
int FreeSurface<dim>::jacobian_prec(Vector<double> &dst,  
			            const Vector<double> &src_yy,
			            const Vector<double> &src)
{

jacobian_prec(current_time,dst,src_yy,current_sol_dot,src,1e7);
//for (unsigned int i=0; i<dst.size();++i)
//    {
//    cout<<" *** "<<i<<" "<<src_yy(i)<<" "<<src(i)<<" "<<dst(i)<<" "<<current_sol_dot(i)<<endl;
//    }
  return 0;
}


				     /** Jacobian inverse preconditioner
					 vector product for DAE. */
template <int dim>
int FreeSurface<dim>::jacobian_prec(const double t,
			    Vector<double> &dst,  
			    const Vector<double> &src_yy,
			    const Vector<double> &src_yp,
			    const Vector<double> &src,
			    const double alpha)    
{

cout<<"alpha (using) "<<alpha<<endl;
SparseDirectUMFPACK prec_direct;
prec_direct.initialize(jacobian_preconditioner_matrix);
prec_direct.vmult(dst,src);

Vector<double> residual(dst.size());
jacobian_preconditioner_matrix.vmult(residual,dst);
residual *= -1.0;
residual.add(src);



  //SolverGMRES<Vector<double> > solver (solver_control,
  //SolverGMRES<Vector<double> >::AdditionalData(100));
  //solver.solve (jacobian_preconditioner_matrix, dst, src, preconditioner_preconditioner_matrix);


//  cout<<"----Rullo di tamburi--------"<<endl;
//  SolverGMRES<Vector<double> > solver (solver_control,
//  SolverGMRES<Vector<double> >::AdditionalData(100));
//  solver.solve (*this, dst, src, jacobian_preconditioner_matrix);
//  cout<<"----Respirate, prego--------"<<endl; 

cout<<"Inverting preconditioner"<<endl;

return 0;
}    

				     /** Setup Jacobian preconditioner for Newton. */
template <int dim>
int FreeSurface<dim>::setup_jacobian_prec(const Vector<double> &src_yy)
{
setup_jacobian_prec(current_time,src_yy,current_sol_dot,0.0);
//cout<<current_sol_dot<<endl;
return 0;
}   
				     /** Setup Jacobian preconditioner for DAE. */
template <int dim>
int FreeSurface<dim>::setup_jacobian_prec(const double t,
				  const Vector<double> &src_yy,
				  const Vector<double> &src_yp,
				  const double alpha)   
{
cout<<"alpha (building) "<<alpha<<endl;
jacobian_preconditioner_matrix = 0;
preconditioner_preconditioner_matrix = 0;
//int preconditioner_band = 100;
for (unsigned int i=0; i<this->n_dofs(); ++i)
    for (SparsityPattern::iterator col=jacobian_sparsity_pattern.begin(i); col!=jacobian_sparsity_pattern.end(i); ++col)
        {
        unsigned int j = col->column();
        //cout<<" "<<i<<" "<<j<<" "<<jacobian_matrix(i,j)+alpha*jacobian_dot_matrix(i,j)<<endl;
        //if ( (j > (int)i-preconditioner_band) &&
        //     (j > (int)i+preconditioner_band)   ) 
           jacobian_preconditioner_matrix.set(i,j,jacobian_matrix(i,j)+alpha*jacobian_dot_matrix(i,j));
        }

for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
    {
    if ( (!constraints.is_constrained(i)) &&
         ((comp_dom.flags[i] & water) || (comp_dom.flags[i] & pressure)) &&
         (!(comp_dom.flags[i] & transom_on_water))   )
       {
       jacobian_preconditioner_matrix.set(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs(),
                                          i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs(),
                                          1.0);
       //cout<<"   "<<i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()<<endl;
       }
    else if ( (!constraints.is_constrained(i)) &&
            (!(comp_dom.flags[i] & water)) )
       {
       jacobian_preconditioner_matrix.set(i+comp_dom.vector_dh.n_dofs(),
                                          i+comp_dom.vector_dh.n_dofs(),
                                          1.0);
       //cout<<"   "<<i+comp_dom.vector_dh.n_dofs()<<endl;
       }
    else if (comp_dom.flags[i] & transom_on_water)
       {
       jacobian_preconditioner_matrix.set(i+comp_dom.vector_dh.n_dofs(),
                                          i+comp_dom.vector_dh.n_dofs(),
                                          1.0);
       //cout<<"   "<<i+comp_dom.vector_dh.n_dofs()<<endl;
       }
    }



cout<<"Building preconditioner"<<endl;

/*
//if (fabs(t-0.08)<0.00001)
{
for (unsigned int i=0; i<this->n_dofs(); ++i)
    for (SparsityPattern::iterator c=jacobian_sparsity_pattern.begin(i); c!=jacobian_sparsity_pattern.end(i); ++c)
        {
        unsigned int j = c->column();
        cout<<" "<<i<<" "<<j<<" "<<jacobian_preconditioner_matrix(i,j)<<" "<<jacobian_matrix(i,j)<<endl;
        }
}
//*/
return 0;
} 


				     /** And an identification of the
					 differential components. This
					 has to be 1 if the
					 corresponding variable is a
					 differential component, zero
					 otherwise.  */
template <int dim>
Vector<double> & FreeSurface<dim>::differential_components()
{
  if(diff_comp.size() != this->n_dofs())
    {

      diff_comp.reinit(this->n_dofs());
      alg_comp.reinit(this->n_dofs());
      sys_comp.reinit(this->n_dofs());

				       // 0: X algebraic (horizontal coors, vertical constrained coors,
                                       //                 vertical coors not on free surface)
                                       // 1: X differential (vertical non constrained free surface coors) 
				       // 2: phi algebraic (on neumann boundaries, non constrained nodes)
                                       // 3: phi differential (non constrained nodes on dirichlet pressure boundaries / free surface)
				       // 4: dphi/dn on boat or inflow coming from L2 projection (algebraic)
                                       // 5: dphi/dn algebraic (constrained nodes)
                                       // 6: phi algebraic (on constrained nodes)
                                       // 7: dphi/dn algebraic on water/pressure coming from BEM solution (non constrained nodes)
                                       // 8: X algebraic coming from x position smoothing
                                       // 9: X algebraic coming from y position smoothing
                                       // 10: phi algebraic (imposed equal to water counterpart on pressure side of pressure/water interface)
                                       // 11: rigid body velocities and displacements (differential)  
      
      //let's start with the nodes positions dofs
      for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
	  {
          if (comp_dom.vector_flags[3*i] & water)
             // if node is on free surface, nodes have first two components set to algebraic,
             // third set to differential (and sys_comp = 1)
             {
             // we distinguish the cases of internal (x and y given by smoothing)
             // and edge (x and y given by geometric treatment) nodes
             if ( !(comp_dom.vector_flags[3*i] & edge) )
                {
                // first component is algebraic and marked with 8 in sys_comp
                sys_comp(3*i) = 8;
                diff_comp(3*i) = 0;
                alg_comp(3*i) = 1;
                // second component is algebraic and marked with 9 in sys_comp
                sys_comp(3*i+1) = 9;
                diff_comp(3*i+1) = 0;
                alg_comp(3*i+1) = 1;
                // only third component is differential and marked 1 in sys_comp
                sys_comp(3*i+2) = 1;
                diff_comp(3*i+2) = 1;
                alg_comp(3*i+2) = 0;
                }
             // all other are edges dofs and have first and second algebraic components, third is differential
             else
                {// only exception is when node is a transom_on_water and near_inflow nodes. in such case it's all algebraic
                if (comp_dom.vector_flags[3*i] & transom_on_water ||
                    comp_dom.vector_flags[3*i] & near_inflow)
                   {
                   for (unsigned int j=0; j<dim; ++j)
                       {
                       sys_comp(3*i+j) = 0;
                       diff_comp(3*i+j) = 0;
                       alg_comp(3*i+j) = 1;
                       }
                   }
                else
                   {
                   // first component is algebraic and marked with 0 in sys_comp
                   sys_comp(3*i) = 0;
                   diff_comp(3*i) = 0;
                   alg_comp(3*i) = 1;
                   // second component is algebraic and marked with 0 in sys_comp
                   sys_comp(3*i+1) = 0;
                   diff_comp(3*i+1) = 0;
                   alg_comp(3*i+1) = 1;
                   // only third component is differential and marked 1 in sys_comp
                   sys_comp(3*i+2) = 1;
                   diff_comp(3*i+2) = 1;
                   alg_comp(3*i+2) = 0;
                   }
                }
             }
          else
             // if node is not on free surface, it's algebraic for sure
             {
             for (unsigned int j=0; j<dim-1; ++j)
                 {
                 sys_comp(3*i+j) = 0;
                 diff_comp(3*i+j) = 0;
                 alg_comp(3*i+j) = 1;
                 }
             }
          }


      // let's take care of the potential dofs
      for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
	  {
          if (comp_dom.flags[i] & water)
             // if node is on free surface, component is differential
             {// unless node is a transom_on_water node. then it's algebraic
             if (comp_dom.flags[i] & transom_on_water)
                {
                sys_comp(i+comp_dom.vector_dh.n_dofs()) = 2;
                diff_comp(i+comp_dom.vector_dh.n_dofs()) = 0;
                alg_comp(i+comp_dom.vector_dh.n_dofs()) = 1; 
                }
             else
                {
                diff_comp(i+comp_dom.vector_dh.n_dofs()) = 1;
                alg_comp(i+comp_dom.vector_dh.n_dofs()) = 0;
                sys_comp(i+comp_dom.vector_dh.n_dofs()) = 3; 
                }
             }
          else
          // if node is not on free surface, it's algebraic for sure, unless it's a pressure node
          // in that case potential is differential 
             {
             if ( comp_dom.flags[i] & pressure )
                {
                diff_comp(i+comp_dom.vector_dh.n_dofs()) = 1;
                alg_comp(i+comp_dom.vector_dh.n_dofs()) = 0;
                sys_comp(i+comp_dom.vector_dh.n_dofs()) = 3; 
                // pressure node on the free surface must be treated differently
                if (comp_dom.flags[i] & near_water)
                   {
                   diff_comp(i+comp_dom.vector_dh.n_dofs()) = 0;
                   alg_comp(i+comp_dom.vector_dh.n_dofs()) = 1;
                   sys_comp(i+comp_dom.vector_dh.n_dofs()) = 10; 
                   }
                }
             else
                {
                diff_comp(i+comp_dom.vector_dh.n_dofs()) = 0;
                alg_comp(i+comp_dom.vector_dh.n_dofs()) = 1;
                sys_comp(i+comp_dom.vector_dh.n_dofs()) = 2;
                }
             } 
          }

       // let's take care of the potential normal derivative dofs
      for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
	  {// in this case all dofs are algebraic components
          if (!(comp_dom.flags[i] & water) && !(comp_dom.flags[i] & pressure))
             {// neumann boundary condition on boat surface is treated separately
             if (constraints.is_constrained(i) )
	        {// if node is constrained, it is component 5
                diff_comp(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()) = 0;
                alg_comp(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()) = 1;
                sys_comp(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()) = 5;             
                }
             else
                {// otherwise, it is component 4
                diff_comp(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()) = 0;
                alg_comp(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()) = 1;
                sys_comp(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()) = 4;
                }
             }
          else
             {
             if (constraints.is_constrained(i) )
	        {// if node is constrained, it is component 5 again
                diff_comp(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()) = 0;
                alg_comp(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()) = 1;
                sys_comp(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()) = 5;             
                }
             else
                {// otherwise, it is component 7
                diff_comp(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()) = 0;
                alg_comp(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()) = 1;
                sys_comp(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()) = 7;
                }
             }         
          }

	
      // all hanging (constrained) nodes are algebraic
      for (unsigned int i=0; i<comp_dom.vector_dh.n_dofs(); i++)
	{
	  if (vector_constraints.is_constrained(i))
	    {
	      diff_comp(i) = 0;
	      sys_comp(i) = 0;
              alg_comp(i) = 1;
	    }
	}

      // all hanging (constrained) nodes are algebraic
      for (unsigned int i=0; i<comp_dom.dh.n_dofs(); i++)
	{
	  if (constraints.is_constrained(i))
	    {
	      diff_comp(i+comp_dom.vector_dh.n_dofs()) = 0;
	      sys_comp(i+comp_dom.vector_dh.n_dofs()) = 6;
              alg_comp(i+comp_dom.vector_dh.n_dofs()) = 1;
	    }
	}

  // finally, we set the final 13 components, needed to displace the boat hull
  for (unsigned int i=0; i<13; ++i)
     {
     diff_comp(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()) = 1;
     sys_comp(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()) = 11;
     alg_comp(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()) = 0;     
     }


  std::vector<unsigned int> row_lengths(dofs_number);

  for (unsigned int i=0; i<dofs_number; ++i)
      row_lengths[i] = 5*comp_dom.dh.max_couplings_between_dofs();

  row_lengths[comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()] = dofs_number;
  row_lengths[comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+1] = dofs_number;
  row_lengths[comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+2] = dofs_number;
  row_lengths[comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+3] = 2;
  row_lengths[comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+4] = 2;
  row_lengths[comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+5] = 2;
  row_lengths[comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+6] = dofs_number;
  row_lengths[comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+7] = dofs_number;
  row_lengths[comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+8] = dofs_number;
  row_lengths[comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+9] = 7;
  row_lengths[comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+10] = 7;
  row_lengths[comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+11] = 7;
  row_lengths[comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+12] = 7;

  

  //ConstraintMatrix constr;
  //constr.clear();
  //constr.close();   

                                      // we initialize here the jacobian
                                      // sparsity patterns

  jacobian_sparsity_pattern.reinit(dofs_number,
                                   dofs_number,
                                   row_lengths);
  
  cell_it
  cell = comp_dom.dh.begin_active(),
  endc = comp_dom.dh.end();

  const unsigned int  dofs_per_cell   = comp_dom.fe.dofs_per_cell;
  std::vector<unsigned int> local_dof_indices(dofs_per_cell);
   
  for (; cell!=endc; ++cell)
      {
      cell->get_dof_indices(local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
          unsigned int ii = local_dof_indices[i];
          if (comp_dom.flags[ii] & transom_on_water)
             {
             std::set<unsigned int> duplicates = comp_dom.double_nodes_set[ii];
             jacobian_sparsity_pattern.add(3*ii,3*ii);
             jacobian_sparsity_pattern.add(3*ii+1,3*ii+1);
             jacobian_sparsity_pattern.add(3*ii+2,3*ii+2);
             for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos!=duplicates.end();++pos)
                 {
                 unsigned int jj=*pos;
                 jacobian_sparsity_pattern.add(ii+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs(),jj+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs());
                 }
             }

          for (unsigned int j=0; j<dofs_per_cell; ++j)
              {
              unsigned int jj = local_dof_indices[j];
              jacobian_sparsity_pattern.add(3*ii,3*jj);
              jacobian_sparsity_pattern.add(3*ii,3*jj+1);
              jacobian_sparsity_pattern.add(3*ii,3*jj+2);
              jacobian_sparsity_pattern.add(3*ii,jj+comp_dom.vector_dh.n_dofs());
              jacobian_sparsity_pattern.add(3*ii,jj+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs());

              jacobian_sparsity_pattern.add(3*ii+1,3*jj);
              jacobian_sparsity_pattern.add(3*ii+1,3*jj+1);
              jacobian_sparsity_pattern.add(3*ii+1,3*jj+2);
              jacobian_sparsity_pattern.add(3*ii+1,jj+comp_dom.vector_dh.n_dofs());
              jacobian_sparsity_pattern.add(3*ii+1,jj+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs());

              jacobian_sparsity_pattern.add(3*ii+2,3*jj);
              jacobian_sparsity_pattern.add(3*ii+2,3*jj+1);
              jacobian_sparsity_pattern.add(3*ii+2,3*jj+2);
              jacobian_sparsity_pattern.add(3*ii+2,jj+comp_dom.vector_dh.n_dofs());
              jacobian_sparsity_pattern.add(3*ii+2,jj+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs());

              jacobian_sparsity_pattern.add(ii+comp_dom.vector_dh.n_dofs(),3*jj);
              jacobian_sparsity_pattern.add(ii+comp_dom.vector_dh.n_dofs(),3*jj+1);
              jacobian_sparsity_pattern.add(ii+comp_dom.vector_dh.n_dofs(),3*jj+2);
              jacobian_sparsity_pattern.add(ii+comp_dom.vector_dh.n_dofs(),jj+comp_dom.vector_dh.n_dofs());
              jacobian_sparsity_pattern.add(ii+comp_dom.vector_dh.n_dofs(),jj+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs());

              jacobian_sparsity_pattern.add(ii+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs(),3*jj);
              jacobian_sparsity_pattern.add(ii+comp_dom.dh.n_dofs()+comp_dom.vector_dh.n_dofs(),3*jj+1);
              jacobian_sparsity_pattern.add(ii+comp_dom.dh.n_dofs()+comp_dom.vector_dh.n_dofs(),3*jj+2);
              jacobian_sparsity_pattern.add(ii+comp_dom.dh.n_dofs()+comp_dom.vector_dh.n_dofs(),jj+comp_dom.vector_dh.n_dofs());
              jacobian_sparsity_pattern.add(ii+comp_dom.dh.n_dofs()+comp_dom.vector_dh.n_dofs(),jj+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()); 
              //cout<<ii+comp_dom.dh.n_dofs()+comp_dom.vector_dh.n_dofs()<<","<<3*jj<<endl;
              //cout<<ii+comp_dom.dh.n_dofs()+comp_dom.vector_dh.n_dofs()<<","<<3*jj+1<<endl;
              //cout<<ii+comp_dom.dh.n_dofs()+comp_dom.vector_dh.n_dofs()<<","<<3*jj+2<<endl;             
              }
          // this is for the dependence of positions, phi and dphi_dn on the rigid modes 
          for (unsigned int d=0; d<13; ++d)
              {
              for (unsigned int k=0; k<3; ++k)
                  {
                  jacobian_sparsity_pattern.add(3*ii+k,d+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs());
                  }
              jacobian_sparsity_pattern.add(ii+comp_dom.vector_dh.n_dofs(),
                                            d+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs());
              jacobian_sparsity_pattern.add(ii+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs(),
                                            d+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs());
              }
          }
      }

    for (unsigned int i=0; i<comp_dom.vector_dh.n_dofs(); ++i)
      {
      if (vector_constraints.is_constrained(i))
         {
         jacobian_sparsity_pattern.add(i,i);
         std::vector< std::pair< unsigned int, double > > entries = *vector_constraints.get_constraint_entries(i);
         for (unsigned int k=0; k<entries.size(); ++k)
             {
             jacobian_sparsity_pattern.add(i,entries[k].first);
             }
         }
      }

    for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
      {
      if (constraints.is_constrained(i))
         {
         jacobian_sparsity_pattern.add(i+comp_dom.vector_dh.n_dofs(),i+comp_dom.vector_dh.n_dofs());
         jacobian_sparsity_pattern.add(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs(),i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs());
         std::vector< std::pair< unsigned int, double > > entries = *constraints.get_constraint_entries(i);
         for (unsigned int k=0; k<entries.size(); ++k)
             {
             jacobian_sparsity_pattern.add(i+comp_dom.vector_dh.n_dofs(),entries[k].first+comp_dom.vector_dh.n_dofs());
             jacobian_sparsity_pattern.add(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs(),
                                           entries[k].first+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs());
             }
         }
      }

     //this takes care of the left water line nodes projection (without smoothing)
     if (!comp_dom.no_boat)
        {
        for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
            { 
            if ( (comp_dom.flags[i] & water) &&
                 (comp_dom.flags[i] & near_boat) &&
                 (comp_dom.flags[i] & left_side) &&
                 (comp_dom.moving_point_ids[0] != i) &&
                 (comp_dom.moving_point_ids[1] != i) &&
                 (comp_dom.moving_point_ids[2] != i) ) // to avoid the bow and stern node
               {//cout<<"**** "<<i<<endl;
   
               jacobian_sparsity_pattern.add(3*i,3*i);
               jacobian_sparsity_pattern.add(3*i,3*i+1);
               jacobian_sparsity_pattern.add(3*i,3*i+2);
               jacobian_sparsity_pattern.add(3*i+1,3*i+1);
               jacobian_sparsity_pattern.add(3*i+1,3*i+2);
               jacobian_sparsity_pattern.add(3*i+1,3*i);
               //cout<<i<<"   "<<temp_src(3*i+1)<<"   ("<<comp_dom.iges_normals[i]<<")"<<endl;
               }              
            }
        //this takes care of the bow and stern nodes 
        for (unsigned int k=0; k<3; ++k)
            { 
            unsigned int i = comp_dom.moving_point_ids[k];

            jacobian_sparsity_pattern.add(3*i,3*i);
            jacobian_sparsity_pattern.add(3*i,3*i+2); 
            jacobian_sparsity_pattern.add(3*i+1,3*i+1);
            jacobian_sparsity_pattern.add(3*i+1,3*i+2);
               //cout<<i<<" (point) "<<comp_dom.support_points[i]<<endl;
               //cout<<i<<" (edges_tangents) "<<comp_dom.edges_tangents(3*i)<<","<<comp_dom.edges_tangents(3*i+1)<<","<<comp_dom.edges_tangents(3*i+2)<<endl;
            }              
         }

     // this cycle hooks the boat and far field double nodes
     // to their water twins that have been moved
     for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
         {
         if ( (comp_dom.flags[i] & water) &&
              (comp_dom.flags[i] & edge) )
            {//cout<<"c "<<i<<endl;
            std::set<unsigned int> duplicates = comp_dom.double_nodes_set[i];
            duplicates.erase(i); 
            for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
                {
                for (unsigned int k=0; k<3; ++k)
                    {
                    jacobian_sparsity_pattern.add(3*(*pos)+k,3*(*pos)+k);
                    if (comp_dom.flags[i] & transom_on_water)
                       {
                       for (unsigned int d=0; d<3; ++d)
                           jacobian_sparsity_pattern.add(3*(*pos)+k,3*i+d);
                       for (unsigned int d=0; d<3; ++d)
                           jacobian_sparsity_pattern.add(3*(*pos)+k,comp_dom.vector_dh.n_dofs()+
                                                                    comp_dom.dh.n_dofs()+
                                                                    comp_dom.dh.n_dofs()+3+d);
                       for (unsigned int d=0; d<4; ++d)
                           jacobian_sparsity_pattern.add(3*(*pos)+k,comp_dom.vector_dh.n_dofs()+
                                                                    comp_dom.dh.n_dofs()+
                                                                    comp_dom.dh.n_dofs()+9+d);
                       }
                    else
                       {
                       jacobian_sparsity_pattern.add(3*(*pos)+k,3*i+k);
                       }
                    }
                }
            }
         }

     // this cycle is to impose the value of potential on the pressure
     // side of pressure/free surface interface 
     for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
         {
         if ( (comp_dom.flags[i] & pressure) &&
              (comp_dom.flags[i] & near_water) )
            {//cout<<"c "<<i<<endl;
            std::set<unsigned int> duplicates = comp_dom.double_nodes_set[i];
            duplicates.erase(i); 
            for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
                {
                if (comp_dom.flags[*pos] & water)
                   {
                   jacobian_sparsity_pattern.add(comp_dom.vector_dh.n_dofs()+i,comp_dom.vector_dh.n_dofs()+ *pos);
                   jacobian_sparsity_pattern.add(comp_dom.vector_dh.n_dofs()+i,comp_dom.vector_dh.n_dofs()+i);
                   }
                }
            }
         }

  // here's the rigid body part of the dofs jacobian sparsity pattern
  for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
      {
      if ( comp_dom.flags[i] & boat )
         {
         for (unsigned int k=0; k<3; ++k)
             {
             jacobian_sparsity_pattern.add(k+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                           3*i);
             jacobian_sparsity_pattern.add(k+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                           3*i+1);
             jacobian_sparsity_pattern.add(k+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                           3*i+2);
             jacobian_sparsity_pattern.add(k+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                           i+comp_dom.vector_dh.n_dofs());
             jacobian_sparsity_pattern.add(k+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                           i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs());
             jacobian_sparsity_pattern.add(k+6+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                           3*i);
             jacobian_sparsity_pattern.add(k+6+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                           3*i+1);
             jacobian_sparsity_pattern.add(k+6+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                           3*i+2);
             jacobian_sparsity_pattern.add(k+6+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                           i+comp_dom.vector_dh.n_dofs());
             jacobian_sparsity_pattern.add(k+6+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                           i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs());

             }
         }
      } 

   for (unsigned int k=0; k<3; ++k)
       for (unsigned int d=0; d<13; ++d)
           {
           jacobian_sparsity_pattern.add(k+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                         d+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs());
           jacobian_sparsity_pattern.add(k+6+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                         d+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs());
           
           }

  for (unsigned int k=0; k<3; ++k)
      {
      jacobian_sparsity_pattern.add(k+3+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                    k+3+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs());
      jacobian_sparsity_pattern.add(k+3+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                    k+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs());
      }

  for (unsigned int k=0; k<4; ++k)
      {
      for (unsigned int d=0; d<7; ++d)
           {
           jacobian_sparsity_pattern.add(k+9+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),
                                         d+6+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs());
           }
      }

  jacobian_sparsity_pattern.compress();
  jacobian_matrix.reinit(jacobian_sparsity_pattern);
  jacobian_dot_matrix.reinit(jacobian_sparsity_pattern);
  jacobian_preconditioner_matrix.reinit(jacobian_sparsity_pattern);
  preconditioner_preconditioner_matrix.reinit(jacobian_sparsity_pattern); 


//*/

    }
  
  return diff_comp;
}






template <int dim>
void FreeSurface<dim>::prepare_bem_vectors(double time,
					   Vector<double> &bem_bc,
					   Vector<double> &dphi_dn) const
{

   // we first compute the vectors surface_nodes and other_nodes
   // with flags for the Dirichlet (free surface) and Neumann (other)
   // nodes for the bem_problem class
   cell_it
   cell = comp_dom.dh.begin_active(),
   endc = comp_dom.dh.end();
     
   comp_dom.surface_nodes.reinit(comp_dom.dh.n_dofs());
   comp_dom.other_nodes.reinit(comp_dom.dh.n_dofs());
   comp_dom.other_nodes.add(1);
   std::vector<unsigned int> dofs(comp_dom.fe.dofs_per_cell);
  
   for (; cell != endc; ++cell)
       {
       if (cell->material_id() == comp_dom.free_sur_ID1 ||
           cell->material_id() == comp_dom.free_sur_ID2 ||
           cell->material_id() == comp_dom.free_sur_ID3)
          {
          // This is a free surface node.
          cell->get_dof_indices(dofs);
          for (unsigned int i=0; i<comp_dom.fe.dofs_per_cell; ++i)
              {
	      comp_dom.surface_nodes(dofs[i]) = 1;
	      comp_dom.other_nodes(dofs[i]) = 0;
              }
          }
       else
          {
          for (unsigned int i=0; i<comp_dom.fe.dofs_per_cell; ++i)
              {
              cell->get_dof_indices(dofs);
              }
          }
       } 

   
   cell = comp_dom.dh.begin_active();


  std::vector<Point<dim> > support_points(comp_dom.dh.n_dofs());
  DoFTools::map_dofs_to_support_points<dim-1, dim>( *comp_dom.mapping, comp_dom.dh, support_points);
			   
  const unsigned int  dofs_per_cell   = comp_dom.fe.dofs_per_cell;
  std::vector<unsigned int> local_dof_indices (dofs_per_cell);   
    for(cell = comp_dom.dh.begin_active(); cell != endc; ++cell)
         {
	  cell->get_dof_indices(local_dof_indices);
  	  for(unsigned int j=0; j<comp_dom.fe.dofs_per_cell; ++j) 
	      {	 
	      if (cell->material_id() != comp_dom.free_sur_ID1 &&
	          cell->material_id() != comp_dom.free_sur_ID2 &&
		  cell->material_id() != comp_dom.free_sur_ID3 )
	          {
		  if (cell->material_id() == comp_dom.wall_sur_ID1 ||
		      cell->material_id() == comp_dom.wall_sur_ID2 ||
		      cell->material_id() == comp_dom.wall_sur_ID3)
		     {
                     Vector<double> wind_value(dim);
                     wind.vector_value(support_points[local_dof_indices[j]],wind_value);
                     Point<dim> Vinf;
	             for (unsigned int i = 0; i < dim; i++)
 	                 Vinf(i) = wind_value(i);
                     ////////////////////////////////////
                     // needs to be changed
                     //Point<3> original(support_points[local_dof_indices[j]](0),
                     //             fabs(support_points[local_dof_indices[j]](1)),
                     //                  support_points[local_dof_indices[j]](2));
                     //Point<3> projection;
                     //Point<3> normal;
                     //double mean_curvature;
                     //comp_dom.boat_model.boat_water_line_right->axis_projection_and_diff_forms(projection, normal, mean_curvature, original);
                     //double b = support_points[local_dof_indices[j]](1) < 0.0 ? -1.0 : 1.0; 
                     //projection(1)*=b;
                     //normal(1)*=b;
                     
                     //BoatSurface<3> boat_surface;
                     //Point<3> normal2 = boat_surface.HullNormal(support_points[local_dof_indices[j]]);
                     //std::cout<<std::endl;
                     //std::cout<<local_dof_indices[j]<<"  Point "<<support_points[local_dof_indices[j]]<<std::endl;
                     //std::cout<<"NormalNew "<<comp_dom.iges_normals[local_dof_indices[j]]<<std::endl;
                     //std::cout<<"NormalExact "<<normal2<<std::endl;
                     //std::cout<<"Error "<<normal2.distance(comp_dom.iges_normals[local_dof_indices[j]])<<std::endl;
                     /////////////////////////////////////
                     //bem_bc(local_dof_indices[j]) = -comp_dom.iges_normals[local_dof_indices[j]]*Vinf;
                     bem_bc(local_dof_indices[j]) = -comp_dom.node_normals[local_dof_indices[j]]*Vinf;
                     dphi_dn(local_dof_indices[j]) = bem_bc(local_dof_indices[j]);
                     //std::cout<<bem_bc(local_dof_indices[j])<<" "<<dphi_dn(local_dof_indices[j])<<" "<<Vinf<<std::endl;
		     }
		  else if (cell->material_id() == comp_dom.inflow_sur_ID1 ||
		           cell->material_id() == comp_dom.inflow_sur_ID2 ||
		           cell->material_id() == comp_dom.inflow_sur_ID3)
		     {
		     bem_bc(local_dof_indices[j]) = inflow_norm_potential_grad.value(support_points[local_dof_indices[j]]);
                     dphi_dn(local_dof_indices[j]) = inflow_norm_potential_grad.value(support_points[local_dof_indices[j]]);
                     cout<<inflow_norm_potential_grad.value(support_points[local_dof_indices[j]])<<endl;
		     }
                  else 
		     {
		     bem_bc(local_dof_indices[j]) = 0;
                     dphi_dn(local_dof_indices[j]) = 0;
		     }
	          }
	      else
	          {
                  
		  }

	      }

	
	 }

      // trying a fix for transom stern nodes
      Vector<double> wind_value(dim);
      wind.vector_value(Point<3>(0.0,0.0,0.0),wind_value);
      Point<dim> Vinf;
      for (unsigned int i = 0; i < dim; i++)
 	   Vinf(i) = wind_value(i);
      for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
          {
          if ( comp_dom.flags[i] & transom_on_water )
             {
	     comp_dom.surface_nodes(i) = 0;
	     comp_dom.other_nodes(i) = 1;
             std::set<unsigned int> duplicates = comp_dom.double_nodes_set[i];
             duplicates.erase(i); 
             bem_bc(i) = 0;
             for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
                 bem_bc(i) += bem_dphi_dn(*pos)/duplicates.size();
             dphi_dn(i) = bem_bc(i);
             }
          }   


}


				 // @sect4{BEMProblem::compute_gradient}

				 // The next function simply solves
				 // computes the gradient of the potential
                                 // (velocity field in fluid dynamics).
template <int dim>
void FreeSurface<dim>::compute_DXDt_and_DphiDt(double time,
                                               const Vector<double> & phi,
                                               const Vector<double> & dphi_dn,
                                               const Vector<double> & nodes_velocities)
{
   ConstraintMatrix constr;
   ConstraintMatrix vector_constr;

   constr.clear();
   vector_constr.clear();


   constr.close();
   vector_constr.close();
	
   std::vector<unsigned int> vector_dofs(comp_dom.vector_fe.dofs_per_cell); 
   std::vector<Point<dim> > vector_support_points(comp_dom.vector_dh.n_dofs());
   DoFTools::map_dofs_to_support_points<dim-1, dim>( *comp_dom.mapping, comp_dom.vector_dh, vector_support_points);

   DXDt_and_DphiDt_vector.reinit(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs());

   Vector<double> elevations(comp_dom.dh.n_dofs());  
   for (unsigned int i=0; i<comp_dom.dh.n_dofs(); i++)
       {
       elevations(i) = vector_support_points[dim*i](dim-1);
       }

   //for (unsigned int i=0; i<comp_dom.dh.n_dofs(); i++)
   //    {
   //    cout<<i<<"  "<<dphi_dn(i)<<endl;
   //    }

   DphiDt_sys_matrix.reinit (DphiDt_sparsity_pattern);
   DphiDt_sys_matrix_2.reinit (DphiDt_sparsity_pattern);
   vector_sys_matrix.reinit (vector_sparsity_pattern);
 
   DphiDt_sys_solution.reinit (comp_dom.dh.n_dofs());
   DphiDt_sys_solution_2.reinit (comp_dom.dh.n_dofs());
   DphiDt_sys_solution_3.reinit (comp_dom.dh.n_dofs());
   break_wave_press.reinit (comp_dom.dh.n_dofs());
   vector_sys_solution.reinit (comp_dom.vector_dh.n_dofs());
   vector_sys_solution_2.reinit (comp_dom.vector_dh.n_dofs());
   
   DphiDt_sys_rhs.reinit (comp_dom.dh.n_dofs()); 
   DphiDt_sys_rhs_2.reinit (comp_dom.dh.n_dofs()); 
   DphiDt_sys_rhs_3.reinit (comp_dom.dh.n_dofs()); 
   DphiDt_sys_rhs_4.reinit (comp_dom.dh.n_dofs()); 
   vector_sys_rhs.reinit (comp_dom.vector_dh.n_dofs());
   vector_sys_rhs_2.reinit (comp_dom.vector_dh.n_dofs());

   FEValues<dim-1,dim> vector_fe_v(*comp_dom.mapping, comp_dom.vector_fe, *comp_dom.quadrature,
			           update_values | update_gradients | 
			           update_cell_normal_vectors |
			           update_quadrature_points |
			           update_JxW_values);

   FEValues<dim-1,dim> fe_v(*comp_dom.mapping, comp_dom.fe, *comp_dom.quadrature,
			    update_values | update_gradients |
			    update_cell_normal_vectors |
			    update_quadrature_points |
			    update_JxW_values);

   const unsigned int DphiDt_n_q_points = fe_v.n_quadrature_points;
   const unsigned int  DphiDt_dofs_per_cell   = comp_dom.fe.dofs_per_cell;
   std::vector<unsigned int> DphiDt_local_dof_indices (DphiDt_dofs_per_cell);

   FullMatrix<double>   local_DphiDt_matrix (DphiDt_dofs_per_cell, DphiDt_dofs_per_cell);
   Vector<double>       local_DphiDt_rhs (DphiDt_dofs_per_cell);
   Vector<double>       local_DphiDt_rhs_2 (DphiDt_dofs_per_cell);
   
   const unsigned int vector_n_q_points = vector_fe_v.n_quadrature_points;
   const unsigned int   vector_dofs_per_cell   = comp_dom.vector_fe.dofs_per_cell;
   std::vector<unsigned int> vector_local_dof_indices (vector_dofs_per_cell);

   FullMatrix<double>   local_vector_matrix (vector_dofs_per_cell, vector_dofs_per_cell);
   Vector<double>       local_vector_rhs (vector_dofs_per_cell);
   Vector<double>       local_vector_rhs_2 (vector_dofs_per_cell);
   Vector<double>       local_vector_gradients (vector_dofs_per_cell);
   Vector<double>       local_vector_normals (vector_dofs_per_cell);


   std::vector< Tensor<1,dim> > vector_phi_surf_grads(DphiDt_n_q_points);
   std::vector< Tensor<1,dim> > vector_eta_surf_grads(DphiDt_n_q_points);
   std::vector<double> vector_phi_norm_grads(vector_n_q_points);
   std::vector<Vector<double> > quad_nodes_velocities(vector_n_q_points,
						  Vector<double>(dim));
   
   std::vector< Tensor<1,dim> > DphiDt_phi_surf_grads(vector_n_q_points);
   std::vector< Tensor<1,dim> > DphiDt_eta_surf_grads(vector_n_q_points);
   std::vector<double> DphiDt_phi_norm_grads(DphiDt_n_q_points);



   cell_it
   vector_cell = comp_dom.vector_dh.begin_active(),
   vector_endc = comp_dom.vector_dh.end();

   cell_it
   cell = comp_dom.dh.begin_active(),
   endc = comp_dom.dh.end();
   
         

   std::vector<Point<dim> > support_points(comp_dom.dh.n_dofs());
   DoFTools::map_dofs_to_support_points<dim-1, dim>( *comp_dom.mapping, comp_dom.dh, support_points);
   

   for (; cell!=endc,vector_cell!=vector_endc; ++cell,++vector_cell)
    {
     Assert(cell->index() == vector_cell->index(), ExcInternalError());
     //cout<<"Original "<<cell<<endl;   
     vector_fe_v.reinit (vector_cell);
     fe_v.reinit(cell);
     local_DphiDt_matrix = 0;
     local_DphiDt_rhs = 0;
     local_DphiDt_rhs_2 = 0;
     local_vector_matrix = 0;
     local_vector_rhs = 0;
     local_vector_rhs_2 = 0;
     //cell->get_dof_indices (DphiDt_local_dof_indices);
     //for (unsigned int i=0; i<DphiDt_dofs_per_cell; ++i)
         //{
         //phis[i] = phi(local_dof_indices[i]);
         //phis[i].diff(i+3*dofs_per_cell,5*dofs_per_cell);
         //dphi_dns[i] = dphi_dn(local_dof_indices[i]);
         //dphi_dns[i].diff(i+4*dofs_per_cell,5*dofs_per_cell);
         //std::cout<<i<<"--> "<<DphiDt_local_dof_indices[i]<<"--------->"<<phi(DphiDt_local_dof_indices[i])<<"  "<<dphi_dn(DphiDt_local_dof_indices[i])<<endl;
         //}
    
     
     fe_v.get_function_gradients((const Vector<double>&)phi, DphiDt_phi_surf_grads);
     fe_v.get_function_values(dphi_dn, DphiDt_phi_norm_grads);
     fe_v.get_function_gradients((const Vector<double>&)phi, vector_phi_surf_grads);
     fe_v.get_function_gradients(elevations, vector_eta_surf_grads);
     fe_v.get_function_values(dphi_dn, vector_phi_norm_grads);
     vector_fe_v.get_function_values(nodes_velocities, quad_nodes_velocities);

     const std::vector<Point<dim> > &DphiDt_node_normals = fe_v.get_normal_vectors();
     const std::vector<Point<dim> > &DphiDt_node_positions = fe_v.get_quadrature_points();
     const std::vector<Point<dim> > &vector_node_normals = vector_fe_v.get_normal_vectors();
     const std::vector<Point<dim> > &vector_node_positions = vector_fe_v.get_quadrature_points();
     
     double g = 9.81;
     
     std::vector<Vector<double> > DphiDt_wind_values(DphiDt_n_q_points);
     for (unsigned int q=0; q<DphiDt_n_q_points; ++q)
         {
	 DphiDt_wind_values[q].reinit(dim,true);
	 wind.vector_value(DphiDt_node_positions[q],DphiDt_wind_values[q]);
	 }
     std::vector<Vector<double> > vector_wind_values(DphiDt_n_q_points);
     for (unsigned int q=0; q<vector_n_q_points; ++q)
         {
	 vector_wind_values[q].reinit(dim,true);
	 wind.vector_value(vector_node_positions[q], vector_wind_values[q]);
         }

     unsigned int comp_i, comp_j;
     Point<dim> ez(0,0,1);
      for (unsigned int q=0; q<vector_n_q_points; ++q)
        {
	Point<dim> gradient = vector_node_normals[q]*vector_phi_norm_grads[q] + vector_phi_surf_grads[q];
        
   //std::cout<<gradient<<std::endl;
        //double nu_0 = 1.0;
        double nu_0 = 0;
        double damping = nu_0 *
                       pow((fmax(abs(vector_node_positions[q](0)),20)-20)/(31.169752-20),2) *
                       vector_phi_norm_grads[q];
//        double damping = nu_0 *
//                         pow((fmax(abs(vector_node_positions[q](0)),25)-25)/(25),2) *
//                         vector_node_positions[q](2);


	Point<dim> Vinf;
	for (unsigned int i = 0; i < dim; i++)
	    Vinf(i) = vector_wind_values[q](i);
	Point<dim> eta_grad = Point<dim>(0,0,0);
        //eta_grad(0) = -vector_node_normals[q](0)/vector_node_normals[q](2);
        //eta_grad(1) = -vector_node_normals[q](1)/vector_node_normals[q](2);
        eta_grad = eta_grad + vector_eta_surf_grads[q];
	eta_grad(dim-1) = 0;

        double fluid_vel_norm = (gradient+Vinf).norm();
        if (fluid_vel_norm < 1e-3)
            fluid_vel_norm = -8.0e+05*pow(fluid_vel_norm,3.0) + 1.7e+03*pow(fluid_vel_norm,2.0) + 0.0001;
        Point<dim> velocity_unit_vect = (gradient+Vinf)/fluid_vel_norm;
        
        const double delta = cell->diameter()/sqrt(2);		     
        for (unsigned int i=0; i<vector_dofs_per_cell; ++i)
            {
            double supg_shape_fun = 1.0*(delta)*vector_fe_v.shape_grad(i, q)*velocity_unit_vect;
	    comp_i = comp_dom.vector_fe.system_to_component_index(i).first;
            for (unsigned int j=0; j<vector_dofs_per_cell; ++j)
                {
	        comp_j = comp_dom.vector_fe.system_to_component_index(j).first;
                if (comp_i == comp_j) 
	           local_vector_matrix(i,j) += ( (vector_fe_v.shape_value(i, q)+supg_shape_fun)*
                                                   vector_fe_v.shape_value(j, q) ) *
                                                 vector_fe_v.JxW(q);
	        }
        
	    if ( (cell->material_id() == comp_dom.free_sur_ID1 ||
	          cell->material_id() == comp_dom.free_sur_ID2 ||
		  cell->material_id() == comp_dom.free_sur_ID3 )   )
	       {
               Point<dim> u(quad_nodes_velocities[q](0),quad_nodes_velocities[q](1),quad_nodes_velocities[q](2)); 
               double dEta_dt = gradient[dim-1] + (u-gradient-Vinf)*eta_grad-damping;
               //cout<<q<<"  "<<i<<"   "<<gradient[dim-1]<<"    "<<eta_grad<<"   "<<u-gradient-Vinf<<endl;
               //cout<<q<<"  "<<i<<"   "<<dEta_dt<<endl;
               Point<dim> uu = gradient; //(gradient[dim-1],(u-gradient-Vinf)*eta_grad,Vinf(0));

	       // nodes displacement velocity: x is 0 for now, z is devised to follow the wave, y to make
	       // u parallel to the boat surface
               u(2) = dEta_dt;				     
	       local_vector_rhs(i) += (vector_fe_v.shape_value(i, q)+supg_shape_fun) *
                                        u(comp_i) * vector_fe_v.JxW(q);
               local_vector_rhs_2(i) += (vector_fe_v.shape_value(i, q)+supg_shape_fun) *
                                         uu(comp_i) * vector_fe_v.JxW(q);	 
               }
	    else
	       {
	       local_vector_rhs(i) += (vector_fe_v.shape_value(i, q)+supg_shape_fun) *
                                        0 * vector_fe_v.JxW(q);
               local_vector_rhs_2(i) += (vector_fe_v.shape_value(i, q)+supg_shape_fun) *
                                          0 * vector_fe_v.JxW(q);
               }
            } 
       }

       vector_cell->get_dof_indices (vector_local_dof_indices);
       vector_constr.distribute_local_to_global
       (local_vector_matrix,
	local_vector_rhs,
	vector_local_dof_indices,
	vector_sys_matrix,
	vector_sys_rhs);

       vector_constr.distribute_local_to_global
       (local_vector_rhs_2,
	vector_local_dof_indices,
	vector_sys_rhs_2);
     
       for (unsigned int q=0; q<DphiDt_n_q_points; ++q)
        {
        Point<dim> gradient = DphiDt_node_normals[q]*DphiDt_phi_norm_grads[q] + DphiDt_phi_surf_grads[q];
	Point<dim> Vinf;
	for (unsigned int i = 0; i < dim; i++)
	    Vinf(i) = DphiDt_wind_values[q](i);

	Point<dim> eta_grad = Point<dim>(0,0,0);
        eta_grad = eta_grad + vector_eta_surf_grads[q];
	eta_grad(dim-1) = 0;

	Point<dim> phi_surf_grad_corrected = Point<dim>(0,0,0);
        phi_surf_grad_corrected(0) = DphiDt_phi_surf_grads[q][0] -
                                     DphiDt_phi_surf_grads[q][2]*DphiDt_node_normals[q][0]/DphiDt_node_normals[q][2];
        phi_surf_grad_corrected(1) = DphiDt_phi_surf_grads[q][1] -
                                     DphiDt_phi_surf_grads[q][2]*DphiDt_node_normals[q][1]/DphiDt_node_normals[q][2];
	//phi_surf_grad_corrected(dim-1) = 0;
        //phi_surf_grad_corrected(0) = gradient(0);
        //phi_surf_grad_corrected(1) = gradient(1);
        //phi_surf_grad_corrected(2) = gradient(2);

        double fluid_vel_norm = (gradient+Vinf).norm();
        if (fluid_vel_norm < 1e-3)
            fluid_vel_norm = -8.0e+05*pow(fluid_vel_norm,3.0) + 1.7e+03*pow(fluid_vel_norm,2.0) + 0.0001;
        Point<dim> velocity_unit_vect = (gradient+Vinf)/fluid_vel_norm;

        Point<dim> eta_grad_unit_vect = Point<3>(0,0,0);
        if (sqrt((eta_grad).square()) > 1e-3)
           eta_grad_unit_vect = (eta_grad)/sqrt((eta_grad).square());
        double Fr = sqrt((gradient+Vinf).square())/sqrt(g*comp_dom.Lx_boat);
        double eta = DphiDt_node_positions[q](dim-1);
        double delta_p = 0;
        if ( (eta_grad*(gradient+Vinf) > 0) &&
             (sqrt(eta*eta+eta_grad.square()*pow(Fr,4.0)) > 0.1725*Fr*Fr) &&
             (cell->material_id() == comp_dom.free_sur_ID1 ||
	          cell->material_id() == comp_dom.free_sur_ID2 ||
		  cell->material_id() == comp_dom.free_sur_ID3 ) )
           {
           delta_p = 0.0*eta_grad.square()/Fr/Fr*(velocity_unit_vect*eta_grad_unit_vect)*
                     (eta+0.117*Fr*Fr);
           //std::cout<<"eta "<<DphiDt_node_positions[q]<<"    |grad_eta|^2 "<<eta_grad.square()<<"    Fr "<<Fr<<std::endl;
           //sstd::cout<<sqrt(eta*eta+eta_grad.square()*pow(Fr,4.0))<<" vs "<<0.1725*Fr*Fr<<std::endl;
           }
        
	//double nu_0 = 1.0;
        double nu_0 = 0;
	double damping = nu_0 *
	                 pow((fmax(abs(vector_node_positions[q](0)),20)-20)/(31.169752-20),2) *
			 DphiDt_phi_norm_grads[q];

//        double damping = nu_0 *
//                         pow((fmax(abs(DphiDt_node_positions[q](0)),25)-25)/(25),2) *
//                         DphiDt_phi_norm_grads[q];
		 
        const double delta = cell->diameter()/sqrt(2);	
	for (unsigned int i=0; i<DphiDt_dofs_per_cell; ++i)
            {
            double supg_shape_fun = 1.0*(delta)*fe_v.shape_grad(i, q)*velocity_unit_vect;
            for (unsigned int j=0; j<DphiDt_dofs_per_cell; ++j)
              local_DphiDt_matrix(i,j) += ((fe_v.shape_value (i, q)+supg_shape_fun) *
                                           fe_v.shape_value (j, q)) *
                                           fe_v.JxW(q);
              if (cell->material_id() == comp_dom.free_sur_ID1 ||
	          cell->material_id() == comp_dom.free_sur_ID2 ||
		  cell->material_id() == comp_dom.free_sur_ID3 )
                 {
                 Point<dim> u(quad_nodes_velocities[q](0),quad_nodes_velocities[q](1),quad_nodes_velocities[q](2)); 
                 double dEta_dt = (gradient[dim-1]+(u - gradient - Vinf)*eta_grad-damping);
	
	         // nodes displacement velocity: x is 0 for now, z is devised to follow the wave, y to make
	         // u parallel to the boat surface
	         u(2) = dEta_dt;	
                 local_DphiDt_rhs(i) +=   (fe_v.shape_value (i, q)+supg_shape_fun) *
                                          ((gradient * gradient)/2 - g*DphiDt_node_positions[q](dim-1) +
				           phi_surf_grad_corrected*(u - gradient - Vinf) -
                                           damping - delta_p) *
                                          fe_v.JxW(q);
                 local_DphiDt_rhs_2(i) +=   (fe_v.shape_value (i, q)+supg_shape_fun) *
                                            (delta_p) *
                                            fe_v.JxW(q);
                 //cout<<q<<"  "<<i<<"   "<<fe_v.shape_grad(i, q)<<"   "<<velocity_unit_vect<<"   "<<1.0*(delta)<<endl;
                 //cout<<q<<"  "<<i<<"   "<<(fe_v.shape_value (i, q)+supg_shape_fun)<<"   "<<((gradient * gradient)/2 - g*DphiDt_node_positions[q](dim-1) +
		//		           phi_surf_grad_corrected*(u - gradient - Vinf) -
                  //                         damping - delta_p)<<"   "<<fe_v.JxW(q)<<endl;
                 }
	      else
	         local_DphiDt_rhs(i) += 0;
              			  
          }
        }

     cell->get_dof_indices (DphiDt_local_dof_indices);
     constr.distribute_local_to_global
       (local_DphiDt_matrix,
	local_DphiDt_rhs,
	DphiDt_local_dof_indices,
	DphiDt_sys_matrix,
	DphiDt_sys_rhs);

     constr.distribute_local_to_global
       (local_DphiDt_rhs_2,
	DphiDt_local_dof_indices,
	DphiDt_sys_rhs_2);
   }


   rhs_evaluations_counter++;
 
   
/*   //checking the volume conservation	
   double volume = 0; 
   
   cell = comp_dom.dh.begin_active();
   for (; cell!=endc; ++cell)
        {
	 fe_v.reinit(cell);
         if (cell->material_id() == comp_dom.free_sur_ID1 ||
	     cell->material_id() == comp_dom.free_sur_ID2 ||
             cell->material_id() == comp_dom.free_sur_ID3 )
	    {
	    const std::vector<Point<dim> > &q_points = fe_v.get_quadrature_points();
	    const std::vector<Point<dim> > &normals = fe_v.get_normal_vectors();
	    
	    for (unsigned int q=0; q<DphiDt_n_q_points; ++q){
	    //std::cout<<q_points[q](0)<<" "<<q_points[q](1)<<std::endl;
	              volume += q_points[q](dim-1) * normals[q](dim-1) * fe_v.JxW(q);}
	     }
	 }
	 
    std::cout<<"Volume: "<<volume<<std::endl;//*/
	 

}   


template <int dim>
void FreeSurface<dim>::compute_potential_gradients(Vector<double> &complete_potential_gradients,
                                                   const Vector<double> & phi,
                                                   const Vector<double> & dphi_dn)
{

   std::vector<unsigned int> vector_dofs(comp_dom.vector_fe.dofs_per_cell); 
   std::vector<Point<dim> > vector_support_points(comp_dom.vector_dh.n_dofs());
   DoFTools::map_dofs_to_support_points<dim-1, dim>( *comp_dom.mapping, comp_dom.vector_dh, vector_support_points);

   vector_sys_matrix.reinit (vector_sparsity_pattern);
   Vector<double> complete_gradient_sys_rhs(comp_dom.vector_dh.n_dofs());


  FEValues<dim-1,dim> vector_fe_v(*comp_dom.mapping, comp_dom.vector_fe, *comp_dom.quadrature,
			             update_values | update_gradients | 
			             update_cell_normal_vectors |
			             update_quadrature_points |
			             update_JxW_values);

   FEValues<dim-1,dim> fe_v(*comp_dom.mapping, comp_dom.fe, *comp_dom.quadrature,
			    update_values | update_gradients | 
			    update_cell_normal_vectors |
			    update_quadrature_points |
			    update_JxW_values);

   
   const unsigned int vector_n_q_points = vector_fe_v.n_quadrature_points;
   const unsigned int   vector_dofs_per_cell   = comp_dom.vector_fe.dofs_per_cell;
   std::vector<unsigned int> vector_local_dof_indices (vector_dofs_per_cell);

   FullMatrix<double>   local_vector_matrix (vector_dofs_per_cell, vector_dofs_per_cell);
   Vector<double>       local_complete_gradient_rhs (vector_dofs_per_cell);

   std::vector< Tensor<1,dim> > vector_phi_surf_grads(vector_n_q_points);
   std::vector<double> vector_phi_norm_grads(vector_n_q_points);


   cell_it
   vector_cell = comp_dom.vector_dh.begin_active(),
   vector_endc = comp_dom.vector_dh.end();
   cell_it
   cell = comp_dom.dh.begin_active(),
   endc = comp_dom.dh.end();

   for (; cell!=endc,vector_cell!=vector_endc; ++cell,++vector_cell)
    {
     Assert(cell->index() == vector_cell->index(), ExcInternalError());
        
     vector_fe_v.reinit (vector_cell);
     fe_v.reinit(cell);     

     local_vector_matrix = 0;
     local_complete_gradient_rhs = 0;
     
     fe_v.get_function_gradients((const Vector<double>&)phi, vector_phi_surf_grads);
     fe_v.get_function_values(dphi_dn, vector_phi_norm_grads);


     const std::vector<Point<dim> > &vector_node_normals = vector_fe_v.get_normal_vectors();
     const std::vector<Point<dim> > &vector_node_positions = vector_fe_v.get_quadrature_points();
     
     double g = 9.81;
     
     std::vector<Vector<double> > vector_wind_values(vector_n_q_points);
     for (unsigned int q=0; q<vector_n_q_points; ++q)
         {
	 vector_wind_values[q].reinit(dim,true);
	 wind.vector_value(vector_node_positions[q], vector_wind_values[q]);
         }

     unsigned int comp_i, comp_j;
     for (unsigned int q=0; q<vector_n_q_points; ++q)
        {
        Point<dim> gradient = vector_node_normals[q]*vector_phi_norm_grads[q] + vector_phi_surf_grads[q];

        for (unsigned int i=0; i<vector_dofs_per_cell; ++i)
            {
	    comp_i = comp_dom.vector_fe.system_to_component_index(i).first;
            for (unsigned int j=0; j<vector_dofs_per_cell; ++j)
                {
	        comp_j = comp_dom.vector_fe.system_to_component_index(j).first;
                if (comp_i == comp_j) 
	           {
                   local_vector_matrix(i,j) += ( (vector_fe_v.shape_value(i, q))*
                                                   vector_fe_v.shape_value(j, q) ) *
                                                   vector_fe_v.JxW(q);

                   }
	        }			     
	    local_complete_gradient_rhs(i) += (vector_fe_v.shape_value(i, q)) *
                                               gradient(comp_i) * vector_fe_v.JxW(q);
            } 
       }

       vector_cell->get_dof_indices (vector_local_dof_indices);
       vector_constraints.distribute_local_to_global
       (local_vector_matrix,
	local_complete_gradient_rhs,
	vector_local_dof_indices,
	vector_sys_matrix,
	complete_gradient_sys_rhs);

    }
  
   SparseDirectUMFPACK vector_direct;
   vector_direct.initialize(vector_sys_matrix);
   vector_direct.vmult(complete_potential_gradients,complete_gradient_sys_rhs);
   vector_constraints.distribute(complete_potential_gradients);

   for (unsigned int i=0; i<comp_dom.vector_dh.n_dofs(); ++i)
       {
       std::set <unsigned int> doubles = comp_dom.vector_double_nodes_set[i];
       if ( doubles.size() > 1 )
          {
          double average = 0.0;
          for (std::set<unsigned int>::iterator pos=doubles.begin(); pos != doubles.end(); ++pos)
              {
              average += complete_potential_gradients(*pos);
              }
          average /= doubles.size();
          for (std::set<unsigned int>::iterator pos=doubles.begin(); pos != doubles.end(); ++pos)
              {
              complete_potential_gradients(*pos) = average;
              }
          }
       }

   vector_constraints.distribute(complete_potential_gradients);
}


 



template <int dim>
void FreeSurface<dim>::compute_pressure(Vector<double> & press, 
                                        Vector<double> & comp_1, Vector<double> & comp_2, Vector<double> & comp_3, Vector<double> & comp_4, 
                                        const double t,
                                        const Vector<double> & solution,
                                        const Vector<double> & solution_dot)
{

   std::cout<<"Computing pressure... "<<std::endl;

   double g = 9.81;
   double rho = 1025.1;

   comp_dom.update_support_points();
   press.reinit(comp_dom.dh.n_dofs());


   comp_1.reinit(comp_dom.dh.n_dofs());
   comp_2.reinit(comp_dom.dh.n_dofs());
   comp_3.reinit(comp_dom.dh.n_dofs());
   comp_4.reinit(comp_dom.dh.n_dofs());
   Vector<double> rhs_comp_1(comp_dom.dh.n_dofs());
   Vector<double> rhs_comp_2(comp_dom.dh.n_dofs());
   Vector<double> rhs_comp_3(comp_dom.dh.n_dofs());
   Vector<double> rhs_comp_4(comp_dom.dh.n_dofs());
   rhs_comp_1 = 0;
   rhs_comp_2 = 0;
   rhs_comp_3 = 0;
   rhs_comp_4 = 0;

   DphiDt_sys_rhs.reinit (comp_dom.dh.n_dofs());
   DphiDt_sys_matrix.reinit(DphiDt_sparsity_pattern); 

   VectorView<double> phi(comp_dom.dh.n_dofs(),solution.begin()+comp_dom.vector_dh.n_dofs());
   VectorView<double> dphi_dn(comp_dom.dh.n_dofs(),solution.begin()+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs());
   VectorView<double> DphiDt(comp_dom.dh.n_dofs(),solution_dot.begin()+comp_dom.vector_dh.n_dofs());
   VectorView<double> node_vels(comp_dom.vector_dh.n_dofs(),solution_dot.begin());
   AssertDimension(DphiDt.size(), node_vels.size()/dim);

   Point<3> hull_baricenter_displ;
   Point<3> hull_baricenter_vel;
   Point<3> hull_ang_vel;
   for (unsigned int d=0; d<3; ++d)
       {
       hull_baricenter_vel(d) = solution(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+d);
       hull_baricenter_displ(d) = solution(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+3+d);
       hull_ang_vel(d) = solution(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+6+d);
       }

   Point<3> baricenter_pos(hull_baricenter_displ(0)+comp_dom.boat_model.reference_hull_baricenter(0),
                           hull_baricenter_displ(1)+comp_dom.boat_model.reference_hull_baricenter(1),
                           hull_baricenter_displ(2)+comp_dom.boat_model.reference_hull_baricenter(2));

   Vector<double> v_x(comp_dom.dh.n_dofs());
   Vector<double> v_y(comp_dom.dh.n_dofs());
   Vector<double> v_z(comp_dom.dh.n_dofs());

   for (unsigned int i = 0; i < comp_dom.dh.n_dofs(); i++)
       {
       v_x(i) = node_vels(i*dim)+hull_ang_vel(1)*(comp_dom.support_points[i](2)-baricenter_pos(2))-
                                 hull_ang_vel(2)*(comp_dom.support_points[i](1)-baricenter_pos(1))+hull_baricenter_vel(0);
       v_y(i) = node_vels(i*dim+1)+hull_ang_vel(2)*(comp_dom.support_points[i](0)-baricenter_pos(0))-
                                 hull_ang_vel(0)*(comp_dom.support_points[i](2)-baricenter_pos(2))+hull_baricenter_vel(1);
       v_z(i) = node_vels(i*dim+2)+hull_ang_vel(0)*(comp_dom.support_points[i](1)-baricenter_pos(1))-
                                 hull_ang_vel(1)*(comp_dom.support_points[i](0)-baricenter_pos(0))+hull_baricenter_vel(2);
       }

   FEValues<dim-1,dim> fe_v(*comp_dom.mapping, comp_dom.fe, *comp_dom.quadrature,
			    update_values | update_gradients |
			    update_cell_normal_vectors |
			    update_quadrature_points |
			    update_JxW_values);

   const unsigned int DphiDt_n_q_points = fe_v.n_quadrature_points;
   const unsigned int  DphiDt_dofs_per_cell   = comp_dom.fe.dofs_per_cell;
   std::vector<unsigned int> DphiDt_local_dof_indices (DphiDt_dofs_per_cell);

   FullMatrix<double>   local_DphiDt_matrix (DphiDt_dofs_per_cell, DphiDt_dofs_per_cell);
   Vector<double>       local_DphiDt_rhs (DphiDt_dofs_per_cell);

   Vector<double>       local_DphiDt_rhs_comp_1 (DphiDt_dofs_per_cell);
   Vector<double>       local_DphiDt_rhs_comp_2 (DphiDt_dofs_per_cell);
   Vector<double>       local_DphiDt_rhs_comp_3 (DphiDt_dofs_per_cell);
   Vector<double>       local_DphiDt_rhs_comp_4 (DphiDt_dofs_per_cell);
      
   std::vector< Tensor<1,dim> > DphiDt_phi_surf_grads(DphiDt_n_q_points);
   std::vector<double> DphiDt_phi_norm_grads(DphiDt_n_q_points);
   
   std::vector<double> DphiDt_v_x_values(DphiDt_n_q_points);
   std::vector<double> DphiDt_v_y_values(DphiDt_n_q_points);
   std::vector<double> DphiDt_v_z_values(DphiDt_n_q_points);
   std::vector<double> DphiDt_DphiDt_values(DphiDt_n_q_points);

   cell_it
   cell = comp_dom.dh.begin_active(),
   endc = comp_dom.dh.end();
   Point<dim> press_force_test_1;
   Point<dim> press_force_test_2;
   Point<dim> press_moment;
   
   for (; cell!=endc; ++cell)
     {
     fe_v.reinit(cell);
     local_DphiDt_rhs = 0;
     local_DphiDt_matrix = 0;

     local_DphiDt_rhs_comp_1 = 0;
     local_DphiDt_rhs_comp_2 = 0;
     local_DphiDt_rhs_comp_3 = 0;
     local_DphiDt_rhs_comp_4 = 0;

     fe_v.get_function_gradients((const Vector<double>&)phi, DphiDt_phi_surf_grads);
     fe_v.get_function_values((const Vector<double>&)dphi_dn, DphiDt_phi_norm_grads);
     fe_v.get_function_values(v_x, DphiDt_v_x_values);
     fe_v.get_function_values(v_y, DphiDt_v_y_values);
     fe_v.get_function_values(v_z, DphiDt_v_z_values);
     fe_v.get_function_values((const Vector<double>&)DphiDt, DphiDt_DphiDt_values);
     const std::vector<Point<dim> > &DphiDt_node_normals = fe_v.get_normal_vectors();
     const std::vector<Point<dim> > &DphiDt_node_positions = fe_v.get_quadrature_points();
     
     std::vector<Vector<double> > DphiDt_wind_values(DphiDt_n_q_points);
     for (unsigned int q=0; q<DphiDt_n_q_points; ++q)
        {
	DphiDt_wind_values[q].reinit(dim,true);
	wind.vector_value(DphiDt_node_positions[q],DphiDt_wind_values[q]);
	}

     
     for (unsigned int q=0; q<DphiDt_n_q_points; ++q)
        {
        Point<dim> gradient = DphiDt_node_normals[q]*DphiDt_phi_norm_grads[q] + DphiDt_phi_surf_grads[q];
	Point<dim> Vinf;
	for (unsigned int i = 0; i < dim; i++)
	    Vinf(i) = DphiDt_wind_values[q](i);

	Point<dim> phi_surf_grad = Point<dim>(0,0,0);
        phi_surf_grad = phi_surf_grad + DphiDt_phi_surf_grads[q];
	phi_surf_grad(dim-1) = 0;

        Point<dim> node_vel_vect;
        node_vel_vect(0) = DphiDt_v_x_values[q];
        node_vel_vect(1) = DphiDt_v_y_values[q];
        node_vel_vect(2) = DphiDt_v_z_values[q];

        Point<dim> velocity_unit_vect;
        if (sqrt((gradient+Vinf).square()) > 1e-3)
           velocity_unit_vect = (gradient+Vinf)/sqrt((gradient+Vinf).square());

	double press = rho*((Vinf)*(Vinf))/2 -
                       rho*((gradient+Vinf)*(gradient+Vinf))/2 -
                       rho*g*DphiDt_node_positions[q](dim-1) -
                       rho*(DphiDt_DphiDt_values[q]-node_vel_vect*gradient);
        double press2 = -rho*((Vinf)*(gradient)) 
                        -rho*((gradient)*(gradient))/2 
                        -rho*g*DphiDt_node_positions[q](dim-1) 
                        -rho*(DphiDt_DphiDt_values[q]-node_vel_vect*gradient);	
        /*std::cout<<"cell "<<cell<<"  node "<<q<<std::endl;
        std::cout<<"press "<<press<<std::endl;
        std::cout<<"term1 "<<rho*((Vinf)*(Vinf))/2<<std::endl;
        std::cout<<"term2 "<<-rho*((gradient+Vinf)*(gradient+Vinf))/2<<std::endl;
        std::cout<<"term3 "<<-rho*g*DphiDt_node_positions[q](dim-1)<<std::endl;
        std::cout<<"term4 "<<-rho*(DphiDt_DphiDt_values[q]-node_vel_vect*phi_surf_grad)/2<<std::endl;
        std::cout<<std::endl;//*/

	for (unsigned int i=0; i<DphiDt_dofs_per_cell; ++i)
            {
	    for (unsigned int j=0; j<DphiDt_dofs_per_cell; ++j)
              local_DphiDt_matrix(i,j) += ((fe_v.shape_value (i, q)) *
                                           fe_v.shape_value (j, q)) *
                                           fe_v.JxW(q);
            local_DphiDt_rhs(i) +=   (fe_v.shape_value (i, q)) *
                                     (press) *
                                     fe_v.JxW(q);
            local_DphiDt_rhs_comp_1(i) +=   (fe_v.shape_value (i, q)) *
                                     (rho*((Vinf)*(Vinf))/2) *
                                     fe_v.JxW(q);
            local_DphiDt_rhs_comp_2(i) +=   (fe_v.shape_value (i, q)) *
                                     (-rho*((gradient+Vinf)*(gradient+Vinf))/2) *
                                     fe_v.JxW(q);
            local_DphiDt_rhs_comp_3(i) +=   (fe_v.shape_value (i, q)) *
                                     (-rho*g*DphiDt_node_positions[q](dim-1)) *
                                     fe_v.JxW(q);
            local_DphiDt_rhs_comp_4(i) +=   (fe_v.shape_value (i, q)) *
                                     (-rho*(DphiDt_DphiDt_values[q]-node_vel_vect*phi_surf_grad)) *
                                     fe_v.JxW(q);
            
            }
       if ((cell->material_id() == comp_dom.wall_sur_ID1 ||
            cell->material_id() == comp_dom.wall_sur_ID2 ||
            cell->material_id() == comp_dom.wall_sur_ID3 ))
          {
          press_force_test_1 += press*DphiDt_node_normals[q]*fe_v.JxW(q);
          press_force_test_2 += press2*DphiDt_node_normals[q]*fe_v.JxW(q);
          press_moment += press*Point<3>((DphiDt_node_positions[q](1)-baricenter_pos(1))*DphiDt_node_normals[q](2)-
                                         (DphiDt_node_positions[q](2)-baricenter_pos(2))*DphiDt_node_normals[q](1),
                                         (DphiDt_node_positions[q](2)-baricenter_pos(2))*DphiDt_node_normals[q](0)-
                                         (DphiDt_node_positions[q](0)-baricenter_pos(0))*DphiDt_node_normals[q](2),
                                         (DphiDt_node_positions[q](0)-baricenter_pos(0))*DphiDt_node_normals[q](1)-
                                         (DphiDt_node_positions[q](1)-baricenter_pos(1))*DphiDt_node_normals[q](0))*fe_v.JxW(q);
          }
       }

    cell->get_dof_indices (DphiDt_local_dof_indices);
    constraints.distribute_local_to_global
       (local_DphiDt_matrix,
	local_DphiDt_rhs,
	DphiDt_local_dof_indices,
	DphiDt_sys_matrix,
	DphiDt_sys_rhs);

     constraints.distribute_local_to_global
       (local_DphiDt_rhs_comp_1,
	DphiDt_local_dof_indices,
	rhs_comp_1);
     constraints.distribute_local_to_global
       (local_DphiDt_rhs_comp_2,
	DphiDt_local_dof_indices,
	rhs_comp_2);
     constraints.distribute_local_to_global
       (local_DphiDt_rhs_comp_3,
	DphiDt_local_dof_indices,
	rhs_comp_3);
     constraints.distribute_local_to_global
       (local_DphiDt_rhs_comp_4,
	DphiDt_local_dof_indices,
	rhs_comp_4);

     

    }

   //DphiDt_sys_matrix.print(std::cout);
   SparseDirectUMFPACK pressure_direct;
   pressure_direct.initialize(DphiDt_sys_matrix);
   pressure_direct.vmult(press, DphiDt_sys_rhs);
   constraints.distribute(press);

   pressure_direct.vmult(comp_1, rhs_comp_1);
   constraints.distribute(comp_1);
   pressure_direct.vmult(comp_2, rhs_comp_2);
   constraints.distribute(comp_2);
   pressure_direct.vmult(comp_3, rhs_comp_3);
   constraints.distribute(comp_3);
   pressure_direct.vmult(comp_4, rhs_comp_4);
   constraints.distribute(comp_4);


   Vector<double> complete_potential_gradients(comp_dom.vector_dh.n_dofs()); 
   compute_potential_gradients(complete_potential_gradients,phi,dphi_dn);
 

   wind.set_time(t);
   Vector<double> instantVelValue(dim);
   Point<dim> zzero(0,0,0);
   wind.vector_value(zzero,instantVelValue);
   Point<3> V_inf(instantVelValue(0),
                  instantVelValue(1),
                  instantVelValue(2));

   for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
       {
       Point<3> gradient(complete_potential_gradients(3*i),
                         complete_potential_gradients(3*i+1),
                         complete_potential_gradients(3*i+2));
       Point<3> v(node_vels(3*i),node_vels(3*i+1),node_vels(3*i+2));
       press(i) = 0.5*rho*V_inf*V_inf -
                     0.5*rho*(V_inf+gradient)*(V_inf+gradient) -
                     rho*g*comp_dom.support_points[i](2) -
                     rho*(DphiDt(i)-v*gradient);
       }

   // preparing iges normal vectors vector for pressure computation
   Vector<double> iges_normals_x_values(comp_dom.dh.n_dofs());
   Vector<double> iges_normals_y_values(comp_dom.dh.n_dofs());
   Vector<double> iges_normals_z_values(comp_dom.dh.n_dofs());
   for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
       {
       iges_normals_x_values(i) = comp_dom.iges_normals[i](0);
       iges_normals_y_values(i) = comp_dom.iges_normals[i](1);
       iges_normals_z_values(i) = comp_dom.iges_normals[i](2);
       }


   // pressure force computation
   Point<dim> press_force;
   double wet_surface = 0;

   cell = comp_dom.dh.begin_active(),
   endc = comp_dom.dh.end();

   std::vector<double> pressure_quad_values(DphiDt_n_q_points);
   std::vector<double> n_x_quad_values(DphiDt_n_q_points);
   std::vector<double> n_y_quad_values(DphiDt_n_q_points);
   std::vector<double> n_z_quad_values(DphiDt_n_q_points);

   for (; cell!=endc; ++cell)
       {
       if ((cell->material_id() == comp_dom.wall_sur_ID1 ||
            cell->material_id() == comp_dom.wall_sur_ID2 ||
            cell->material_id() == comp_dom.wall_sur_ID3 ))
          {
          fe_v.reinit(cell);
          fe_v.get_function_values(press, pressure_quad_values);
          fe_v.get_function_values(iges_normals_x_values, n_x_quad_values);
          fe_v.get_function_values(iges_normals_y_values, n_y_quad_values);
          fe_v.get_function_values(iges_normals_z_values, n_z_quad_values);

          const std::vector<Point<dim> > &DphiDt_node_normals = fe_v.get_normal_vectors();
          for (unsigned int q=0; q<DphiDt_n_q_points; ++q)
            {
            Point<3> normal(n_x_quad_values[q],n_y_quad_values[q],n_z_quad_values[q]);
            //Point<dim> local_press_force = pressure_quad_values[q]*normal;
            Point<dim> local_press_force = pressure_quad_values[q]*DphiDt_node_normals[q];
            press_force += (local_press_force) * fe_v.JxW(q);
            wet_surface += 1.0 * fe_v.JxW(q);
            }
          }
       }

   // breaking wave additional pressure computation
   Point<3> break_press_force(0.0,0.0,0.0);
   for (tria_it elem=comp_dom.tria.begin_active(); elem!= comp_dom.tria.end();++elem)
       {
       if ((elem->material_id() == comp_dom.free_sur_ID1 ||
            elem->material_id() == comp_dom.free_sur_ID2 ||
            elem->material_id() == comp_dom.free_sur_ID3 )) 
	  {
	  if (elem->at_boundary())
             {
             for (unsigned int f=0; f<GeometryInfo<2>::faces_per_cell; ++f)
                 {
		 if ( elem->face(f)->boundary_indicator() == 27 ||
                      elem->face(f)->boundary_indicator() == 29 ) // left side
                    {
                    std::vector<Point<3> > vertices(4);
                    std::vector<CellData<2> > cells(1);
                    Vector<double> pressure_vect(4);

                    unsigned int index_0 = comp_dom.find_point_id(elem->face(f)->vertex(0),comp_dom.ref_points);
                    std::set<unsigned int> duplicates = comp_dom.double_nodes_set[index_0]; 
                    for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
                        {
                        //cout<<i<<" mpid"<<*pos<<"  is in?"<<boundary_dofs[i][3*(*pos)]<<endl;
                        if (comp_dom.flags[*pos] & water)
	                   {
                           index_0 = *pos;
                           break;
                           }
                        }
                    vertices[0] = comp_dom.support_points[index_0];
                    pressure_vect(0) = fmax(break_wave_press(index_0),1e-4)*rho;

                    unsigned int index_1 = comp_dom.find_point_id(elem->face(f)->vertex(1),comp_dom.ref_points);
                    duplicates = comp_dom.double_nodes_set[index_1]; 
                    for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
                        {
                        //cout<<i<<" mpid"<<*pos<<"  is in?"<<boundary_dofs[i][3*(*pos)]<<endl;
                        if (comp_dom.flags[*pos] & water)
	                   {
                           index_1 = *pos;
                           break;
                           }
                        }
                    vertices[1] = comp_dom.support_points[index_1];
                    pressure_vect(1) = fmax(break_wave_press(index_1),1e-4)*rho;
                    //cout<<"********* "<<pressure_vect(0)<<" "<<pressure_vect(1)<<endl;
                    Point<3> temp_point(vertices[0](0),vertices[0](1),vertices[0](2)+pressure_vect(0)/g/rho);
                    // if point we're dealing with is an internal one, we're projecting it on boat in y direction 
                    if ( (comp_dom.moving_point_ids[3] != index_0) &&
                         (comp_dom.moving_point_ids[4] != index_0) &&
                         (comp_dom.moving_point_ids[5] != index_0) &&
                         (comp_dom.moving_point_ids[6] != index_0) )
                         {
                         //cout<<index_0<<" "<<comp_dom.moving_point_ids[3]<<" "<<comp_dom.moving_point_ids[4]<<" ";
                         //cout<<comp_dom.moving_point_ids[5]<<" "<<comp_dom.moving_point_ids[6]<<endl;
                         comp_dom.boat_model.boat_water_line_left->axis_projection(vertices[2],temp_point);  // y axis projection
                         }
                    //otherwise it needs to be properly projected on the curve (to be implemented yet)
                    else
                         vertices[2] = vertices[0];
                    temp_point = Point<3>(vertices[1](0),vertices[1](1),vertices[1](2)+pressure_vect(1)/g/rho);
                     
                    if ( (comp_dom.moving_point_ids[3] != index_1) &&
                         (comp_dom.moving_point_ids[4] != index_1) &&
                         (comp_dom.moving_point_ids[5] != index_1) &&
                         (comp_dom.moving_point_ids[6] != index_1) )
                         {
                         //cout<<index_0<<" "<<comp_dom.moving_point_ids[3]<<" "<<comp_dom.moving_point_ids[4]<<" ";
                         //cout<<comp_dom.moving_point_ids[5]<<" "<<comp_dom.moving_point_ids[6]<<endl;
                         comp_dom.boat_model.boat_water_line_left->axis_projection(vertices[3],temp_point);  // y axis projection
                         }
                    //otherwise it needs to be properly projected on the curve (to be implemented yet)
                    else
                         vertices[3] = vertices[1];
                    //cout<<"########## "<<vertices[0]<<" "<<vertices[1]<<endl;
                    //cout<<"########## "<<vertices[2]<<" "<<vertices[3]<<endl;
                    pressure_vect(2) = 0.0;
                    pressure_vect(3) = 0.0;
                    if (vertices[1](0) < vertices[0](0))
                       {
                       cells[0].vertices[0]=0;
                       cells[0].vertices[1]=1;
                       cells[0].vertices[2]=3;
                       cells[0].vertices[3]=2;
                       //cout<<"l* "<<elem<<endl;
                       }
                    else
                       {
                       cells[0].vertices[0]=1;
                       cells[0].vertices[1]=0;
                       cells[0].vertices[2]=2;
                       cells[0].vertices[3]=3;
                       //cout<<"l** "<<elem<<endl;
                       }
                    SubCellData subcelldata;
                    Triangulation<dim-1, dim> break_tria;
                    GridTools::delete_unused_vertices (vertices, cells, subcelldata);
                    GridReordering<2,3>::reorder_cells (cells);
                    break_tria.create_triangulation_compatibility(vertices, cells, subcelldata );
                    FE_Q<dim-1,dim> break_fe(1);
                    DoFHandler<dim-1,dim> break_dh(break_tria);
                    break_dh.distribute_dofs(break_fe);

                    FEValues<dim-1,dim> break_fe_v(break_fe, *comp_dom.quadrature,
		              	                   update_values | update_gradients |
		              	                   update_cell_normal_vectors |
		              	                   update_quadrature_points |
		              	                   update_JxW_values);

                    const unsigned int break_n_q_points = break_fe_v.n_quadrature_points;
                    std::vector<double> break_pressure_quad_values(break_n_q_points);
                    cell_it cell = break_dh.begin_active();
                    break_fe_v.reinit(cell);
                    break_fe_v.get_function_values(pressure_vect, break_pressure_quad_values);
                    const std::vector<Point<dim> > &break_node_normals = break_fe_v.get_normal_vectors();
                    const std::vector<Point<dim> > &break_quad_nodes = break_fe_v.get_quadrature_points();

                    for (unsigned int q=0; q<break_n_q_points; ++q)
                        {
                        Point<dim> local_press_force = break_pressure_quad_values[q]*break_node_normals[q];
                        break_press_force += (local_press_force) * break_fe_v.JxW(q);
                        //cout<<q<<"  F = ("<<local_press_force<<")    n = ("<<break_node_normals[q]<<")"<<"  "<<break_fe_v.JxW(q);
                        //cout<<"  p = ("<<break_quad_nodes[q]<<")"<<endl;
                        }


                    }
                 if ( elem->face(f)->boundary_indicator() == 26 ||
                      elem->face(f)->boundary_indicator() == 28 ) // right side
                    {
                    std::vector<Point<3> > vertices(4);
                    std::vector<CellData<2> > cells(1);
                    Vector<double> pressure_vect(4);

                    unsigned int index_0 = comp_dom.find_point_id(elem->face(f)->vertex(0),comp_dom.ref_points);
                    std::set<unsigned int> duplicates = comp_dom.double_nodes_set[index_0]; 
                    for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
                        {
                        //cout<<i<<" mpid"<<*pos<<"  is in?"<<boundary_dofs[i][3*(*pos)]<<endl;
                        if (comp_dom.flags[*pos] & water)
	                   {
                           index_0 = *pos;
                           break;
                           }
                        }
                    vertices[0] = comp_dom.support_points[index_0];
                    pressure_vect(0) = fmax(break_wave_press(index_0),1e-4)*rho;

                    unsigned int index_1 = comp_dom.find_point_id(elem->face(f)->vertex(1),comp_dom.ref_points);
                    duplicates = comp_dom.double_nodes_set[index_1]; 
                    for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
                        {
                        //cout<<i<<" mpid"<<*pos<<"  is in?"<<boundary_dofs[i][3*(*pos)]<<endl;
                        if (comp_dom.flags[*pos] & water)
	                   {
                           index_1 = *pos;
                           break;
                           }
                        }
                    vertices[1] = comp_dom.support_points[index_1];
                    pressure_vect(1) = fmax(break_wave_press(index_1),1e-4)*rho;
                    Point<3> temp_point(vertices[0](0),vertices[0](1),vertices[0](2)+pressure_vect(0)/g/rho);
                    // if point we're dealing with is an internal one, we're projectiong it on boat in y direction 
                    if ( (comp_dom.moving_point_ids[3] != index_0) &&
                         (comp_dom.moving_point_ids[4] != index_0) &&
                         (comp_dom.moving_point_ids[5] != index_0) &&
                         (comp_dom.moving_point_ids[6] != index_0) )
                         {
                         //cout<<index_0<<" "<<comp_dom.moving_point_ids[3]<<" "<<comp_dom.moving_point_ids[4]<<" ";
                         //cout<<comp_dom.moving_point_ids[5]<<" "<<comp_dom.moving_point_ids[6]<<endl;
                         comp_dom.boat_model.boat_water_line_right->axis_projection(vertices[2],temp_point);  // y axis projection
                         }
                    //otherwise it needs to be properly projected on the curve
                    else
                         vertices[2] = vertices[0];
                    temp_point = Point<3>(vertices[1](0),vertices[1](1),vertices[1](2)+pressure_vect(1)/g/rho);
                    if ( (comp_dom.moving_point_ids[3] != index_1) &&
                         (comp_dom.moving_point_ids[4] != index_1) &&
                         (comp_dom.moving_point_ids[5] != index_1) &&
                         (comp_dom.moving_point_ids[6] != index_1) )
                         {
                         //cout<<index_1<<" "<<comp_dom.moving_point_ids[3]<<" "<<comp_dom.moving_point_ids[4]<<" ";
                         //cout<<comp_dom.moving_point_ids[5]<<" "<<comp_dom.moving_point_ids[6]<<endl;
                         comp_dom.boat_model.boat_water_line_right->axis_projection(vertices[3],temp_point);  // y axis projection
                         }
                    //otherwise it needs to be properly projected on the curve
                    else
                         vertices[3] = vertices[1];

                    pressure_vect(2) = 0.0;
                    pressure_vect(3) = 0.0;
                    if (vertices[1](0) > vertices[0](0))
                       {
                       cells[0].vertices[0]=0;
                       cells[0].vertices[1]=1;
                       cells[0].vertices[2]=3;
                       cells[0].vertices[3]=2;
                       //cout<<"* "<<elem<<endl;
                       }
                    else
                       {
                       cells[0].vertices[0]=1;
                       cells[0].vertices[1]=0;
                       cells[0].vertices[2]=2;
                       cells[0].vertices[3]=3;
                       //cout<<"** "<<elem<<endl;
                       }
                    SubCellData subcelldata;
                    Triangulation<dim-1, dim> break_tria;
                    GridTools::delete_unused_vertices (vertices, cells, subcelldata);
                    GridReordering<2,3>::reorder_cells (cells);
                    break_tria.create_triangulation_compatibility(vertices, cells, subcelldata );
                    FE_Q<dim-1,dim> break_fe(1);
                    DoFHandler<dim-1,dim> break_dh(break_tria);
                    break_dh.distribute_dofs(break_fe);
                    FEValues<dim-1,dim> break_fe_v(break_fe, *comp_dom.quadrature,
		              	                   update_values | update_gradients |
		              	                   update_cell_normal_vectors |
		              	                   update_quadrature_points |
		              	                   update_JxW_values);

                    const unsigned int break_n_q_points = break_fe_v.n_quadrature_points;
                    std::vector<double> break_pressure_quad_values(break_n_q_points);
                    cell_it cell = break_dh.begin_active();
                    break_fe_v.reinit(cell);
                    break_fe_v.get_function_values(pressure_vect, break_pressure_quad_values);
                    const std::vector<Point<dim> > &break_node_normals = break_fe_v.get_normal_vectors();
                    const std::vector<Point<dim> > &break_quad_nodes = break_fe_v.get_quadrature_points();

                    for (unsigned int q=0; q<break_n_q_points; ++q)
                        {
                        Point<dim> local_press_force = break_pressure_quad_values[q]*break_node_normals[q];
                        break_press_force += (local_press_force) * break_fe_v.JxW(q);
                        if (local_press_force.distance(Point<3>(0.0,0.0,0.0)) > 1e-3)
                           {
                           //cout<<q<<"  F = ("<<local_press_force<<")    n = ("<<break_node_normals[q]<<")"<<"  "<<break_fe_v.JxW(q);
                           //cout<<"  p = ("<<break_quad_nodes[q]<<")"<<endl;
                           }
                        }


                    }
                }
             }
          }
       }

   // transom pressure force computation
   wind.set_time(t);
   Vector<double> instantWindValue(dim);
   Point<dim> zero(0,0,0);
   wind.vector_value(zero,instantWindValue);
   double Vinf = instantWindValue(0);
   Point<dim> transom_press_force;
   double transom_wet_surface = 0;

   if (!comp_dom.no_boat && comp_dom.boat_model.is_transom)
      {
      //double transom_draft = ref_transom_wet_surface/(fabs(comp_dom.boat_model.PointLeftTransom(1))+
      //                                                     fabs(comp_dom.boat_model.PointRightTransom(1)));
      double transom_draft = fabs(comp_dom.boat_model.CurrentPointCenterTransom(2));
      double transom_aspect_ratio = (fabs(comp_dom.boat_model.CurrentPointLeftTransom(1))+
                                     fabs(comp_dom.boat_model.CurrentPointRightTransom(1)))/transom_draft;
      double FrT = sqrt(Vinf*Vinf)/sqrt(9.81*transom_draft);
      double ReT = sqrt(9.81*pow(transom_draft,3.0))/1.307e-6;
      double eta_dry=fmin(0.05*pow(FrT,2.834)*pow(transom_aspect_ratio,0.1352)*pow(ReT,0.01338),1.0);
      //cout<<"eta_dry "<<eta_dry<<"    mean transom draft "<<mean_transom_draft<<endl;
      cout<<"eta_dry "<<eta_dry<<endl;
      std::vector<Point<3> > vertices;
      std::vector<CellData<2> > cells;
      std::vector<double> transom_pressure_vect;
      for (tria_it elem=comp_dom.tria.begin_active(); elem!= comp_dom.tria.end();++elem)
          {
	  if ((elem->material_id() == comp_dom.wall_sur_ID1 ||
	       elem->material_id() == comp_dom.wall_sur_ID2 ||
	       elem->material_id() == comp_dom.wall_sur_ID3 )) 
	     {
	     if (elem->at_boundary())
                {
                for (unsigned int f=0; f<GeometryInfo<2>::faces_per_cell; ++f)
		    if ( elem->face(f)->boundary_indicator() == 32 ||
                         elem->face(f)->boundary_indicator() == 37 )
                       {
                       unsigned int index_0 = comp_dom.find_point_id(elem->face(f)->vertex(0),comp_dom.ref_points);
                       double pressure_0 = (1-eta_dry)*comp_dom.ref_points[3*index_0](2)*rho*g;
                       unsigned int index_1 = comp_dom.find_point_id(elem->face(f)->vertex(1),comp_dom.ref_points);
                       double pressure_1 = (1-eta_dry)*comp_dom.ref_points[3*index_1](2)*rho*g;
                       if (comp_dom.ref_points[3*index_1](1) < comp_dom.ref_points[3*index_0](1))
                          {
                          vertices.push_back(comp_dom.ref_points[3*index_0]);
                          //cout<<comp_dom.ref_points[3*index_0]<<" / "<<elem->face(f)->vertex(0)<<endl;
                          transom_pressure_vect.push_back(pressure_0);
                          vertices.push_back(comp_dom.ref_points[3*index_1]);
                          //cout<<comp_dom.ref_points[3*index_1]<<" / "<<elem->face(f)->vertex(1)<<endl;
                          transom_pressure_vect.push_back(pressure_1);
                          vertices.push_back(comp_dom.ref_points[3*index_1]+
                                             -1.0*Point<3>(0.0,0.0,(1-eta_dry)*comp_dom.ref_points[3*index_1](2)));
                          //cout<<comp_dom.ref_points[3*index_1]-
                          //                   Point<3>(0.0,0.0,(1-eta_dry)*comp_dom.ref_points[3*index_1](2))<<endl;
                          transom_pressure_vect.push_back(0.0);
                          vertices.push_back(comp_dom.ref_points[3*index_0]+
                                              -1.0*Point<3>(0.0,0.0,(1-eta_dry)*comp_dom.ref_points[3*index_0](2)));
                          //cout<<comp_dom.ref_points[3*index_0]-
                          //                   Point<3>(0.0,0.0,(1-eta_dry)*comp_dom.ref_points[3*index_0](2))<<endl;
                          transom_pressure_vect.push_back(0.0);
                          cells.resize(cells.size()+1);
                          cells[cells.size()-1].vertices[0]=4*(cells.size()-1)+0;
                          cells[cells.size()-1].vertices[1]=4*(cells.size()-1)+1;
                          cells[cells.size()-1].vertices[2]=4*(cells.size()-1)+2;
                          cells[cells.size()-1].vertices[3]=4*(cells.size()-1)+3;
                          }
                       else
                          {
                          vertices.push_back(comp_dom.ref_points[3*index_1]);
                          transom_pressure_vect.push_back(pressure_1);
                          vertices.push_back(comp_dom.ref_points[3*index_0]);
                          transom_pressure_vect.push_back(pressure_0);
                          vertices.push_back(comp_dom.ref_points[3*index_0]+
                                             -1.0*Point<3>(0.0,0.0,(1-eta_dry)*comp_dom.ref_points[3*index_0](2)));
                          transom_pressure_vect.push_back(0.0);
                          vertices.push_back(comp_dom.ref_points[3*index_1]+
                                             -1.0*Point<3>(0.0,0.0,(1-eta_dry)*comp_dom.ref_points[3*index_1](2)));
                          transom_pressure_vect.push_back(0.0);
                          cells.resize(cells.size()+1);
                          cells[cells.size()-1].vertices[0]=4*(cells.size()-1)+0;
                          cells[cells.size()-1].vertices[1]=4*(cells.size()-1)+1;
                          cells[cells.size()-1].vertices[2]=4*(cells.size()-1)+2;
                          cells[cells.size()-1].vertices[3]=4*(cells.size()-1)+3;
                          }                          
                       }
	        }
             }
          }
      Vector<double> transom_pressure(transom_pressure_vect.size());
      for (unsigned int i=0; i<transom_pressure_vect.size(); ++i)
          {
          transom_pressure(i) = transom_pressure_vect[i];
          //cout<<i<<" " <<transom_pressure(i)<<" "<<vertices[i]<<endl;
          }

      SubCellData subcelldata;
      Triangulation<dim-1, dim> transom_tria;
      GridTools::delete_unused_vertices (vertices, cells, subcelldata);
      GridReordering<2,3>::reorder_cells (cells);   
      transom_tria.create_triangulation_compatibility(vertices, cells, subcelldata );
      FE_Q<dim-1,dim> transom_fe(1);
      DoFHandler<dim-1,dim> transom_dh(transom_tria);
      transom_dh.distribute_dofs(transom_fe);

      FEValues<dim-1,dim> transom_fe_v(transom_fe, *comp_dom.quadrature,
			               update_values | update_gradients |
			               update_cell_normal_vectors |
			               update_quadrature_points |
			               update_JxW_values);

      const unsigned int transom_n_q_points = transom_fe_v.n_quadrature_points;
      //const unsigned int  DphiDt_dofs_per_cell   = comp_dom.fe.dofs_per_cell;

      if (eta_dry < 1.0)
         {
         std::vector<double> transom_pressure_quad_values(transom_n_q_points);
         for (cell_it cell = transom_dh.begin_active(); cell!=transom_dh.end(); ++cell)
             {
             transom_fe_v.reinit(cell);
             transom_fe_v.get_function_values(transom_pressure, transom_pressure_quad_values);
             const std::vector<Point<dim> > &transom_node_normals = transom_fe_v.get_normal_vectors();
             const std::vector<Point<dim> > &transom_quad_nodes = transom_fe_v.get_quadrature_points();
             for (unsigned int q=0; q<transom_n_q_points; ++q)
                 {
                 Point<dim> local_press_force = transom_pressure_quad_values[q]*transom_node_normals[q];
                 transom_press_force += (local_press_force) * transom_fe_v.JxW(q);
                 transom_wet_surface += 1 * transom_fe_v.JxW(q);
                 //cout<<q<<"  F = ("<<local_press_force<<")    n = ("<<transom_node_normals[q]<<")"<<"  "<<transom_fe_v.JxW(q);
                 //cout<<"  p = ("<<transom_quad_nodes[q]<<")"<<endl;
                 }
             }
          }
      else
          {
          transom_press_force = Point<3>(0.0,0.0,0.0);
          transom_wet_surface = 0;
          }
      //std::string filename = ( output_file_name + "_" +
      //                         Utilities::int_to_string(0) +
      //                         ".vtu" );

      DataOut<dim-1, DoFHandler<dim-1, dim> > dataout;
      dataout.attach_dof_handler(transom_dh);
      dataout.add_data_vector(transom_pressure, "transom_pressure");
      dataout.build_patches(StaticMappingQ1<2,3>::mapping,transom_fe.degree,
                            DataOut<dim-1, DoFHandler<dim-1, dim> >::curved_inner_cells);     
      std::ofstream file("transom.vtu");
      dataout.write_vtu(file);
      }


  // viscous force computed with flat plate or ITTC 1957 model
  
  //water kinematic viscosity: 1.307e-6 m^2/s at 10°C; 1.004e-6 m^2/s at 20°C 
  double Re_l = fmax(1e-5,Vinf*comp_dom.Lx_boat/1.307e-6);
  //double Cd = 1.3282/sqrt(Re_l); // for flat plate
  double Cd = 0.075/(pow(log10(Re_l)-2., 2));
  double viscous_drag = 0.5*rho*Vinf*Vinf*comp_dom.boat_model.boatWetSurface*Cd;
  Point<3> visc_force(viscous_drag,0.0,0.0);

  std::string press_filename = ( output_file_name + "_force.txt" );

  ofstream myfile;
  if ( fabs(t-initial_time) < 1e-5 )
     myfile.open (press_filename.c_str());
  else
     myfile.open (press_filename.c_str(),ios::app);
  myfile << t <<" "<<press_force<<" "<<transom_press_force<<" "<<visc_force<<" "<<break_press_force<<" "<<press_moment<<" "<<wet_surface+transom_wet_surface<<" "<<press_force_test_1<<" "<<press_force_test_2<<" \n";
  myfile.close();

  std::cout<<"Total Force: "<<press_force-transom_press_force+visc_force+break_press_force<<std::endl;
  std::cout<<"Pressure Force: "<<press_force<<std::endl;
  std::cout<<"Breaking Wave Pressure Force: "<<break_press_force<<std::endl;
  std::cout<<"Transom Pressure Force: "<<transom_press_force<<std::endl;
  std::cout<<"Viscous Force: "<<visc_force<<std::endl;
  std::cout<<"Pressure Moment: "<<press_moment<<std::endl;
  std::cout<<"Wet Surface: "<<wet_surface<<"   Transom Wet Surface: "<<transom_wet_surface<<" ("<<ref_transom_wet_surface<<")"<<std::endl;

  double max_err = 0.0; 
  for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
      {
      if (comp_dom.flags[i] & pressure)
         {
         double k=0.62994; double omega=2.4835; double h=5.5;  double a=0.05;
         double time_factor = 1.0;
         double time_factor_deriv = 0.0;
         double ramp_length = 20.0;
         if (t<ramp_length)
            {
            time_factor = 0.5*sin(3.141592654*(t)/ramp_length-3.141592654/2)+0.5;
            time_factor_deriv = 0.5*3.141592654/ramp_length*cos(3.141592654*(t)/ramp_length-3.141592654/2);
            }                  
         double dphi_dt = omega*omega*a/k*cosh(k*(comp_dom.support_points[i](2)+h))/sinh(k*h)*cos(k*comp_dom.support_points[i](0)+omega*t)*time_factor +
                          omega*a/k*cosh(k*(comp_dom.support_points[i](2)+h))/sinh(k*h)*sin(k*comp_dom.support_points[i](0)+omega*t)*time_factor_deriv;
         Point<3> grad_phi(omega*a*cosh(k*(comp_dom.support_points[i](2)+h))/sinh(k*h)*cos(k*comp_dom.support_points[i](0)+omega*t)*time_factor,
                           0.0*time_factor,
                           omega*a*sinh(k*(comp_dom.support_points[i](2)+h))/sinh(k*h)*sin(k*comp_dom.support_points[i](0)+omega*t)*time_factor);



         double ex_pressure = ( -dphi_dt - 0.5*(grad_phi*grad_phi) - 
                                grad_phi(0)*instantWindValue(0)-grad_phi(1)*instantWindValue(1)-grad_phi(2)*instantWindValue(2) -
                                comp_dom.support_points[i](2)*g ) * rho;
         cout<<i<<":  P=("<<comp_dom.support_points[i]<<")   p_ex="<<ex_pressure<<" vs p="<<press(i)<<"   err="<<fabs(ex_pressure-press(i))<<endl;
         max_err = fmax(fabs(ex_pressure-press(i)),max_err);
         }
      }
  cout<<"Max Err: "<<max_err<<endl;
/*
// test of drag computation from control box momentum balance
  Vector<double> x_coors(comp_dom.dh.n_dofs());
  Vector<double> y_coors(comp_dom.dh.n_dofs());
  Vector<double> z_coors(comp_dom.dh.n_dofs());  
  for (unsigned int i=0; i<comp_dom.dh.n_dofs(); i++)
      {
      x_coors(i) = comp_dom.support_points[i](0);
      y_coors(i) = comp_dom.support_points[i](1);
      z_coors(i) = comp_dom.support_points[i](2);
      }

  dealii::Functions::FEFieldFunction<3,DoFHandler<2,3>, Vector <double> > fe_function(comp_dom.dh, x_coors, *comp_dom.mapping);
  
//void Functions::FEFieldFunction< dim, DH, VECTOR >::value_list	(	const std::vector< Point< dim > > & 	points,
//std::vector< double > & 	values,
//const unsigned int 	component = 0 
//)		 const
*/





	 
  std::cout<<"...done computing pressure. "<<std::endl;

}


template <int dim>
void FreeSurface<dim>::output_results(const std::string filename,
                                      const double t, 
                                      const Vector<double> & solution,
                                      const Vector<double> & solution_dot)
{ 


   VectorView<double> phi(comp_dom.dh.n_dofs(),solution.begin()+comp_dom.vector_dh.n_dofs());
   VectorView<double> phi_dot(comp_dom.dh.n_dofs(),solution_dot.begin()+comp_dom.vector_dh.n_dofs());
   VectorView<double> dphi_dn(comp_dom.dh.n_dofs(),solution.begin()+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs());
   VectorView<double> nodes_vel(comp_dom.vector_dh.n_dofs(),solution_dot.begin());

   wind.set_time(t);
   Vector<double> instantWindValue(dim);
   Point<dim> zero(0,0,0);
   wind.vector_value(zero,instantWindValue);
   Point<dim> Vinf;
   for (unsigned int i = 0; i < dim; i++)
       Vinf(i) = instantWindValue(i);

   Point<3> hull_lin_vel;
   Point<3> hull_ang_vel;
   Point<3> hull_lin_displ;
  for (unsigned int d=0; d<3;++d)
      {
      hull_lin_vel(d) = solution(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+d);
      hull_ang_vel(d) = solution(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+6+d);
      hull_lin_displ(d) = solution(comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()+3+d);
      }

   Point<3> baricenter_pos(hull_lin_displ(0)+comp_dom.boat_model.reference_hull_baricenter(0),
                           hull_lin_displ(1)+comp_dom.boat_model.reference_hull_baricenter(1),
                           hull_lin_displ(2)+comp_dom.boat_model.reference_hull_baricenter(2));



   Vector<double> fluid_vel_x(comp_dom.dh.n_dofs());
   Vector<double> fluid_vel_y(comp_dom.dh.n_dofs());
   Vector<double> fluid_vel_z(comp_dom.dh.n_dofs());
   Vector<double> nodes_vel_x(comp_dom.dh.n_dofs());
   Vector<double> nodes_vel_y(comp_dom.dh.n_dofs());
   Vector<double> nodes_vel_z(comp_dom.dh.n_dofs());
   //Vector<double> map_init_x(comp_dom.dh.n_dofs());
   //Vector<double> map_init_y(comp_dom.dh.n_dofs());
   //Vector<double> map_init_z(comp_dom.dh.n_dofs());

   Vector<double> complete_potential_gradients(comp_dom.vector_dh.n_dofs()); 
   compute_potential_gradients(complete_potential_gradients,phi,dphi_dn);

   for (unsigned int i = 0; i < comp_dom.dh.n_dofs(); i++)
       {
       fluid_vel_x(i) = complete_potential_gradients(i*dim)+Vinf(0);
       fluid_vel_y(i) = complete_potential_gradients(i*dim+1)+Vinf(1);
       fluid_vel_z(i) = complete_potential_gradients(i*dim+2)+Vinf(2);
       if ( ((comp_dom.flags[i] & boat) &&
            !(comp_dom.flags[i] & near_water) ) ||
             (comp_dom.flags[i] & transom_on_water) )
          {

          Point<3> rig_vel(hull_ang_vel(1)*(comp_dom.support_points[i](2)-hull_lin_displ(2))-
                           hull_ang_vel(2)*(comp_dom.support_points[i](1)-hull_lin_displ(1))+hull_lin_vel(0),
                           hull_ang_vel(2)*(comp_dom.support_points[i](0)-hull_lin_displ(0))-
                           hull_ang_vel(0)*(comp_dom.support_points[i](2)-hull_lin_displ(2))+hull_lin_vel(1),
                           hull_ang_vel(0)*(comp_dom.support_points[i](1)-hull_lin_displ(1))-
                           hull_ang_vel(1)*(comp_dom.support_points[i](0)-hull_lin_displ(0))+hull_lin_vel(2));

          nodes_vel_x(i) = nodes_vel(i*dim)+rig_vel(0);
          nodes_vel_y(i) = nodes_vel(i*dim+1)+rig_vel(1);
          nodes_vel_z(i) = nodes_vel(i*dim+2)+rig_vel(2);
          }
       else
          {
          nodes_vel_x(i) = nodes_vel(i*dim);
          nodes_vel_y(i) = nodes_vel(i*dim+1);
          nodes_vel_z(i) = nodes_vel(i*dim+2);
          }
       //map_init_x(i) = comp_dom.initial_map_points(i*dim);
       //map_init_y(i) = comp_dom.initial_map_points(i*dim+1);
       //map_init_z(i) = comp_dom.initial_map_points(i*dim+2);       
       }


  Vector<double> elevations(comp_dom.dh.n_dofs());  
  for (unsigned int i=0; i<comp_dom.dh.n_dofs(); i++)
      {
      elevations(i) = solution(dim*i+dim-1);
      }

  Vector<double> pressure(comp_dom.dh.n_dofs()); 
  Vector<double> comp_1(comp_dom.dh.n_dofs());
  Vector<double> comp_2(comp_dom.dh.n_dofs());
  Vector<double> comp_3(comp_dom.dh.n_dofs());
  Vector<double> comp_4(comp_dom.dh.n_dofs());
  compute_pressure(pressure,comp_1,comp_2,comp_3,comp_4,t,solution,solution_dot);

   // preparing iges normal vectors vector for pressure computation
   Vector<double> iges_normals_x_values(comp_dom.dh.n_dofs());
   Vector<double> iges_normals_y_values(comp_dom.dh.n_dofs());
   Vector<double> iges_normals_z_values(comp_dom.dh.n_dofs());
   for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
       {
       iges_normals_x_values(i) = comp_dom.iges_normals[i](0);
       iges_normals_y_values(i) = comp_dom.iges_normals[i](1);
       iges_normals_z_values(i) = comp_dom.iges_normals[i](2);
       }

  {
  DataOut<dim-1, DoFHandler<dim-1, dim> > dataout;
  
  dataout.attach_dof_handler(comp_dom.dh);


  dataout.add_data_vector((const Vector<double>&)phi, "phi");
  dataout.add_data_vector((const Vector<double>&)phi_dot, "phi_dot");
  dataout.add_data_vector(elevations, "elevations");
  dataout.add_data_vector(pressure, "pressure");
  dataout.add_data_vector((const Vector<double>&)dphi_dn, "dphi_dn");
  dataout.add_data_vector(fluid_vel_x, "fluid_vel_x");
  dataout.add_data_vector(fluid_vel_y, "fluid_vel_y");
  dataout.add_data_vector(fluid_vel_z, "fluid_vel_z");
  dataout.add_data_vector(nodes_vel_x, "nodes_vel_x");
  dataout.add_data_vector(nodes_vel_y, "nodes_vel_y");
  dataout.add_data_vector(nodes_vel_z, "nodes_vel_z");
  //cout<<"####### "<<break_wave_press.size()<<endl;
  //cout<<"#-----# "<<nodes_vel_z.size()<<endl;
  dataout.add_data_vector(break_wave_press, "break_wave_press");
  dataout.add_data_vector(comp_1, "comp_1");
  dataout.add_data_vector(comp_2, "comp_2");
  dataout.add_data_vector(comp_3, "comp_3");
  dataout.add_data_vector(comp_4, "comp_4");
  dataout.add_data_vector(iges_normals_x_values, "iges_normals_x");
  dataout.add_data_vector(iges_normals_y_values, "iges_normals_y");
  dataout.add_data_vector(iges_normals_z_values, "iges_normals_z");
  //dataout.add_data_vector(map_init_x, "map_init_x");
  //dataout.add_data_vector(map_init_y, "map_init_y");
  //dataout.add_data_vector(map_init_z, "map_init_z");
  //dataout.add_data_vector(DphiDt_sys_solution_2, "damping");

  dataout.build_patches(*comp_dom.mapping,
			comp_dom.vector_fe.degree,
			DataOut<dim-1, DoFHandler<dim-1, dim> >::curved_inner_cells);     

  std::ofstream file(filename.c_str());
  
  dataout.write_vtu(file);
  }

  std::string hull_motions_filename = ( output_file_name + "_hull_motions.txt" );

  ofstream myfile;
  if ( fabs(t-initial_time) < 1e-5 )
     myfile.open (hull_motions_filename.c_str());
  else
     myfile.open (hull_motions_filename.c_str(),ios::app);
  myfile << t <<" "<<baricenter_pos<<" "<<comp_dom.boat_model.yaw_angle<<" "<<-comp_dom.boat_model.pitch_angle+comp_dom.boat_model.initial_trim<<" "<<-comp_dom.boat_model.roll_angle<<" "<<hull_lin_vel<<" "<<hull_ang_vel<<" \n";
  myfile.close();


/*
Standard_Integer NbPointConstraint=1;
// Initialize a BuildPlateSurface
GeomPlate_BuildPlateSurface BPSurf(4,15,4);
// Point constraints
for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
    {
    if (comp_dom.flags[i] & water)
       {
       Handle(GeomPlate_PointConstraint) PCont= new GeomPlate_PointConstraint(Pnt(comp_dom.support_points[i]),0);
       BPSurf.Add(PCont);
       }
    }
// Compute the Plate surface
BPSurf.Perform();



// Approximation of the Plate surface
Standard_Integer MaxSeg=100;
Standard_Integer MaxDegree=8;
Standard_Integer CritOrder=0;
Standard_Real dmax,Tol;
Handle(GeomPlate_Surface) PSurf = BPSurf.Surface();
//dmax = Max(0.0001,10*BPSurf.G0Error());
dmax = 0.00001;
Tol=0.00001;
GeomPlate_MakeApprox
Mapp(PSurf,Tol,MaxSeg,MaxDegree,dmax,CritOrder);
Handle (Geom_Surface) Surf (Mapp.Surface());
// create a face corresponding to the approximated Plate Surface
Standard_Real Umin, Umax, Vmin, Vmax;
PSurf->Bounds( Umin, Umax, Vmin, Vmax);
BRepBuilderAPI_MakeFace MF(Surf,Umin, Umax, Vmin, Vmax,1e-4);

TopoDS_Shape cut_edge;
intersect_plane(MF,cut_edge,1.0,0.0,0.0,-2.32422);
// These lines can be used to dump the surface on an .igs file

IGESControl_Controller::Init();
IGESControl_Writer ICW ("MM", 0);
Standard_Boolean ok = ICW.AddShape (MF);
ICW.AddShape (cut_edge);
//Standard_Boolean ok = ICW.AddShape (Pnt(comp_dom.support_points[0]));
ICW.ComputeModel();
Standard_Boolean OK = ICW.Write ("free_surf.igs");
*/

//compute_internal_velocities(phi,dphi_dn);


}


template <int dim>
void FreeSurface<dim>::compute_internal_velocities(const Vector<double> &phi,
                                                   const Vector<double> &dphi_dn)
{

   unsigned int n_points;

   // Create streamobject
   ifstream infile;
   infile.open("points.txt");

   // Exit if file opening failed
   if (!infile.is_open())
      {
      cerr<<"Opening failed"<<endl;
      exit(1);
      }

    // Get number of points
    infile >> n_points;

    std::vector< Point<3> > points(n_points,Point<3>(0.0,0.0,0.0));
    std::vector< Point<3> > velocities(n_points,Point<3>(0.0,0.0,0.0));

    unsigned int count = 0;
    while (!infile.eof())
          {
          infile >> points[count](0) >> points[count](1) >> points[count](2);
          ++count;
          }
    infile.close();
    
 


   comp_dom.update_support_points();


   FEValues<dim-1,dim> fe_v(*comp_dom.mapping, comp_dom.fe, *comp_dom.quadrature,
			    update_values | update_gradients |
			    update_cell_normal_vectors |
			    update_quadrature_points |
			    update_JxW_values);


   const unsigned int n_q_points = fe_v.n_quadrature_points;
   const unsigned int  dofs_per_cell   = comp_dom.fe.dofs_per_cell;
   std::vector<unsigned int> local_dof_indices (dofs_per_cell);
      
   std::vector<double> q_phi(n_q_points);
   std::vector<double> q_dphi_dn(n_q_points);
   


   cell_it
   cell = comp_dom.dh.begin_active(),
   endc = comp_dom.dh.end();

   for (; cell!=endc; ++cell)
     {
     fe_v.reinit(cell);

     fe_v.get_function_values(dphi_dn, q_dphi_dn);
     fe_v.get_function_values(phi, q_phi);
     const std::vector<Point<dim> > &quad_nodes = fe_v.get_quadrature_points();
     const std::vector<Point<dim> > &quad_nodes_normals = fe_v.get_normal_vectors();

     for (unsigned int i=0; i<n_points; ++i)
         {
         fad_double x,y,z;
         x = points[i](0);
         y = points[i](1);
         z = points[i](2);

         x.diff(0,3);  
         y.diff(1,3);
         z.diff(2,3);

         for (unsigned int q=0; q<n_q_points; ++q)
             {

             Point <dim, fad_double > r(quad_nodes[q](0)-x,quad_nodes[q](1)-y,quad_nodes[q](2)-z);

             fad_double G = fad_double(1.0/(4.0*numbers::PI))/r.norm();

             fad_double dG_dn = -(r(0)*fad_double(quad_nodes_normals[q](0))+
                                  r(1)*fad_double(quad_nodes_normals[q](1))+
                                  r(2)*fad_double(quad_nodes_normals[q](2)))/(fad_double(4.0*numbers::PI)*pow(r.norm(),3.0));

     
             //cout<<"G: "<<G.val()<<"  "<<G_x_plus<<"  "<<G_z_plus<<endl;
         
             Point<dim> grad_G(G.fastAccessDx(0),G.fastAccessDx(1),G.fastAccessDx(2));
             Point<dim> grad_dG_dn(dG_dn.fastAccessDx(0),dG_dn.fastAccessDx(1),dG_dn.fastAccessDx(2));

             velocities[i] += (q_dphi_dn[q]*grad_G - q_phi[q]*grad_dG_dn)*fe_v.JxW(q);

             }

         }

     }


  for (unsigned int i=0; i<n_points; ++i)
      cout<<i<<"  P("<<points[i]<<")   grad_phi("<<velocities[i]<<") "<<endl;

  ofstream myfile;
  myfile.open ("velocities.txt");
  for (unsigned int i=0; i<n_points; ++i)
      myfile<<velocities[i]<<endl;
  myfile.close();
}

template <int dim>
void FreeSurface<dim>::compute_constraints(ConstraintMatrix &c,
					   ConstraintMatrix &cc) {
   std::vector<Point<dim> > supp_points(comp_dom.vector_dh.n_dofs());
   DoFTools::map_dofs_to_support_points<dim-1, dim>( *comp_dom.mapping, 
 						    comp_dom.vector_dh, 
 						    supp_points);

  // we start clearing the constraint matrices
  c.clear();
  cc.clear();


  // here we prepare the constraint matrices so as to account for the presence hanging
  // nodes

  DoFTools::make_hanging_node_constraints (comp_dom.dh,c);
  DoFTools::make_hanging_node_constraints (comp_dom.vector_dh,cc);
/*
   for(unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
      {
      std::cout << i << std::endl;
      if (c.is_constrained(i))
         {
         std::cout << "Constraining " << i << std::endl;
       cc.add_line(i);
       c.add_line(i);
       cc.set_inhomogeneity(i, 0);
         }
      }
  //*/

  cc.close();
  c.close();
}


template <int dim>
void FreeSurface<dim>::dump_solution(const Vector<double> &y,
				     const Vector<double> &yp,
				     const std::string fname) const
{
  std::cout << "Dumping solution: " << fname << std::endl;
  std::ofstream ofile((fname+"_y.dat").c_str());
  y.block_write(ofile);
  ofile.close();
  ofile.open((fname+"_yp.dat").c_str());
  yp.block_write(ofile);
  comp_dom.dump_tria(fname+"_tria.dat");
  ofile.close();
}



template <int dim>
void FreeSurface<dim>::restore_solution(Vector<double> &y,
					Vector<double> &yp,
					const std::string fname)
{
  std::cout << "Restoring solution: " << fname << std::endl;

  std::ifstream ifile((fname+"_y.dat").c_str());
  y.block_read(ifile);
  ifile.close();
  ifile.open((fname+"_yp.dat").c_str());
  yp.block_read(ifile);
  ifile.close();
  comp_dom.restore_tria(fname+"_tria.dat");
  comp_dom.update_mapping(y);


}


template <int dim>
void FreeSurface<dim>::enforce_full_geometry_constraints()
{
  for(unsigned int i=0; i<comp_dom.vector_dh.n_dofs(); ++i)
    {
      if(sys_comp(i) == 0)
	{
	  std::set<unsigned int> doubles = comp_dom.vector_double_nodes_set[i];
	  for (std::set<unsigned int>::iterator it = doubles.begin() ; it != doubles.end(); it++ )
	    {
					       //if(sys_comp(*it) == 1)
	      comp_dom.map_points(i) = comp_dom.map_points(*it);
					       //else
					       //  map_points(*it) = map_points(i);   
	    }
	}
    }

  comp_dom.full_mesh_treatment(); 
  vector_constraints.distribute(comp_dom.map_points);
}

template <int dim>
void FreeSurface<dim>::enforce_partial_geometry_constraints(const double blend_factor)
{
  for(unsigned int i=0; i<comp_dom.vector_dh.n_dofs(); ++i)
    {
      if(sys_comp(i) == 0)
	{
	  std::set<unsigned int> doubles = comp_dom.vector_double_nodes_set[i];
	  for (std::set<unsigned int>::iterator it = doubles.begin() ; it != doubles.end(); it++ )
	    {
					       //if(sys_comp(*it) == 1)
	      comp_dom.map_points(i) = comp_dom.map_points(*it);
					       //else
					       //  map_points(*it) = map_points(i);   
	    }
	}
    }
  comp_dom.partial_mesh_treatment(blend_factor);
  vector_constraints.distribute(comp_dom.map_points);
}




template <int dim>
void FreeSurface<dim>::compute_keel_smoothing(Vector<double> & map_points)
{

}


  

template class FreeSurface<3>;

