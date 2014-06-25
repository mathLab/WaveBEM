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
//    Authors: Luca Heltai, Cataldo Manigrasso, Andrea Mola
//
//----------------------------  step-34.cc  ---------------------------

#define TOLL 0.001
#define MAXELEMENTSPERBLOCK 1

#include "../include/numerical_towing_tank.h"
#include "../include/boat_surface.h"
#include <numerics/matrix_tools.h>
#include <grid/grid_refinement.h>

using namespace dealii;
using namespace OpenCascade;
using namespace std;

NumericalTowingTank::NumericalTowingTank(const unsigned int fe_degree,
					 const unsigned int mapping_degree) :
		ComputationalDomain<3>(fe_degree, mapping_degree),
                restart_surface_smoother(NULL),
                surface_smoother(NULL),
                line_smoothers(7,NULL),
                mapping_degree(mapping_degree)
{}

void NumericalTowingTank::full_mesh_treatment()

{
std::cout<<"Performing full mesh treatment"<<std::endl;

perform_line_smoothing(7);

perform_water_line_nodes_projection();
// for (unsigned int i=0; i<vector_dh.n_dofs(); ++i)

//perform_surface_projection();
perform_smoothing(true,1);
perform_surface_projection();

vector_constraints.distribute(map_points);
// writing boat/free surface nodes to a file
update_support_points();
//std::ofstream myfile;
//myfile.open ("waterLine.txt",std::ios::trunc);
//for (std::set <unsigned int>::iterator pos = free_surf_and_boat_nodes.begin(); pos != free_surf_and_boat_nodes.end();  pos++)
//    {
//       myfile << support_points[*pos]<<" \n";       
//    }    
//myfile.close();

std::cout<<"done full mesh treatment"<<std::endl;
}

void NumericalTowingTank::partial_mesh_treatment(const double blend_factor)

{
std::cout<<"Performing partial mesh treatment"<<std::endl;

perform_line_smoothing(3);
perform_water_line_nodes_projection();
perform_smoothing(false,blend_factor);
vector_constraints.distribute(map_points);

std::cout<<"done partial mesh treatment"<<std::endl;
}

void NumericalTowingTank::apply_curvatures(const Vector<double> &curvatures, const vector<bool> boundary_dofs)

{
std::cout<<"Computing geometrically conforming mesh"<<std::endl;

smoothing_map_points = map_points;

restart_surface_smoother->apply_curvatures(curvatures,boundary_dofs);

for (unsigned int i=0; i<vector_dh.n_dofs()/3;++i) 
    {
    if ((flags[i] & water) == 0)
       {
       map_points(3*i) = smoothing_map_points(3*i);
       map_points(3*i+1) = smoothing_map_points(3*i+1);
       map_points(3*i+2) = smoothing_map_points(3*i+2);  
       }
    }


std::cout<<"Done computing geometrically conforming mesh"<<std::endl;
}


void NumericalTowingTank::update_mapping(const Vector<double> &new_disp)
{
  for (unsigned int i=0; i<vector_dh.n_dofs(); ++i)
       map_points(i) = new_disp(i);
}

void NumericalTowingTank::compute_curvatures(Vector<double> &curvatures)

{
std::cout<<"Computing curvatures for geometrically conforming mesh refinement"<<std::endl;

smoothing_map_points = map_points;

restart_surface_smoother->compute_curvatures(curvatures);


std::cout<<"Done computing curvatures for geometrically conforming mesh"<<std::endl;
}

void NumericalTowingTank::refine_and_resize()
{
  std::cout<<"Refining and resizing mesh as required"<<std::endl;  
  
  dh.distribute_dofs(fe);
  vector_dh.distribute_dofs(vector_fe);
  
  std::cout<<"Total number of dofs of first mesh: "<<dh.n_dofs()<<std::endl;
  
  map_points.reinit(vector_dh.n_dofs());
  smoothing_map_points.reinit(vector_dh.n_dofs());
  old_map_points.reinit(vector_dh.n_dofs()); 
  ref_points.resize(vector_dh.n_dofs());
  DoFTools::map_dofs_to_support_points<2,3>(StaticMappingQ1<2,3>::mapping,
					    vector_dh, ref_points);
 
  mapping = new MappingQEulerian<2, Vector<double>, 3>
	    (mapping_degree, map_points, vector_dh);
  
  generate_double_nodes_set();
  full_mesh_treatment();
// anisotropic refinement
  remove_mesh_anisotropy(tria);
  refine_global_on_boat(init_global_boat_refs);

  Vector<double> positions(vector_dh.n_dofs());
  std::vector< Point<3> > normals(vector_dh.n_dofs()/3);
  Vector<float> estimated_error_per_cell(tria.n_active_cells());

  update_support_points();
  for (unsigned int i=0; i<vector_dh.n_dofs()/3; ++i)
      {
      for(unsigned int j=0; j<3; ++j) 
       	 {
       	 positions(3*i+j) = vector_support_points[3*i][j];
       	 }
      }


  Point<3> projection;

  for (unsigned int k=0; k<init_adaptive_boat_refs; ++k)
      {
      std::cout<<"Adaptive refinement cycle "<<k+1<<" of "<<init_adaptive_boat_refs<<std::endl;
      estimated_error_per_cell.reinit(tria.n_active_cells());    
      /*
      QGauss<1> quad(2);

      KellyErrorEstimator<2,3>::estimate (vector_dh,
					  quad,
					  FunctionMap<3>::type(),
					  positions,
					  estimated_error_per_cell);

      */
      QGauss<2> quad(1);

      FEValues<2,3> fe_v(*mapping, fe, quad,
	                 update_values | update_cell_normal_vectors | update_quadrature_points | 
		         update_JxW_values);

      const unsigned int n_q_points = fe_v.n_quadrature_points;
      unsigned int cell_count = 0; 
      cell_it
      cell = dh.begin_active(),
      endc = dh.end();
   
      for (; cell!=endc; ++cell)
          {
          fe_v.reinit (cell);
          const std::vector<Point<3> > &node_normals = fe_v.get_normal_vectors();
          const std::vector<Point<3> > &quadrature_points = fe_v.get_quadrature_points();
          if ((cell->material_id() == wall_sur_ID1 ||
               cell->material_id() == wall_sur_ID2 ||
               cell->material_id() == wall_sur_ID3 ) )
             {
             Point<3> proj_node;
             if (quadrature_points[0](1) > 0) //right side
                {
                boat_model.boat_water_line_right->assigned_axis_projection(proj_node,
                                                                           quadrature_points[0],
                                                                           node_normals[0]);  // for projection in mesh normal direction
                //cout<<cell<<"  center:"<<quadrature_points[0]<<" normal: " <<node_normals[0]<<endl;
                //cout<<"Projection: "<<proj_node<<" distance: "<<proj_node.distance(quadrature_points[0])<<endl;
                //cout<<endl;
                estimated_error_per_cell(cell_count) = proj_node.distance(quadrature_points[0]);
                }
             else  // left side
                {
                boat_model.boat_water_line_left->assigned_axis_projection(proj_node,
                                                                           quadrature_points[0],
                                                                           node_normals[0]);  // for projection in mesh normal direction
                //cout<<cell<<"  center:"<<quadrature_points[0]<<" normal: " <<node_normals[0]<<endl;
                //cout<<"Projection: "<<proj_node<<" distance: "<<proj_node.distance(quadrature_points[0])<<endl;
                //cout<<endl;
                estimated_error_per_cell(cell_count) = proj_node.distance(quadrature_points[0]);
                }
             }
          ++cell_count;
          }



      GridRefinement::refine_and_coarsen_fixed_fraction	(tria,
                                                         estimated_error_per_cell,
                                                         init_adaptive_boat_refs_fraction,
							 0.0,
							 5000);

      tria.prepare_coarsening_and_refinement();
      tria.execute_coarsening_and_refinement();

      dh.distribute_dofs(fe);
      vector_dh.distribute_dofs(vector_fe);  
      map_points.reinit(vector_dh.n_dofs());
      smoothing_map_points.reinit(vector_dh.n_dofs());
      old_map_points.reinit(vector_dh.n_dofs());
      ref_points.resize(vector_dh.n_dofs());
      DoFTools::map_dofs_to_support_points<2,3>(StaticMappingQ1<2,3>::mapping,
					        vector_dh, ref_points);
      generate_double_nodes_set();
      make_edges_conformal(true);
      make_edges_conformal(true);
      full_mesh_treatment();

      update_support_points();
      positions.reinit(vector_dh.n_dofs());

      for (unsigned int i=0; i<vector_dh.n_dofs()/3; ++i)
          for(unsigned int j=0; j<3; ++j) 
       	     positions(3*i+j) = vector_support_points[3*i][j];
      }

      remove_transom_hanging_nodes();


      min_diameter = 10000;
      Triangulation<2,3>::active_cell_iterator
      cell = tria.begin_active(), endc = tria.end();

      for ( ; cell != endc; ++cell)
          {
          if (cell->material_id() == free_sur_ID1 ||
	      cell->material_id() == free_sur_ID2 ||
	      cell->material_id() == free_sur_ID3   )
          min_diameter = std::min(min_diameter,cell->diameter());
          }  
      std::cout << "Min diameter: << " << min_diameter << std::endl;

        

  std::cout<<"Total number of dofs after whole refinement: "<<dh.n_dofs()<<std::endl;

  std::cout<<"...done refining and resizing mesh"<<std::endl;



}
  
void NumericalTowingTank::read_domain ()
{
  //tria.set_mesh_smoothing (Triangulation<2,3>::do_not_produce_unrefined_islands );
  //coarse_tria.set_mesh_smoothing (Triangulation<2,3>::do_not_produce_unrefined_islands );

  boat_model.start_iges_model(iges_file_name,
                              1.0/1000.0,
                              boat_displacement,
                              assigned_sink,
                              assigned_trim,
                              back_keel_length, 
                              front_keel_length,
                              middle_keel_length,
                              number_of_transom_edges);
  //boat_model.start_iges_model(iges_file_name,1.0/1000.0,boat_displacement,assigned_sink,assigned_trim,0.068979, 0.075, .5); //con kcs.iges
  //boat_model.start_iges_model(iges_file_name, 1.0/1000.0,boat_displacement,assigned_sink,assigned_trim,0.13, 0.05, .5); //con goteborg.iges
  //boat_model.start_iges_model(iges_file_name, 1.0/1000.0,boat_displacement,assigned_sink,assigned_trim,0.13, .05, .47); //con goteborgLow.iges
  //boat_model.start_iges_model(iges_file_name, 1.0/1000.0,boat_displacement,assigned_sink,assigned_trim,0.13, .075, .5); //con goteborgTransom.iges
  //boat_model.start_iges_model(iges_file_name, 1.0/1000.0,boat_displacement,assigned_sink,assigned_trim,0.055556, 0.055556, .5); //con wigley.iges
  //boat_model.start_iges_model(iges_file_name, 1.0/1000.0,boat_displacement,assigned_sink,assigned_trim,0.068979, 0.049807, .5); //con series_60.iges
  //boat_model.start_iges_model(iges_file_name, 1.0/1000.0,boat_displacement,assigned_sink,assigned_trim,0.35, 0.2, .5); //con vela_enave.iges

  cout<<"FrontBot "<<boat_model.PointFrontBot<<endl;
  cout<<"BackBot "<<boat_model.PointBackBot<<endl;
  create_initial_mesh(boat_model.PointFrontTop,
                      boat_model.PointFrontBot,
                      boat_model.PointMidTop,
                      boat_model.PointMidBot,
                      boat_model.PointBackTop,
                      boat_model.PointBackBot,
                      boat_model.PointLeftTransom,
                      boat_model.PointRightTransom,
                      boat_model.PointCenterTransom,
                      coarse_tria);

  tria.restore();
 
/*
  Triangulation<2,3>::active_cell_iterator
    cell = tria.begin_active(), endc = tria.end();

  for (cell=tria.begin_active(); cell!= endc;++cell)
      {
      if ((cell->material_id() == free_sur_ID1 ||
           cell->material_id() == free_sur_ID2 ||
           cell->material_id() == free_sur_ID3 ))
          {
	  if ( cell->at_boundary() )
             for (unsigned int f = 0; f < GeometryInfo<2>::faces_per_cell;++f)
                 {
                 if (cell->face(f)->boundary_indicator()==0)
                    cell->face(f)->set_boundary_indicator(22);
                 }
          }
      }
*/
  tria.set_boundary(5, *boat_model.undist_water_surf);
  tria.set_boundary(6, *boat_model.undist_water_surf);
  //tria.set_boundary(21, *boat_model.boat_surface_right);
  tria.set_boundary(21, *boat_model.water_line_right); //tria.set_boundary(21, *boat_model.boat_water_line_right);
  //tria.set_boundary(23, *boat_model.boat_surface_right);
  tria.set_boundary(23, *boat_model.water_line_right);//tria.set_boundary(23, *boat_model.boat_water_line_right);
  //tria.set_boundary(22, *boat_model.boat_surface_left);
  tria.set_boundary(22, *boat_model.water_line_left);//tria.set_boundary(22, *boat_model.boat_water_line_left);
  //tria.set_boundary(24, *boat_model.boat_surface_left);
  tria.set_boundary(24, *boat_model.water_line_left);//tria.set_boundary(22, *boat_model.boat_water_line_left);
  //tria.set_boundary(26, *boat_model.boat_surface_right);
  tria.set_boundary(26, *boat_model.water_line_right);//tria.set_boundary(26, *boat_model.boat_water_line_right);
  //tria.set_boundary(28, *boat_model.boat_surface_right);
  tria.set_boundary(28, *boat_model.water_line_right);//tria.set_boundary(28, *boat_model.boat_water_line_right);
  //tria.set_boundary(27, *boat_model.boat_surface_left);
  tria.set_boundary(27, *boat_model.water_line_left);//tria.set_boundary(27, *boat_model.boat_water_line_left);
  //tria.set_boundary(29, *boat_model.boat_surface_left);
  tria.set_boundary(29, *boat_model.water_line_left);//tria.set_boundary(29, *boat_model.boat_water_line_left);
  if (boat_model.is_transom)
     tria.set_boundary(32, *boat_model.boat_transom_left);
  else
     tria.set_boundary(32, *boat_model.boat_keel);
  tria.set_boundary(31, *boat_model.boat_keel);
  tria.set_boundary(30, *boat_model.boat_keel);
  if (boat_model.is_transom)
     tria.set_boundary(37, *boat_model.boat_transom_right);
  else
     tria.set_boundary(37, *boat_model.boat_keel);
  tria.set_boundary(36, *boat_model.boat_keel);
  tria.set_boundary(35, *boat_model.boat_keel);
  if (boat_model.is_transom)
     {
     tria.set_boundary(40, *boat_model.boat_transom_left);
     tria.set_boundary(41, *boat_model.boat_transom_right);
     }

  coarse_tria.set_boundary(5, *boat_model.undist_water_surf);
  coarse_tria.set_boundary(6, *boat_model.undist_water_surf);
  //coarse_tria.set_boundary(21, *boat_model.boat_surface_right);
  coarse_tria.set_boundary(21, *boat_model.water_line_right); //coarse_tria.set_boundary(21, *boat_model.boat_water_line_right);
  //coarse_tria.set_boundary(23, *boat_model.boat_surface_right);
  coarse_tria.set_boundary(23, *boat_model.water_line_right);//coarse_tria.set_boundary(23, *boat_model.boat_water_line_right);
  //coarse_tria.set_boundary(22, *boat_model.boat_surface_left);
  coarse_tria.set_boundary(22, *boat_model.water_line_left);//coarse_tria.set_boundary(22, *boat_model.boat_water_line_left);
  //coarse_tria.set_boundary(24, *boat_model.boat_surface_left);
  coarse_tria.set_boundary(24, *boat_model.water_line_left);//coarse_tria.set_boundary(22, *boat_model.boat_water_line_left);
  //coarse_tria.set_boundary(26, *boat_model.boat_surface_right);
  coarse_tria.set_boundary(26, *boat_model.water_line_right);//coarse_tria.set_boundary(26, *boat_model.boat_water_line_right);
  //coarse_tria.set_boundary(28, *boat_model.boat_surface_right);
  coarse_tria.set_boundary(28, *boat_model.water_line_right);//coarse_tria.set_boundary(28, *boat_model.boat_water_line_right);
  //coarse_tria.set_boundary(27, *boat_model.boat_surface_left);
  coarse_tria.set_boundary(27, *boat_model.water_line_left);//coarse_tria.set_boundary(27, *boat_model.boat_water_line_left);
  //coarse_tria.set_boundary(29, *boat_model.boat_surface_left);
  coarse_tria.set_boundary(29, *boat_model.water_line_left);//coarse_tria.set_boundary(29, *boat_model.boat_water_line_left);
    if (boat_model.is_transom)
     coarse_tria.set_boundary(32, *boat_model.boat_transom_left);
  else
     coarse_tria.set_boundary(32, *boat_model.boat_keel);
  coarse_tria.set_boundary(31, *boat_model.boat_keel);
  coarse_tria.set_boundary(30, *boat_model.boat_keel);
    if (boat_model.is_transom)
     coarse_tria.set_boundary(37, *boat_model.boat_transom_right);
  else
     coarse_tria.set_boundary(37, *boat_model.boat_keel);
  coarse_tria.set_boundary(36, *boat_model.boat_keel);
  coarse_tria.set_boundary(35, *boat_model.boat_keel);
  if (boat_model.is_transom)
     {
     coarse_tria.set_boundary(40, *boat_model.boat_transom_left);
     coarse_tria.set_boundary(41, *boat_model.boat_transom_right);
     }
}

 
NumericalTowingTank::~NumericalTowingTank()
{
  for (unsigned int i=0; i<line_smoothers.size(); ++i)
      delete line_smoothers[i];
  if (surface_smoother)
      delete surface_smoother;

  //tria.set_boundary(3);
  //tria.set_boundary(4);
  tria.set_boundary(21);
  tria.set_boundary(23);
  tria.set_boundary(22);
  tria.set_boundary(24);
  tria.set_boundary(26);
  tria.set_boundary(28);
  tria.set_boundary(27);
  tria.set_boundary(29);
  tria.set_boundary(32);
  tria.set_boundary(31);
  tria.set_boundary(30);
  tria.set_boundary(37);
  tria.set_boundary(36);
  tria.set_boundary(35);
  //coarse_tria.set_boundary(3);
  //coarse_tria.set_boundary(4);
  coarse_tria.set_boundary(21);
  coarse_tria.set_boundary(23);
  coarse_tria.set_boundary(22);
  coarse_tria.set_boundary(24);
  coarse_tria.set_boundary(26);
  coarse_tria.set_boundary(28);
  coarse_tria.set_boundary(27);
  coarse_tria.set_boundary(29);
  coarse_tria.set_boundary(32);
  coarse_tria.set_boundary(31);
  coarse_tria.set_boundary(30);
  coarse_tria.set_boundary(37);
  coarse_tria.set_boundary(36);
  coarse_tria.set_boundary(35);
  

}



  
void NumericalTowingTank::declare_parameters (ParameterHandler &prm)
{
  ComputationalDomain<3>::declare_parameters (prm);
// parameters to set domain 2ensions with respect to boat length could
// be asked here to the user		      
  prm.declare_entry("Iges file name", "wigley.iges", Patterns::Anything()); 

  prm.declare_entry("Refinement level on boat", "0", Patterns::Integer());

  prm.declare_entry("Boat displacement (Kg)","0",Patterns::Double());

  prm.declare_entry("Assigned sink (m)","0",Patterns::Double());

  prm.declare_entry("Assigned trim (rad)","0",Patterns::Double());

  prm.declare_entry("Max aspect ratio","1.5",Patterns::Double());

  prm.declare_entry("Front mesh inclination coeff","1.0",Patterns::Double());

  prm.declare_entry("Back mesh inclination coeff","1.0",Patterns::Double());

  prm.declare_entry("Back keel length","0.1",Patterns::Double());

  prm.declare_entry("Front keel length","0.05",Patterns::Double());

  prm.declare_entry("Middle keel length","0.5",Patterns::Double());

  prm.declare_entry("Number of boat initial uniform refinements", "2", Patterns::Integer());

  prm.declare_entry("Number of boat initial curvature based refinements", "0", Patterns::Integer());

  prm.declare_entry("Initial curvature based refinements fraction","0.2",Patterns::Double());

  prm.declare_entry("Number of transom edges", "1", Patterns::Integer());
}

  
void NumericalTowingTank::parse_parameters (ParameterHandler &prm)
{
  ComputationalDomain<3>::parse_parameters (prm);
// parameters to set domain dimensions with respect to boat length could
// be asked here to the user	
  n_cycles = prm.get_integer("Refinement level on boat");
  
  iges_file_name = prm.get("Iges file name");

  boat_displacement = prm.get_double("Boat displacement (Kg)");

  assigned_sink = prm.get_double("Assigned sink (m)");

  assigned_trim = prm.get_double("Assigned trim (rad)");

  max_aspect_ratio = prm.get_double("Max aspect ratio");

  front_mesh_inclination_coeff = prm.get_double("Front mesh inclination coeff");

  back_mesh_inclination_coeff = prm.get_double("Back mesh inclination coeff");

  back_keel_length = prm.get_double("Back keel length");

  front_keel_length = prm.get_double("Front keel length");

  middle_keel_length = prm.get_double("Middle keel length");

  init_global_boat_refs = prm.get_integer("Number of boat initial uniform refinements");

  init_adaptive_boat_refs = prm.get_integer("Number of boat initial curvature based refinements");

  init_adaptive_boat_refs_fraction = prm.get_double("Initial curvature based refinements fraction");

  number_of_transom_edges = prm.get_integer("Number of transom edges");
}


        
 
void NumericalTowingTank::create_initial_mesh(const Point<3> PointFrontTop,
                                              const Point<3> PointFrontBot,
                                              const Point<3> PointMidTop,
                                              const Point<3> PointMidBot,
                                              const Point<3> PointBackTop,
                                              const Point<3> PointBackBot,
                                              const Point<3> PointLeftTransom,
                                              const Point<3> PointRightTransom,
                                              const Point<3> PointCenterTransom,
                                              Triangulation<2,3>  &triangulation)
{


std::vector<Point<3> > vertices;
std::vector<CellData<2> > cells;
SubCellData subcelldata;

double a = front_mesh_inclination_coeff;
double b = back_mesh_inclination_coeff;

if (boat_model.is_transom)
{


Lx_boat = boat_model.boatWetLength;
Lx_domain = Lx_boat*12.0;
Ly_domain = Lx_boat*4.0;
Lz_domain = Lx_boat*2.0;

vertices.resize(88);
      
vertices[0](0)=-Lx_domain/2; vertices[0](1)=-Ly_domain/2 ; vertices[0](2)=0;
vertices[1](0)=-Lx_domain/2; vertices[1](1)=-Ly_domain/2; vertices[1](2)=-Lz_domain;
vertices[2](0)=Lx_domain/2; vertices[2](1)=-Ly_domain/2; vertices[2](2)=-Lz_domain;
vertices[3](0)=Lx_domain/2; vertices[3](1)=-Ly_domain/2; vertices[3](2)=0;
vertices[4](0)=b*PointLeftTransom(0); vertices[4](1)=-Ly_domain/2; vertices[4](2)=0;
vertices[5](0)=a*PointFrontTop(0); vertices[5](1)=-Ly_domain/2; vertices[5](2)=0;
vertices[6](0)=b*PointLeftTransom(0); vertices[6](1)=-Ly_domain/2; vertices[6](2)=-Lz_domain;
vertices[7](0)=0; vertices[7](1)=-Ly_domain/2; vertices[7](2)=-Lz_domain;
vertices[8](0)=a*PointFrontTop(0); vertices[8](1)=-Ly_domain/2; vertices[8](2)=-Lz_domain;
vertices[9](0)=0; vertices[9](1)=-Ly_domain/2; vertices[9](2)=0;
vertices[10](0)=Lx_domain/2; vertices[10](1)=Ly_domain/2; vertices[10](2)=0;
vertices[11](0)=Lx_domain/2; vertices[11](1)=Ly_domain/2; vertices[11](2)=-Lz_domain;
vertices[12](0)=-Lx_domain/2; vertices[12](1)=Ly_domain/2; vertices[12](2)=-Lz_domain;
vertices[13](0)=-Lx_domain/2; vertices[13](1)=Ly_domain/2; vertices[13](2)=0;
vertices[14](0)=a*PointFrontTop(0); vertices[14](1)=Ly_domain/2; vertices[14](2)=0;
vertices[15](0)=b*PointRightTransom(0); vertices[15](1)=Ly_domain/2; vertices[15](2)=0;
vertices[16](0)=a*PointFrontTop(0); vertices[16](1)=Ly_domain/2; vertices[16](2)=-Lz_domain;
vertices[17](0)=0; vertices[17](1)=Ly_domain/2; vertices[17](2)=-Lz_domain;
vertices[18](0)=b*PointRightTransom(0); vertices[18](1)=Ly_domain/2; vertices[18](2)=-Lz_domain;
vertices[19](0)=0; vertices[19](1)=Ly_domain/2; vertices[19](2)=0;
vertices[20](0)=PointFrontTop(0); vertices[20](1)=0; vertices[20](2)=0;
vertices[21](0)=PointRightTransom(0); vertices[21](1)=PointRightTransom(1); vertices[21](2)=0;
vertices[22](0)=PointFrontBot(0); vertices[22](1)=0; vertices[22](2)=PointFrontBot(2);
vertices[23](0)=PointCenterTransom(0); vertices[23](1)=0; vertices[23](2)=PointCenterTransom(2);
vertices[24](0)=0; vertices[24](1)=PointMidTop(1); vertices[24](2)=0;
vertices[25](0)=PointMidBot(0); vertices[25](1)=0; vertices[25](2)=PointMidBot(2);
vertices[26](0)=PointLeftTransom(0); vertices[26](1)=PointLeftTransom(1); vertices[26](2)=0;
vertices[27](0)=PointFrontTop(0); vertices[27](1)=0; vertices[27](2)=0;
vertices[28](0)=PointCenterTransom(0); vertices[28](1)=0; vertices[28](2)=PointCenterTransom(2);
vertices[29](0)=PointFrontBot(0); vertices[29](1)=0; vertices[29](2)=PointFrontBot(2);
vertices[30](0)=0; vertices[30](1)=-PointMidTop(1); vertices[30](2)=0;
vertices[31](0)=PointMidBot(0); vertices[31](1)=0; vertices[31](2)=PointMidBot(2);
vertices[32](0)=PointLeftTransom(0); vertices[32](1)=PointLeftTransom(1); vertices[32](2)=0;
vertices[33](0)=PointFrontTop(0); vertices[33](1)=0; vertices[33](2)=0;
vertices[34](0)=-Lx_domain/2; vertices[34](1)=0; vertices[34](2)=0;
vertices[35](0)=-Lx_domain/2; vertices[35](1)=-Ly_domain/2; vertices[35](2)=0;
vertices[36](0)=a*PointFrontTop(0); vertices[36](1)=-Ly_domain/2; vertices[36](2)=0;
vertices[37](0)=b*PointLeftTransom(0); vertices[37](1)=-Ly_domain/2; vertices[37](2)=0;
vertices[38](0)=Lx_domain/2; vertices[38](1)=-Ly_domain/2; vertices[38](2)=0;
vertices[39](0)=Lx_domain/2; vertices[39](1)=-Ly_domain/4; vertices[39](2)=0;
vertices[40](0)=0; vertices[40](1)=-PointMidTop(1); vertices[40](2)=0;
vertices[41](0)=0; vertices[41](1)=-Ly_domain/2; vertices[41](2)=0;
vertices[42](0)=Lx_domain/2; vertices[42](1)=Ly_domain/2; vertices[42](2)=0;
vertices[43](0)=b*PointRightTransom(0); vertices[43](1)=Ly_domain/2; vertices[43](2)=0;
vertices[44](0)=a*PointFrontTop(0); vertices[44](1)=Ly_domain/2; vertices[44](2)=0;
vertices[45](0)=-Lx_domain/2; vertices[45](1)=Ly_domain/2; vertices[45](2)=0;
vertices[46](0)=0; vertices[46](1)=PointMidTop(1); vertices[46](2)=0;
vertices[47](0)=0; vertices[47](1)=Ly_domain/2; vertices[47](2)=0;
vertices[48](0)=-Lx_domain/2; vertices[48](1)=0; vertices[48](2)=-Lz_domain;
vertices[49](0)=Lx_domain/2; vertices[49](1)=0; vertices[49](2)=-Lz_domain;
vertices[50](0)=-Lx_domain/2; vertices[50](1)=Ly_domain/2; vertices[50](2)=-Lz_domain;
vertices[51](0)=Lx_domain/2; vertices[51](1)=Ly_domain/2; vertices[51](2)=-Lz_domain;
vertices[52](0)=a*PointFrontTop(0); vertices[52](1)=0; vertices[52](2)=-Lz_domain;
vertices[53](0)=0; vertices[53](1)=0; vertices[53](2)=-Lz_domain;
vertices[54](0)=b*(PointRightTransom(0)+PointLeftTransom(0))/2.0; vertices[54](1)=0; vertices[54](2)=-Lz_domain;
vertices[55](0)=a*PointFrontTop(0); vertices[55](1)=Ly_domain/2; vertices[55](2)=-Lz_domain;
vertices[56](0)=0; vertices[56](1)=Ly_domain/2; vertices[56](2)=-Lz_domain;
vertices[57](0)=b*PointRightTransom(0); vertices[57](1)=Ly_domain/2; vertices[57](2)=-Lz_domain;
vertices[58](0)=-Lx_domain/2; vertices[58](1)=0; vertices[58](2)=-Lz_domain;
vertices[59](0)=Lx_domain/2; vertices[59](1)=0; vertices[59](2)=-Lz_domain;
vertices[60](0)=Lx_domain/2; vertices[60](1)=-Ly_domain/2; vertices[60](2)=-Lz_domain;
vertices[61](0)=-Lx_domain/2; vertices[61](1)=-Ly_domain/2; vertices[61](2)=-Lz_domain;
vertices[62](0)=a*PointFrontTop(0); vertices[62](1)=0; vertices[62](2)=-Lz_domain;
vertices[63](0)=0; vertices[63](1)=0; vertices[63](2)=-Lz_domain;
vertices[64](0)=b*(PointRightTransom(0)+PointLeftTransom(0))/2.0; vertices[64](1)=0; vertices[64](2)=-Lz_domain;
vertices[65](0)=b*PointLeftTransom(0); vertices[65](1)=-Ly_domain/2; vertices[65](2)=-Lz_domain;
vertices[66](0)=0; vertices[66](1)=-Ly_domain/2; vertices[66](2)=-Lz_domain;
vertices[67](0)=a*PointFrontTop(0); vertices[67](1)=-Ly_domain/2; vertices[67](2)=-Lz_domain;
vertices[68](0)=-Lx_domain/2; vertices[68](1)=0; vertices[68](2)=0;
vertices[69](0)=-Lx_domain/2; vertices[69](1)=0; vertices[69](2)=-Lz_domain;
vertices[70](0)=-Lx_domain/2; vertices[70](1)=Ly_domain/2; vertices[70](2)=0;
vertices[71](0)=-Lx_domain/2; vertices[71](1)=Ly_domain/2; vertices[71](2)=-Lz_domain;
vertices[72](0)=-Lx_domain/2; vertices[72](1)=0; vertices[72](2)=0;
vertices[73](0)=-Lx_domain/2; vertices[73](1)=0; vertices[73](2)=-Lz_domain;
vertices[74](0)=-Lx_domain/2; vertices[74](1)=-Ly_domain/2; vertices[74](2)=-Lz_domain;
vertices[75](0)=-Lx_domain/2; vertices[75](1)=-Ly_domain/2; vertices[75](2)=0;
vertices[76](0)=Lx_domain/2; vertices[76](1)=-Ly_domain/4; vertices[76](2)=0;
vertices[77](0)=Lx_domain/2; vertices[77](1)=Ly_domain/4; vertices[77](2)=0;
vertices[78](0)=Lx_domain/2; vertices[78](1)=-Ly_domain/2; vertices[78](2)=0;
vertices[79](0)=Lx_domain/2; vertices[79](1)=-Ly_domain/2; vertices[79](2)=-Lz_domain;
vertices[80](0)=Lx_domain/2; vertices[80](1)=0; vertices[80](2)=0;
vertices[81](0)=Lx_domain/2; vertices[81](1)=0; vertices[81](2)=-Lz_domain;
vertices[82](0)=Lx_domain/2; vertices[82](1)=Ly_domain/2; vertices[82](2)=-Lz_domain;
vertices[83](0)=Lx_domain/2; vertices[83](1)=Ly_domain/2; vertices[83](2)=0;
vertices[84](0)=PointRightTransom(0); vertices[84](1)=PointRightTransom(1); vertices[84](2)=0;
vertices[85](0)=Lx_domain/2; vertices[85](1)=Ly_domain/4; vertices[85](2)=0;
vertices[86](0)=Lx_domain/2; vertices[86](1)=0; vertices[86](2)=0;
vertices[87](0)=PointCenterTransom(0); vertices[87](1)=0; vertices[87](2)=PointCenterTransom(2);

cells.resize(35);

cells[0].vertices[0]=0; cells[0].vertices[1]=1; cells[0].vertices[2]=8; cells[0].vertices[3]=5;
cells[1].vertices[0]=5; cells[1].vertices[1] =8; cells[1].vertices[2] =7; cells[1].vertices[3]=9;
cells[2].vertices[0]=9; cells[2].vertices[1] =7; cells[2].vertices[2] =6; cells[2].vertices[3]=4;
cells[3].vertices[0]=4; cells[3].vertices[1] =6; cells[3].vertices[2] =2; cells[3].vertices[3]=3;
cells[4].vertices[0]=10; cells[4].vertices[1] =11; cells[4].vertices[2] =18; cells[4].vertices[3]=15;
cells[5].vertices[0]=15; cells[5].vertices[1] =18; cells[5].vertices[2] =17; cells[5].vertices[3]=19;
cells[6].vertices[0]=19; cells[6].vertices[1] =17; cells[6].vertices[2] =16; cells[6].vertices[3]=14;
cells[7].vertices[0]=14; cells[7].vertices[1] =16; cells[7].vertices[2] =12; cells[7].vertices[3]=13;
cells[8].vertices[0]=21; cells[8].vertices[1] =24; cells[8].vertices[2] =25; cells[8].vertices[3]=23;
cells[9].vertices[0]=24; cells[9].vertices[1] =20; cells[9].vertices[2] =22; cells[9].vertices[3]=25;
cells[10].vertices[0]=27; cells[10].vertices[1] =30; cells[10].vertices[2] =31; cells[10].vertices[3]=29;
cells[11].vertices[0]=30; cells[11].vertices[1] =26; cells[11].vertices[2] =28; cells[11].vertices[3]=31;
cells[12].vertices[0]=34; cells[12].vertices[1] =35; cells[12].vertices[2] =36; cells[12].vertices[3]=33;
cells[13].vertices[0]=33; cells[13].vertices[1] =36; cells[13].vertices[2] =41; cells[13].vertices[3]=40;
cells[14].vertices[0]=40; cells[14].vertices[1] =41; cells[14].vertices[2] =37; cells[14].vertices[3]=32;
cells[15].vertices[0]=32; cells[15].vertices[1] =37; cells[15].vertices[2] =38; cells[15].vertices[3]=39;
cells[16].vertices[0]=85; cells[16].vertices[1] =42; cells[16].vertices[2] =43; cells[16].vertices[3]=84;
cells[17].vertices[0]=84; cells[17].vertices[1] =43; cells[17].vertices[2] =47; cells[17].vertices[3]=46;
cells[18].vertices[0]=46; cells[18].vertices[1] =47; cells[18].vertices[2] =44; cells[18].vertices[3]=33;
cells[19].vertices[0]=33; cells[19].vertices[1] =44; cells[19].vertices[2] =45; cells[19].vertices[3]=34;
cells[20].vertices[0]=49; cells[20].vertices[1] =54; cells[20].vertices[2] =57; cells[20].vertices[3]=51;
cells[21].vertices[0]=54; cells[21].vertices[1] =53; cells[21].vertices[2] =56; cells[21].vertices[3]=57;
cells[22].vertices[0]=53; cells[22].vertices[1] =52; cells[22].vertices[2] =55; cells[22].vertices[3]=56;
cells[23].vertices[0]=52; cells[23].vertices[1] =48; cells[23].vertices[2] =50; cells[23].vertices[3]=55;
cells[24].vertices[0]=58; cells[24].vertices[1] =62; cells[24].vertices[2] =67; cells[24].vertices[3]=61;
cells[25].vertices[0]=62; cells[25].vertices[1] =63; cells[25].vertices[2] =66; cells[25].vertices[3]=67;
cells[26].vertices[0]=63; cells[26].vertices[1] =64; cells[26].vertices[2] =65; cells[26].vertices[3]=66;
cells[27].vertices[0]=64; cells[27].vertices[1] =59; cells[27].vertices[2] =60; cells[27].vertices[3]=65;
cells[28].vertices[0]=69; cells[28].vertices[1] =68; cells[28].vertices[2] =70; cells[28].vertices[3]=71;
cells[29].vertices[0]=72; cells[29].vertices[1] =73; cells[29].vertices[2] =74; cells[29].vertices[3]=75;
cells[30].vertices[0]=81; cells[30].vertices[1] =76; cells[30].vertices[2] =78; cells[30].vertices[3]=79;
cells[31].vertices[0]=77; cells[31].vertices[1] =81; cells[31].vertices[2] =82; cells[31].vertices[3]=83;
cells[32].vertices[0]=87; cells[32].vertices[1] =32; cells[32].vertices[2] =39; cells[32].vertices[3]=86;
cells[33].vertices[0]=86; cells[33].vertices[1] =85; cells[33].vertices[2] =84; cells[33].vertices[3]=87;
cells[34].vertices[0]=81; cells[34].vertices[1] =77; cells[34].vertices[2] =80; cells[34].vertices[3]=76;//*/

cells[0].material_id = 1;
cells[1].material_id = 1;
cells[2].material_id = 1;
cells[3].material_id = 1;
cells[4].material_id = 2;
cells[5].material_id = 2;
cells[6].material_id = 2;
cells[7].material_id = 2;
cells[8].material_id = 3;
cells[9].material_id = 3;
cells[10].material_id = 4;
cells[11].material_id = 4;
cells[12].material_id = 5;
cells[13].material_id = 5;
cells[14].material_id = 5;
cells[15].material_id = 5;
cells[16].material_id = 5;
cells[17].material_id = 5;
cells[18].material_id = 5;
cells[19].material_id = 5;
cells[20].material_id = 7;
cells[21].material_id = 7;
cells[22].material_id = 7;
cells[23].material_id = 7;
cells[24].material_id = 8;
cells[25].material_id = 8;
cells[26].material_id = 8;
cells[27].material_id = 8;
cells[28].material_id = 9;
cells[29].material_id = 10;
cells[30].material_id = 11;
cells[31].material_id = 11;
cells[32].material_id = 5;
cells[33].material_id = 5;
cells[34].material_id = 11;//*/

// waterline (on water) rear left
subcelldata.boundary_lines.push_back (CellData<1>());
subcelldata.boundary_lines.back().vertices[0] = 32; subcelldata.boundary_lines.back().vertices[1] = 40;
subcelldata.boundary_lines.back().material_id = 29;
// waterline (on water) front left
subcelldata.boundary_lines.push_back (CellData<1>());
subcelldata.boundary_lines.back().vertices[0] = 40; subcelldata.boundary_lines.back().vertices[1] = 33;
subcelldata.boundary_lines.back().material_id = 27;
// waterline (on water) front right
subcelldata.boundary_lines.push_back (CellData<1>());
subcelldata.boundary_lines.back().vertices[0] = 33; subcelldata.boundary_lines.back().vertices[1] = 46;
subcelldata.boundary_lines.back().material_id = 26;
// waterline (on water) rear right
subcelldata.boundary_lines.push_back (CellData<1>());
subcelldata.boundary_lines.back().vertices[0] = 46; subcelldata.boundary_lines.back().vertices[1] = 84;
subcelldata.boundary_lines.back().material_id = 28;
// waterline (on boat) right front
subcelldata.boundary_lines.push_back (CellData<1>());
subcelldata.boundary_lines.back().vertices[0] = 20; subcelldata.boundary_lines.back().vertices[1] = 24;
subcelldata.boundary_lines.back().material_id = 21;
// waterline (on boat) right rear
subcelldata.boundary_lines.push_back (CellData<1>());
subcelldata.boundary_lines.back().vertices[0] = 24; subcelldata.boundary_lines.back().vertices[1] = 21;
subcelldata.boundary_lines.back().material_id = 23;
// waterline (on boat) left rear
subcelldata.boundary_lines.push_back (CellData<1>());
subcelldata.boundary_lines.back().vertices[0] = 26; subcelldata.boundary_lines.back().vertices[1] = 30;
subcelldata.boundary_lines.back().material_id = 24;
// waterline (on boat) left front
subcelldata.boundary_lines.push_back (CellData<1>());
subcelldata.boundary_lines.back().vertices[0] = 30; subcelldata.boundary_lines.back().vertices[1] = 27;
subcelldata.boundary_lines.back().material_id = 22;
//front part of the keel (region needed for keel smoothing) left
subcelldata.boundary_lines.push_back (CellData<1>());
subcelldata.boundary_lines.back().vertices[0] = 29; subcelldata.boundary_lines.back().vertices[1] = 27;
subcelldata.boundary_lines.back().material_id = 30;
//front part of the keel (region needed for keel smoothing) right
subcelldata.boundary_lines.push_back (CellData<1>());
subcelldata.boundary_lines.back().vertices[0] = 20; subcelldata.boundary_lines.back().vertices[1] = 22;
subcelldata.boundary_lines.back().material_id = 35;
//central/front part of the keel left
subcelldata.boundary_lines.push_back (CellData<1>());
subcelldata.boundary_lines.back().vertices[0] = 31; subcelldata.boundary_lines.back().vertices[1] = 29;
subcelldata.boundary_lines.back().material_id = 31;
//central/front part of the keel right
subcelldata.boundary_lines.push_back (CellData<1>());
subcelldata.boundary_lines.back().vertices[0] = 22; subcelldata.boundary_lines.back().vertices[1] = 25;
subcelldata.boundary_lines.back().material_id = 36;
//central/rear part of the keel left
subcelldata.boundary_lines.push_back (CellData<1>());
subcelldata.boundary_lines.back().vertices[0] = 28; subcelldata.boundary_lines.back().vertices[1] = 31;
subcelldata.boundary_lines.back().material_id = 31;
//central/rear part of the keel right
subcelldata.boundary_lines.push_back (CellData<1>());
subcelldata.boundary_lines.back().vertices[0] = 25; subcelldata.boundary_lines.back().vertices[1] = 23;
subcelldata.boundary_lines.back().material_id = 36;
//rear part of the keel / transom edge on boat (region needed for keel smoothing) left
subcelldata.boundary_lines.push_back (CellData<1>());
subcelldata.boundary_lines.back().vertices[0] = 26; subcelldata.boundary_lines.back().vertices[1] = 28;
subcelldata.boundary_lines.back().material_id = 32;
//rear part of the keel / transom edge on boat (region needed for keel smoothing) right
subcelldata.boundary_lines.push_back (CellData<1>());
subcelldata.boundary_lines.back().vertices[0] = 23; subcelldata.boundary_lines.back().vertices[1] = 21;
subcelldata.boundary_lines.back().material_id = 37;
//transom edge on water left
subcelldata.boundary_lines.push_back (CellData<1>());
subcelldata.boundary_lines.back().vertices[0] = 32; subcelldata.boundary_lines.back().vertices[1] = 87;
subcelldata.boundary_lines.back().material_id = 40;
//transom edge on water right
subcelldata.boundary_lines.push_back (CellData<1>());
subcelldata.boundary_lines.back().vertices[0] = 84; subcelldata.boundary_lines.back().vertices[1] = 87;
subcelldata.boundary_lines.back().material_id = 41;

}

else
{

Lx_boat = boat_model.boatWetLength;
Lx_domain = Lx_boat*12.0;
Ly_domain = Lx_boat*4.0;
Lz_domain = Lx_boat*2.0;

vertices.resize(84);
      
vertices[0](0)=-Lx_domain/2; vertices[0](1)=-Ly_domain/2 ; vertices[0](2)=0;
vertices[1](0)=-Lx_domain/2; vertices[1](1)=-Ly_domain/2; vertices[1](2)=-Lz_domain;
vertices[2](0)=Lx_domain/2; vertices[2](1)=-Ly_domain/2; vertices[2](2)=-Lz_domain;
vertices[3](0)=Lx_domain/2; vertices[3](1)=-Ly_domain/2; vertices[3](2)=0;
vertices[4](0)=b*PointBackTop(0); vertices[4](1)=-Ly_domain/2; vertices[4](2)=0;
vertices[5](0)=a*PointFrontTop(0); vertices[5](1)=-Ly_domain/2; vertices[5](2)=0;
vertices[6](0)=b*PointBackTop(0); vertices[6](1)=-Ly_domain/2; vertices[6](2)=-Lz_domain;
vertices[7](0)=0; vertices[7](1)=-Ly_domain/2; vertices[7](2)=-Lz_domain;
vertices[8](0)=a*PointFrontTop(0); vertices[8](1)=-Ly_domain/2; vertices[8](2)=-Lz_domain;
vertices[9](0)=0; vertices[9](1)=-Ly_domain/2; vertices[9](2)=0;
vertices[10](0)=Lx_domain/2; vertices[10](1)=Ly_domain/2; vertices[10](2)=0;
vertices[11](0)=Lx_domain/2; vertices[11](1)=Ly_domain/2; vertices[11](2)=-Lz_domain;
vertices[12](0)=-Lx_domain/2; vertices[12](1)=Ly_domain/2; vertices[12](2)=-Lz_domain;
vertices[13](0)=-Lx_domain/2; vertices[13](1)=Ly_domain/2; vertices[13](2)=0;
vertices[14](0)=a*PointFrontTop(0); vertices[14](1)=Ly_domain/2; vertices[14](2)=0;
vertices[15](0)=b*PointBackTop(0); vertices[15](1)=Ly_domain/2; vertices[15](2)=0;
vertices[16](0)=a*PointFrontTop(0); vertices[16](1)=Ly_domain/2; vertices[16](2)=-Lz_domain;
vertices[17](0)=0; vertices[17](1)=Ly_domain/2; vertices[17](2)=-Lz_domain;
vertices[18](0)=b*PointBackTop(0); vertices[18](1)=Ly_domain/2; vertices[18](2)=-Lz_domain;
vertices[19](0)=0; vertices[19](1)=Ly_domain/2; vertices[19](2)=0;
vertices[20](0)=PointFrontTop(0); vertices[20](1)=0; vertices[20](2)=0;
vertices[21](0)=PointBackTop(0); vertices[21](1)=0; vertices[21](2)=0;
vertices[22](0)=PointFrontBot(0); vertices[22](1)=0; vertices[22](2)=PointFrontBot(2);
vertices[23](0)=PointBackBot(0); vertices[23](1)=0; vertices[23](2)=PointBackBot(2);
vertices[24](0)=0; vertices[24](1)=PointMidTop(1); vertices[24](2)=0;
vertices[25](0)=PointMidBot(0); vertices[25](1)=0; vertices[25](2)=PointMidBot(2);
vertices[26](0)=PointBackTop(0); vertices[26](1)=0; vertices[26](2)=0;
vertices[27](0)=PointFrontTop(0); vertices[27](1)=0; vertices[27](2)=0;
vertices[28](0)=PointBackBot(0); vertices[28](1)=0; vertices[28](2)=PointBackBot(2);
vertices[29](0)=PointFrontBot(0); vertices[29](1)=0; vertices[29](2)=PointFrontBot(2);
vertices[30](0)=0; vertices[30](1)=-PointMidTop(1); vertices[30](2)=0;
vertices[31](0)=PointMidBot(0); vertices[31](1)=0; vertices[31](2)=PointMidBot(2);
vertices[32](0)=PointBackTop(0); vertices[32](1)=0; vertices[32](2)=0;
vertices[33](0)=PointFrontTop(0); vertices[33](1)=0; vertices[33](2)=0;
vertices[34](0)=-Lx_domain/2; vertices[34](1)=0; vertices[34](2)=0;
vertices[35](0)=-Lx_domain/2; vertices[35](1)=-Ly_domain/2; vertices[35](2)=0;
vertices[36](0)=a*PointFrontTop(0); vertices[36](1)=-Ly_domain/2; vertices[36](2)=0;
vertices[37](0)=b*PointBackTop(0); vertices[37](1)=-Ly_domain/2; vertices[37](2)=0;
vertices[38](0)=Lx_domain/2; vertices[38](1)=-Ly_domain/2; vertices[38](2)=0;
vertices[39](0)=Lx_domain/2; vertices[39](1)=0; vertices[39](2)=0;
vertices[40](0)=0; vertices[40](1)=-PointMidTop(1); vertices[40](2)=0;
vertices[41](0)=0; vertices[41](1)=-Ly_domain/2; vertices[41](2)=0;
vertices[42](0)=Lx_domain/2; vertices[42](1)=Ly_domain/2; vertices[42](2)=0;
vertices[43](0)=b*PointBackTop(0); vertices[43](1)=Ly_domain/2; vertices[43](2)=0;
vertices[44](0)=a*PointFrontTop(0); vertices[44](1)=Ly_domain/2; vertices[44](2)=0;
vertices[45](0)=-Lx_domain/2; vertices[45](1)=Ly_domain/2; vertices[45](2)=0;
vertices[46](0)=0; vertices[46](1)=PointMidTop(1); vertices[46](2)=0;
vertices[47](0)=0; vertices[47](1)=Ly_domain/2; vertices[47](2)=0;
vertices[48](0)=-Lx_domain/2; vertices[48](1)=0; vertices[48](2)=-Lz_domain;
vertices[49](0)=Lx_domain/2; vertices[49](1)=0; vertices[49](2)=-Lz_domain;
vertices[50](0)=-Lx_domain/2; vertices[50](1)=Ly_domain/2; vertices[50](2)=-Lz_domain;
vertices[51](0)=Lx_domain/2; vertices[51](1)=Ly_domain/2; vertices[51](2)=-Lz_domain;
vertices[52](0)=a*PointFrontTop(0); vertices[52](1)=0; vertices[52](2)=-Lz_domain;
vertices[53](0)=0; vertices[53](1)=0; vertices[53](2)=-Lz_domain;
vertices[54](0)=b*PointBackTop(0); vertices[54](1)=0; vertices[54](2)=-Lz_domain;
vertices[55](0)=a*PointFrontTop(0); vertices[55](1)=Ly_domain/2; vertices[55](2)=-Lz_domain;
vertices[56](0)=0; vertices[56](1)=Ly_domain/2; vertices[56](2)=-Lz_domain;
vertices[57](0)=b*PointBackTop(0); vertices[57](1)=Ly_domain/2; vertices[57](2)=-Lz_domain;
vertices[58](0)=-Lx_domain/2; vertices[58](1)=0; vertices[58](2)=-Lz_domain;
vertices[59](0)=Lx_domain/2; vertices[59](1)=0; vertices[59](2)=-Lz_domain;
vertices[60](0)=Lx_domain/2; vertices[60](1)=-Ly_domain/2; vertices[60](2)=-Lz_domain;
vertices[61](0)=-Lx_domain/2; vertices[61](1)=-Ly_domain/2; vertices[61](2)=-Lz_domain;
vertices[62](0)=a*PointFrontTop(0); vertices[62](1)=0; vertices[62](2)=-Lz_domain;
vertices[63](0)=0; vertices[63](1)=0; vertices[63](2)=-Lz_domain;
vertices[64](0)=b*PointBackTop(0); vertices[64](1)=0; vertices[64](2)=-Lz_domain;
vertices[65](0)=b*PointBackTop(0); vertices[65](1)=-Ly_domain/2; vertices[65](2)=-Lz_domain;
vertices[66](0)=0; vertices[66](1)=-Ly_domain/2; vertices[66](2)=-Lz_domain;
vertices[67](0)=a*PointFrontTop(0); vertices[67](1)=-Ly_domain/2; vertices[67](2)=-Lz_domain;
vertices[68](0)=-Lx_domain/2; vertices[68](1)=0; vertices[68](2)=0;
vertices[69](0)=-Lx_domain/2; vertices[69](1)=0; vertices[69](2)=-Lz_domain;
vertices[70](0)=-Lx_domain/2; vertices[70](1)=Ly_domain/2; vertices[70](2)=0;
vertices[71](0)=-Lx_domain/2; vertices[71](1)=Ly_domain/2; vertices[71](2)=-Lz_domain;
vertices[72](0)=-Lx_domain/2; vertices[72](1)=0; vertices[72](2)=0;
vertices[73](0)=-Lx_domain/2; vertices[73](1)=0; vertices[73](2)=-Lz_domain;
vertices[74](0)=-Lx_domain/2; vertices[74](1)=-Ly_domain/2; vertices[74](2)=-Lz_domain;
vertices[75](0)=-Lx_domain/2; vertices[75](1)=-Ly_domain/2; vertices[75](2)=0;
vertices[76](0)=Lx_domain/2; vertices[76](1)=0; vertices[76](2)=0;
vertices[77](0)=Lx_domain/2; vertices[77](1)=0; vertices[77](2)=-Lz_domain;
vertices[78](0)=Lx_domain/2; vertices[78](1)=-Ly_domain/2; vertices[78](2)=0;
vertices[79](0)=Lx_domain/2; vertices[79](1)=-Ly_domain/2; vertices[79](2)=-Lz_domain;
vertices[80](0)=Lx_domain/2; vertices[80](1)=0; vertices[80](2)=0;
vertices[81](0)=Lx_domain/2; vertices[81](1)=0; vertices[81](2)=-Lz_domain;
vertices[82](0)=Lx_domain/2; vertices[82](1)=Ly_domain/2; vertices[82](2)=-Lz_domain;
vertices[83](0)=Lx_domain/2; vertices[83](1)=Ly_domain/2; vertices[83](2)=0;

cells.resize(32);

cells[0].vertices[0]=0; cells[0].vertices[1]=1; cells[0].vertices[2]=8; cells[0].vertices[3]=5;
cells[1].vertices[0]=5; cells[1].vertices[1] =8; cells[1].vertices[2] =7; cells[1].vertices[3]=9;
cells[2].vertices[0]=9; cells[2].vertices[1] =7; cells[2].vertices[2] =6; cells[2].vertices[3]=4;
cells[3].vertices[0]=4; cells[3].vertices[1] =6; cells[3].vertices[2] =2; cells[3].vertices[3]=3;
cells[4].vertices[0]=10; cells[4].vertices[1] =11; cells[4].vertices[2] =18; cells[4].vertices[3]=15;
cells[5].vertices[0]=15; cells[5].vertices[1] =18; cells[5].vertices[2] =17; cells[5].vertices[3]=19;
cells[6].vertices[0]=19; cells[6].vertices[1] =17; cells[6].vertices[2] =16; cells[6].vertices[3]=14;
cells[7].vertices[0]=14; cells[7].vertices[1] =16; cells[7].vertices[2] =12; cells[7].vertices[3]=13;
cells[8].vertices[0]=21; cells[8].vertices[1] =24; cells[8].vertices[2] =25; cells[8].vertices[3]=23;
cells[9].vertices[0]=24; cells[9].vertices[1] =20; cells[9].vertices[2] =22; cells[9].vertices[3]=25;
cells[10].vertices[0]=27; cells[10].vertices[1] =30; cells[10].vertices[2] =31; cells[10].vertices[3]=29;
cells[11].vertices[0]=30; cells[11].vertices[1] =26; cells[11].vertices[2] =28; cells[11].vertices[3]=31;
cells[12].vertices[0]=34; cells[12].vertices[1] =35; cells[12].vertices[2] =36; cells[12].vertices[3]=33;
cells[13].vertices[0]=33; cells[13].vertices[1] =36; cells[13].vertices[2] =41; cells[13].vertices[3]=40;
cells[14].vertices[0]=40; cells[14].vertices[1] =41; cells[14].vertices[2] =37; cells[14].vertices[3]=32;
cells[15].vertices[0]=32; cells[15].vertices[1] =37; cells[15].vertices[2] =38; cells[15].vertices[3]=39;
cells[16].vertices[0]=39; cells[16].vertices[1] =42; cells[16].vertices[2] =43; cells[16].vertices[3]=32;
cells[17].vertices[0]=32; cells[17].vertices[1] =43; cells[17].vertices[2] =47; cells[17].vertices[3]=46;
cells[18].vertices[0]=46; cells[18].vertices[1] =47; cells[18].vertices[2] =44; cells[18].vertices[3]=33;
cells[19].vertices[0]=33; cells[19].vertices[1] =44; cells[19].vertices[2] =45; cells[19].vertices[3]=34;
cells[20].vertices[0]=49; cells[20].vertices[1] =54; cells[20].vertices[2] =57; cells[20].vertices[3]=51;
cells[21].vertices[0]=54; cells[21].vertices[1] =53; cells[21].vertices[2] =56; cells[21].vertices[3]=57;
cells[22].vertices[0]=53; cells[22].vertices[1] =52; cells[22].vertices[2] =55; cells[22].vertices[3]=56;
cells[23].vertices[0]=52; cells[23].vertices[1] =48; cells[23].vertices[2] =50; cells[23].vertices[3]=55;
cells[24].vertices[0]=58; cells[24].vertices[1] =62; cells[24].vertices[2] =67; cells[24].vertices[3]=61;
cells[25].vertices[0]=62; cells[25].vertices[1] =63; cells[25].vertices[2] =66; cells[25].vertices[3]=67;
cells[26].vertices[0]=63; cells[26].vertices[1] =64; cells[26].vertices[2] =65; cells[26].vertices[3]=66;
cells[27].vertices[0]=64; cells[27].vertices[1] =59; cells[27].vertices[2] =60; cells[27].vertices[3]=65;
cells[28].vertices[0]=69; cells[28].vertices[1] =68; cells[28].vertices[2] =70; cells[28].vertices[3]=71;
cells[29].vertices[0]=72; cells[29].vertices[1] =73; cells[29].vertices[2] =74; cells[29].vertices[3]=75;
cells[30].vertices[0]=77; cells[30].vertices[1] =76; cells[30].vertices[2] =78; cells[30].vertices[3]=79;
cells[31].vertices[0]=80; cells[31].vertices[1] =81; cells[31].vertices[2] =82; cells[31].vertices[3]=83;//*/

cells[0].material_id = 1;
cells[1].material_id = 1;
cells[2].material_id = 1;
cells[3].material_id = 1;
cells[4].material_id = 2;
cells[5].material_id = 2;
cells[6].material_id = 2;
cells[7].material_id = 2;
cells[8].material_id = 3;
cells[9].material_id = 3;
cells[10].material_id = 4;
cells[11].material_id = 4;
cells[12].material_id = 5;
cells[13].material_id = 5;
cells[14].material_id = 5;
cells[15].material_id = 5;
cells[16].material_id = 6;
cells[17].material_id = 6;
cells[18].material_id = 6;
cells[19].material_id = 6;
cells[20].material_id = 7;
cells[21].material_id = 7;
cells[22].material_id = 7;
cells[23].material_id = 7;
cells[24].material_id = 8;
cells[25].material_id = 8;
cells[26].material_id = 8;
cells[27].material_id = 8;
cells[28].material_id = 9;
cells[29].material_id = 10;
cells[30].material_id = 11;
cells[31].material_id = 12;//*/

// waterline (on water) rear left
subcelldata.boundary_lines.push_back (CellData<1>());
subcelldata.boundary_lines.back().vertices[0] = 32; subcelldata.boundary_lines.back().vertices[1] = 40;
subcelldata.boundary_lines.back().material_id = 29;
// waterline (on water) front left
subcelldata.boundary_lines.push_back (CellData<1>());
subcelldata.boundary_lines.back().vertices[0] = 40; subcelldata.boundary_lines.back().vertices[1] = 33;
subcelldata.boundary_lines.back().material_id = 27;
// waterline (on water) front right
subcelldata.boundary_lines.push_back (CellData<1>());
subcelldata.boundary_lines.back().vertices[0] = 33; subcelldata.boundary_lines.back().vertices[1] = 46;
subcelldata.boundary_lines.back().material_id = 26;
// waterline (on water) rear right
subcelldata.boundary_lines.push_back (CellData<1>());
subcelldata.boundary_lines.back().vertices[0] = 46; subcelldata.boundary_lines.back().vertices[1] = 32;
subcelldata.boundary_lines.back().material_id = 28;
// waterline (on boat) right front
subcelldata.boundary_lines.push_back (CellData<1>());
subcelldata.boundary_lines.back().vertices[0] = 20; subcelldata.boundary_lines.back().vertices[1] = 24;
subcelldata.boundary_lines.back().material_id = 21;
// waterline (on boat) right rear
subcelldata.boundary_lines.push_back (CellData<1>());
subcelldata.boundary_lines.back().vertices[0] = 24; subcelldata.boundary_lines.back().vertices[1] = 21;
subcelldata.boundary_lines.back().material_id = 23;
// waterline (on boat) left rear
subcelldata.boundary_lines.push_back (CellData<1>());
subcelldata.boundary_lines.back().vertices[0] = 26; subcelldata.boundary_lines.back().vertices[1] = 30;
subcelldata.boundary_lines.back().material_id = 24;
// waterline (on boat) left front
subcelldata.boundary_lines.push_back (CellData<1>());
subcelldata.boundary_lines.back().vertices[0] = 30; subcelldata.boundary_lines.back().vertices[1] = 27;
subcelldata.boundary_lines.back().material_id = 22;
//front part of the keel (region needed for keel smoothing) left
subcelldata.boundary_lines.push_back (CellData<1>());
subcelldata.boundary_lines.back().vertices[0] = 29; subcelldata.boundary_lines.back().vertices[1] = 27;
subcelldata.boundary_lines.back().material_id = 30;
//front part of the keel (region needed for keel smoothing) right
subcelldata.boundary_lines.push_back (CellData<1>());
subcelldata.boundary_lines.back().vertices[0] = 20; subcelldata.boundary_lines.back().vertices[1] = 22;
subcelldata.boundary_lines.back().material_id = 35;
//central/front part of the keel left
subcelldata.boundary_lines.push_back (CellData<1>());
subcelldata.boundary_lines.back().vertices[0] = 31; subcelldata.boundary_lines.back().vertices[1] = 29;
subcelldata.boundary_lines.back().material_id = 31;
//central/front part of the keel right
subcelldata.boundary_lines.push_back (CellData<1>());
subcelldata.boundary_lines.back().vertices[0] = 22; subcelldata.boundary_lines.back().vertices[1] = 25;
subcelldata.boundary_lines.back().material_id = 36;
//central/rear part of the keel left
subcelldata.boundary_lines.push_back (CellData<1>());
subcelldata.boundary_lines.back().vertices[0] = 28; subcelldata.boundary_lines.back().vertices[1] = 31;
subcelldata.boundary_lines.back().material_id = 31;
//central/rear part of the keel right
subcelldata.boundary_lines.push_back (CellData<1>());
subcelldata.boundary_lines.back().vertices[0] = 25; subcelldata.boundary_lines.back().vertices[1] = 23;
subcelldata.boundary_lines.back().material_id = 36;
//rear part of the keel (region needed for keel smoothing) left
subcelldata.boundary_lines.push_back (CellData<1>());
subcelldata.boundary_lines.back().vertices[0] = 26; subcelldata.boundary_lines.back().vertices[1] = 28;
subcelldata.boundary_lines.back().material_id = 32;
//rear part of the keel (region needed for keel smoothing) right
subcelldata.boundary_lines.push_back (CellData<1>());
subcelldata.boundary_lines.back().vertices[0] = 23; subcelldata.boundary_lines.back().vertices[1] = 21;
subcelldata.boundary_lines.back().material_id = 37;

}


GridTools::delete_unused_vertices (vertices, cells, subcelldata);
GridReordering<2,3>::reorder_cells (cells);

triangulation.create_triangulation_compatibility(vertices, cells, subcelldata );


  std::ofstream logfile("meshResult.inp");
  GridOut grid_out;
  grid_out.write_ucd(triangulation, logfile);
//*/

}

void NumericalTowingTank::generate_double_nodes_set()
{

  ComputationalDomain<3>::generate_double_nodes_set();

  compute_nodes_flags();
  surface_nodes.reinit(dh.n_dofs());
  other_nodes.reinit(dh.n_dofs());
  for (unsigned int i=0; i<dh.n_dofs();++i)
      {
      if (flags[i] & water )
         surface_nodes(i) = 1;
      else
         other_nodes(i) = 1;
      } 
  iges_normals.clear();
  iges_mean_curvatures.clear();
  iges_normals.resize(dh.n_dofs());
  old_iges_normals.resize(dh.n_dofs());
  iges_mean_curvatures.resize(dh.n_dofs());
  edges_tangents.reinit(vector_dh.n_dofs());
  edges_length_ratios.reinit(vector_dh.n_dofs());
  smoothing_curvature_vector.reinit(vector_dh.n_dofs());
  compute_constraints(vector_constraints);
  normals_sparsity_pattern.reinit(vector_dh.n_dofs(),
                                   vector_dh.n_dofs(),
                                   vector_dh.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (vector_dh, normals_sparsity_pattern, vector_constraints);
  normals_sparsity_pattern.compress();
  compute_normals_at_nodes(map_points);
  set_up_smoother();

/*
// Let's put all the pressure quad nodes in memory
  cell_it
  cell = dh.begin_active(),
  endc = dh.end();
  FEValues<2,3> fe_v(*mapping, fe, *quadrature,
	              update_values | update_gradients |
		      update_cell_normal_vectors |
		      update_quadrature_points |
		      update_JxW_values);

   const unsigned int n_q_points = fe_v.n_quadrature_points;
   const unsigned int dofs_per_cell   = fe.dofs_per_cell;
   std::vector<unsigned int> local_dof_indices (dofs_per_cell);


  for (; cell != endc; ++cell)    
      {
      if ((cell->material_id() == wall_sur_ID1 )) // right side of the boat
         {
         fe_v.reinit(cell);
         const std::vector<Point<3> > &node_positions = fe_v.get_quadrature_points();
         const std::vector<Point<dim> > &node_normals = fe_v.get_normal_vectors();
         std::vector<Point<3> > proj_quad_nodes(n_q_points);
         for (unsigned int q=0; q<n_q_points; ++q)
             {
             boat_model.boat_water_line_right->assigned_axis_projection_and_diff_forms(proj_quad_nodes[],
                                                                                   iges_normals[i],
                                                                                   iges_mean_curvatures[i],
                                                                                   vector_support_points[3*i],
                                                                                   node_normals[i]);  // for projection in mesh normal direction
             }
         }
      }     
*/

}

  

void NumericalTowingTank::compute_nodes_flags()
{

  double tol = 1e-8;

  unsigned int free_sur_ID1 = 5;
  unsigned int free_sur_ID2 = 6; 
  unsigned int wall_sur_ID1 = 3;
  unsigned int wall_sur_ID2 = 4; 
  unsigned int inflow_sur_ID1 = 9;
  unsigned int inflow_sur_ID2 = 10;


  cell_it
  vector_cell = vector_dh.begin_active(),
  vector_endc = vector_dh.end();
  
  cell_it
  cell = dh.begin_active(),
  endc = dh.end();
     
  std::vector<unsigned int> dofs(fe.dofs_per_cell);
  std::vector<unsigned int> vector_dofs(vector_fe.dofs_per_cell);

   // mappa che associa ad ogni dof le celle cui esso appartiene
  dof_to_elems.clear();
  
  // mappa che associa ad ogni gradient dof le celle cui esso appartiene
  vector_dof_to_elems.clear();
  
  // vettore che associa ad ogni gradient dof la sua componente
  vector_dof_components.clear();
  vector_dof_components.resize(vector_dh.n_dofs());
  
  // mappa che associa ad ogni cella un set contenente le celle circostanti
  elem_to_surr_elems.clear();
  
  // set che raccoglie i nodi della free surface che stanno sulla barca
  free_surf_and_boat_nodes.clear();
  
  // mappa raccoglie i nodi degli edges della barca
  boat_nodes.clear();  

   
   for (; cell!=endc,vector_cell!=vector_endc; ++cell,++vector_cell)
    {
    Assert(cell->index() == vector_cell->index(), ExcInternalError());
    
    cell->get_dof_indices(dofs);
    for(unsigned int j=0; j<fe.dofs_per_cell; ++j)
        {
	dof_to_elems[dofs[j]].push_back(cell);
	}
    vector_cell->get_dof_indices(vector_dofs);
    for(unsigned int j=0; j<vector_fe.dofs_per_cell; ++j)
        {
	vector_dof_to_elems[vector_dofs[j]].push_back(vector_cell);
	vector_dof_components[vector_dofs[j]] = vector_fe.system_to_component_index(j).first;
	}
    }
  
  // qui viene creata la mappa degli elmenti che circondano ciascun elemento    
  for (cell = dh.begin_active(); cell != endc; ++cell)    
    {
    cell->get_dof_indices(dofs);
    for(unsigned int j=0; j<fe.dofs_per_cell; ++j)
        {
	std::set <unsigned int> duplicates = double_nodes_set[dofs[j]];
	for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
	    {
	    std::vector<cell_it>
	    dof_cell_list = dof_to_elems[*pos];
	    for (unsigned int k=0; k<dof_cell_list.size(); ++k)
                elem_to_surr_elems[cell].insert(dof_cell_list[k]);
	    }
	}
    }


std::vector<unsigned int> face_dofs(fe.dofs_per_face);
for (cell=dh.begin_active(); cell!= endc;++cell)
          {
          if ((cell->material_id() == free_sur_ID1 ||
               cell->material_id() == free_sur_ID2 ||
               cell->material_id() == free_sur_ID3 ))
             {
	     if ( cell->at_boundary() )
                {
                for (unsigned int j = 0; j < GeometryInfo<2>::faces_per_cell; j++)
                    {
                    if (cell->face(j)->boundary_indicator() == free_sur_edge_on_boat_ID )
	               {
                       cell->face(j)->get_dof_indices(face_dofs); 
                       for (unsigned int k=0; k<fe.dofs_per_face; ++k)   
	                   free_surf_and_boat_nodes.insert(face_dofs[k]);
                       }
                    }
                }
             }
          } 


update_support_points();

std::vector<unsigned int> vector_local_dof_indices (vector_fe.dofs_per_cell);
      vector_cell = vector_dh.begin_active();
      vector_endc = vector_dh.end();
      
      for (; vector_cell!=vector_endc; ++vector_cell)
	{//std::cout<<"??1 "<<vector_fe.dofs_per_cell<<std::endl;
	  vector_cell->get_dof_indices(vector_local_dof_indices);
	  if (vector_cell->material_id() == wall_sur_ID1 ||
	      vector_cell->material_id() == wall_sur_ID2 ||
	      vector_cell->material_id() == wall_sur_ID3   )
	    for (unsigned int j=0; j<vector_fe.dofs_per_cell; ++j) 
	      {//std::cout<<"??2"<<std::endl;
		unsigned int id=vector_local_dof_indices[j];
                boat_nodes.insert(id);
                //std::cout<<"??3 "<<id<<" "<<comp_i<<std::endl;
	      }
	}

  edge_cells.clear();
  {
    tria_it
      cell = tria.begin_active(),
      endc = tria.end();

    for (; cell != endc; ++cell)
        if (cell->at_boundary() )
           edge_cells.insert(cell);

  }



  
    //for( std::map<tria_it, tria_it>::iterator it=boat_to_water_edge_cells.begin();
    //	it != boat_to_water_edge_cells.end(); ++it)
    //std::cout << it->first << " -> " << it->second << std::endl;
  //}
//*/


  flags.clear();
  flags.resize(dh.n_dofs());

  vector_flags.clear();
  vector_flags.resize(vector_dh.n_dofs());

  cell_flags.clear();
  cell_flags.resize(tria.n_active_cells());



  unsigned int cell_id=0;
  std::vector<unsigned int> vector_face_dofs(vector_fe.dofs_per_face);
  //std::vector<unsigned int> face_dofs(fe.dofs_per_face);

  vector_cell = vector_dh.begin_active();
  vector_endc = vector_dh.end();
  
  cell = dh.begin_active();
  endc = dh.end();
     
  //std::vector<unsigned int> dofs(fe.dofs_per_cell);
  //std::vector<unsigned int> vector_dofs(vector_fe.dofs_per_cell);

  for (cell = dh.begin_active(), vector_cell = vector_dh.begin_active();
       cell != endc; ++cell, ++vector_cell, ++cell_id)    
    {
    Assert(cell->index() == vector_cell->index(), ExcInternalError());
    cell->get_dof_indices(dofs);
    vector_cell->get_dof_indices(vector_dofs);

				     // left or right
    if (cell->center()(1) > 0)
      {
	for(unsigned int j=0; j<fe.dofs_per_cell; ++j)
	  flags[dofs[j]] |= right_side;
	for(unsigned int j=0; j<vector_fe.dofs_per_cell; ++j)
	  vector_flags[vector_dofs[j]] |= right_side;
      }
    else
      {
	for(unsigned int j=0; j<fe.dofs_per_cell; ++j)
	  flags[dofs[j]] |= left_side;
	for(unsigned int j=0; j<vector_fe.dofs_per_cell; ++j)
	  vector_flags[vector_dofs[j]] |= left_side;
      }
				     // Free surface
    if ((cell->material_id() == free_sur_ID1 ||
	 cell->material_id() == free_sur_ID2 )) 
      {
	for(unsigned int j=0; j<fe.dofs_per_cell; ++j)
	  flags[dofs[j]] |= water;
	for(unsigned int j=0; j<vector_fe.dofs_per_cell; ++j)
	  vector_flags[vector_dofs[j]] |= water;

      }				     // boat surface
    else if ((cell->material_id() == wall_sur_ID1 ||
	      cell->material_id() == wall_sur_ID2 ))
      {
	for(unsigned int j=0; j<fe.dofs_per_cell; ++j)
	  flags[dofs[j]] |= boat;
	for(unsigned int j=0; j<vector_fe.dofs_per_cell; ++j)
	  vector_flags[vector_dofs[j]] |= boat;
      }
    else if ((cell->material_id() == inflow_sur_ID1 ||
	      cell->material_id() == inflow_sur_ID2 ))
      {
	for(unsigned int j=0; j<fe.dofs_per_cell; ++j)
	  flags[dofs[j]] |= inflow;
	for(unsigned int j=0; j<vector_fe.dofs_per_cell; ++j)
	  vector_flags[vector_dofs[j]] |= inflow;
      }
    else
      {
	for(unsigned int j=0; j<fe.dofs_per_cell; ++j)
	  flags[dofs[j]] |= walls;
	for(unsigned int j=0; j<vector_fe.dofs_per_cell; ++j)
	  vector_flags[vector_dofs[j]] |= walls;
      }
    
      
    for(unsigned int f=0; f<GeometryInfo<2>::faces_per_cell; ++f)
      if( cell->face(f)->at_boundary() ) 
	{
	  cell->face(f)->get_dof_indices(face_dofs);
	  vector_cell->face(f)->get_dof_indices(vector_face_dofs);
	  for(unsigned int k=0; k<face_dofs.size(); ++k)
	    flags[face_dofs[k]] |= edge;
	  for(unsigned int k=0; k<vector_face_dofs.size(); ++k)
	    vector_flags[vector_face_dofs[k]] |= edge;
	}
    }
  				   // Now set the relationships...
    for(unsigned int i=0; i<dh.n_dofs(); ++i)
    {
      std::set<unsigned int> doubles = double_nodes_set[i];
      doubles.erase(i);
      for (std::set<unsigned int>::iterator it = doubles.begin() ; it != doubles.end(); it++ )
	{
	  if (flags[*it] & water)
	    {
            flags[i] |= near_water;
            }
	  else if (flags[*it] & boat)
	    {
            flags[i] |= near_boat;
            if(flags[i] & boat)
                {
		flags[i] |= keel;
                }
            }
	  else if (flags[*it] & inflow)
	    {
            flags[i] |= near_inflow;
            }
	  else
	    flags[i] |= near_walls;
	}
    }

  				   // Now set the relationships...
    for(unsigned int i=0; i<vector_dh.n_dofs(); ++i)
    {
      std::set<unsigned int> doubles = vector_double_nodes_set[i];
      doubles.erase(i);
      for (std::set<unsigned int>::iterator it = doubles.begin() ; it != doubles.end(); it++ )
	{
	  if (vector_flags[*it] & water)
	    {
            vector_flags[i] |= near_water;
            }
	  else if (vector_flags[*it] & boat)
	    {
            vector_flags[i] |= near_boat;
            if(vector_flags[i] & boat)
                {
		vector_flags[i] |= keel;
                //z++;
                }
            }
	  else if (vector_flags[*it] & inflow)
	    {
            vector_flags[i] |= near_inflow;
            }
	  else
	    vector_flags[i] |= near_walls;
	}
    }

    for(unsigned int i=0; i<dh.n_dofs(); ++i)
    {
    std::set<unsigned int> doubles = double_nodes_set[i];
    if ((flags[i] & near_boat) &&
        (doubles.size() == 3))
        flags[i] |= keel;
        
    }

    for(unsigned int i=0; i<vector_dh.n_dofs(); ++i)
    {
    std::set<unsigned int> doubles = vector_double_nodes_set[i];
    if ((vector_flags[i] & near_boat) &&
        (doubles.size() == 3))
        vector_flags[i] |= keel;
        
    }

  if (boat_model.is_transom)
     {

     // parameters for rear left water line smoothing
     unsigned int left_boat_transom_point_id = 0;     
     unsigned int right_boat_transom_point_id = 0;     
     unsigned int left_water_transom_point_id = 0;     
     unsigned int right_water_transom_point_id = 0;     

     unsigned int point_id_left = find_point_id(boat_model.PointLeftTransom, ref_points);
     unsigned int point_id_right = find_point_id(boat_model.PointRightTransom, ref_points);
     std::set<unsigned int> duplicates = double_nodes_set[point_id_left]; 
     for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
         {
         //cout<<i<<" mpid"<<*pos<<"  is in?"<<boundary_dofs[i][3*(*pos)]<<endl;
         if ( flags[*pos] & water)
	    {
            left_water_transom_point_id = *pos;
            }
         else
	    {
            left_boat_transom_point_id = *pos;
            }          
         }
     duplicates = double_nodes_set[point_id_right];
     for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
         {
         //cout<<i<<" mpid"<<*pos<<"  is in?"<<boundary_dofs[i][3*(*pos)]<<endl;
         if ( flags[*pos] & water)
	    {
            right_water_transom_point_id = *pos;
            }
         else
	    {
            right_boat_transom_point_id = *pos;
            }          
         } 


     std::vector<unsigned int> face_dofs(fe.dofs_per_face);
     for (cell=dh.begin_active(); cell!= endc;++cell)
         {
	 if ( cell->at_boundary() )
            {
            for (unsigned int j = 0; j < GeometryInfo<2>::faces_per_cell; j++)
                {
                if ((cell->face(j)->boundary_indicator() == 40) ||
                    (cell->face(j)->boundary_indicator() == 41)   )
	           {
                   cell->face(j)->get_dof_indices(face_dofs); 
                   for (unsigned int k=0; k<fe.dofs_per_face; ++k)
                       {
                       if ( (right_water_transom_point_id != face_dofs[k]) &&
                            (left_water_transom_point_id != face_dofs[k]) ) 
                          {  
	                  flags[face_dofs[k]] |= transom_on_water;
                          for (unsigned int p=0;p<3;++p)
                              vector_flags[3*face_dofs[k]+p] |= transom_on_water;
                          }
                       }
                   }
                else if ((cell->face(j)->boundary_indicator() == 37) ||
                         (cell->face(j)->boundary_indicator() == 32)   )
                   {
                   cell->face(j)->get_dof_indices(face_dofs); 
                   for (unsigned int k=0; k<fe.dofs_per_face; ++k)   
                       {
                       if ( (right_boat_transom_point_id != face_dofs[k]) &&
                            (left_boat_transom_point_id != face_dofs[k]) ) 
                          {
	                  flags[face_dofs[k]] |= transom_on_boat;
                          for (unsigned int p=0;p<3;++p)
                              vector_flags[3*face_dofs[k]+p] |= transom_on_boat;
                          }
                       }
                   }
                }
            }

         }

     }


}

void NumericalTowingTank::set_up_smoother()
{
if (surface_smoother)
   {
   surface_smoother->update_reference();
   }
else
   {
   surface_smoother = new SurfaceSmoothing(smoothing_map_points, smoothing_curvature_vector,
   				           vector_dh, StaticMappingQ1<2,3>::mapping);
   }
if (restart_surface_smoother)
   {
   restart_surface_smoother->update_reference();
   }
else
   {
   restart_surface_smoother = new SurfaceSmoothing(smoothing_map_points, smoothing_curvature_vector,
   				           vector_dh, *mapping);
   }
if ( line_smoothers[0] ) 
   {
   update_smoother();
   }
else
   {
   initialize_smoother();
   }

}

void NumericalTowingTank::initialize_smoother()
{  
  



				   // The boundary dofs
  boundary_dofs.resize(7, vector<bool>(vector_dh.n_dofs(), false));
				   // The boundary ids
  boundary_ids.resize(7);
    
  base_points.resize(7);
  moving_points.resize(7);
  base_point_ids.resize(7);
  moving_point_ids.resize(7);

  curves.resize(7);
  on_curve_option.resize(7);
  
  // parameters for front keel smoothing
  base_points[0]   = boat_model.PointFrontBot;
  moving_points[0] = boat_model.PointFrontTop;
  curves[0] = boat_model.equiv_keel_bspline;
  boundary_ids[0] = 30;
  on_curve_option[0] = true;

  // parameters for rear keel/left transom smoothing
  if (boat_model.is_transom)
     {
     base_points[1] = boat_model.PointCenterTransom;
     moving_points[1] = boat_model.PointLeftTransom;
     curves[1] = boat_model.left_transom_bspline;
     boundary_ids[1] = 32;
     on_curve_option[1] = true;
     }
  else
     {
     base_points[1] = boat_model.PointBackBot;
     moving_points[1] = boat_model.PointBackTop;
     curves[1] = boat_model.equiv_keel_bspline;
     boundary_ids[1] = 32;
     on_curve_option[1] = true;
     }

  // parameters for rear keel/right transom smoothing
  if (boat_model.is_transom)
     {
     base_points[2] = boat_model.PointCenterTransom;
     moving_points[2] = boat_model.PointRightTransom;
     curves[2] = boat_model.right_transom_bspline;
     boundary_ids[2] = 37;
     on_curve_option[2] = true;
     }
  else
     {
     base_points[2] = boat_model.PointBackBot;
     moving_points[2] = boat_model.PointBackTop;
     curves[2] = boat_model.equiv_keel_bspline;
     boundary_ids[2] = 32;
     on_curve_option[2] = true;
     }

  // parameters for front right water line smoothing
  base_points[3]   = boat_model.PointMidTop;
  moving_points[3] = boat_model.PointFrontTop;
  curves[3] = boat_model.right_undisturbed_waterline_curve;
  boundary_ids[3] = 26;
  on_curve_option[3] = false;

  // parameters for front left water line smoothing
  base_points[4]   = Point<3>(boat_model.PointMidTop(0),-boat_model.PointMidTop(1),boat_model.PointMidTop(2));
  moving_points[4] = boat_model.PointFrontTop;
  curves[4] = boat_model.left_undisturbed_waterline_curve;
  boundary_ids[4] = 27;
  on_curve_option[4] = false;

  // parameters for rear right water line smoothing
  base_points[5]   = boat_model.PointMidTop;
  moving_points[5] = moving_points[2];
  curves[5] = boat_model.right_undisturbed_waterline_curve;
  boundary_ids[5] = 28;
  on_curve_option[5] = false;

  // parameters for rear left water line smoothing
  base_points[6]   = Point<3>(boat_model.PointMidTop(0),-boat_model.PointMidTop(1),boat_model.PointMidTop(2));
  moving_points[6] = moving_points[1];
  curves[6] = boat_model.left_undisturbed_waterline_curve;
  boundary_ids[6] = 29;
  on_curve_option[6] = false;

  
  for (unsigned int i=0; i<7; ++i) 
      {
      extract_boundary_dofs(boundary_dofs[i], boundary_ids[i], vector_dh);
      base_point_ids[i] = find_point_id(base_points[i], ref_points);
      moving_point_ids[i] = find_point_id(moving_points[i], ref_points);
      std::set<unsigned int> duplicates = double_nodes_set[moving_point_ids[i]]; 
      for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
          {
          //cout<<i<<" mpid"<<*pos<<"  is in?"<<boundary_dofs[i][3*(*pos)]<<endl;
          if (boundary_dofs[i][3*(*pos)] == 1)
	     {
             moving_point_ids[i] = *pos;
             break;
             }
          }
      duplicates = double_nodes_set[base_point_ids[i]]; 
      for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
          {
          //cout<<i<<" bpid "<<*pos<<"  is in? "<<boundary_dofs[i][3*(*pos)]<<endl;
          if (boundary_dofs[i][3*(*pos)] == 1)
	     {
             base_point_ids[i] = *pos;
             break;
             }
          }
      //cout<<"Base point "<<i<<"  "<<base_points[i]<<"  base_point_ids "<<base_point_ids[i]<<"  "<<ref_points[3*base_point_ids[i]]<<endl;
      }
  

  for (unsigned int i=0; i<7; ++i)
      {
      double smoother_tolerance = boat_model.boatWetLength*1e-3;
      line_smoothers[i] = new LineSmoothing(smoothing_map_points,
	  				    curves[i],
                                            vector_dh,
					    boundary_dofs[i],
					    base_point_ids[i],
					    moving_point_ids[i],
                                            smoother_tolerance);
      }

}



void NumericalTowingTank::update_smoother()
{ 


  // parameters for front keel smoothing
  base_points[0]   = boat_model.PointFrontBot;
  moving_points[0] = boat_model.PointFrontTop;
  curves[0] = boat_model.equiv_keel_bspline;
  boundary_ids[0] = 30;
  on_curve_option[0] = true;

  // parameters for rear keel/left transom smoothing
  if (boat_model.is_transom)
     {
     base_points[1] = boat_model.PointCenterTransom;
     moving_points[1] = boat_model.PointLeftTransom;
     curves[1] = boat_model.left_transom_bspline;
     boundary_ids[1] = 32;
     on_curve_option[1] = true;
     }
  else
     {
     base_points[1] = boat_model.PointBackBot;
     moving_points[1] = boat_model.PointBackTop;
     curves[1] = boat_model.equiv_keel_bspline;
     boundary_ids[1] = 32;
     on_curve_option[1] = true;
     }

  // parameters for rear keel/right transom smoothing
  if (boat_model.is_transom)
     {
     base_points[2] = boat_model.PointCenterTransom;
     moving_points[2] = boat_model.PointRightTransom;
     curves[2] = boat_model.right_transom_bspline;
     boundary_ids[2] = 37;
     on_curve_option[2] = true;
     }
  else
     {
     base_points[2] = boat_model.PointBackBot;
     moving_points[2] = boat_model.PointBackTop;
     curves[2] = boat_model.equiv_keel_bspline;
     boundary_ids[2] = 32;
     on_curve_option[2] = true;
     }

  // parameters for front right water line smoothing
  base_points[3]   = boat_model.PointMidTop;
  moving_points[3] = boat_model.PointFrontTop;
  curves[3] = boat_model.right_undisturbed_waterline_curve;
  boundary_ids[3] = 26;
  on_curve_option[3] = false;

  // parameters for front left water line smoothing
  base_points[4]   = Point<3>(boat_model.PointMidTop(0),-boat_model.PointMidTop(1),boat_model.PointMidTop(0));
  moving_points[4] = boat_model.PointFrontTop;
  curves[4] = boat_model.left_undisturbed_waterline_curve;
  boundary_ids[4] = 27;
  on_curve_option[4] = false;

  // parameters for rear right water line smoothing
  base_points[5]   = boat_model.PointMidTop;
  moving_points[5] = moving_points[2];
  curves[5] = boat_model.right_undisturbed_waterline_curve;
  boundary_ids[5] = 28;
  on_curve_option[5] = false;

  // parameters for rear left water line smoothing
  base_points[6]   = Point<3>(boat_model.PointMidTop(0),-boat_model.PointMidTop(1),boat_model.PointMidTop(0));
  moving_points[6] = moving_points[1];
  curves[6] = boat_model.left_undisturbed_waterline_curve;
  boundary_ids[6] = 29;
  on_curve_option[6] = false;
  
  
  for (unsigned int i=0; i<7; ++i) 
      { 
      //cout<<"Base point "<<i<<"  "<<base_points[i]<<endl;
      extract_boundary_dofs(boundary_dofs[i], boundary_ids[i], vector_dh);
      base_point_ids[i] = find_point_id(base_points[i], ref_points);
      moving_point_ids[i] = find_point_id(moving_points[i], ref_points);
      std::set<unsigned int> duplicates = double_nodes_set[moving_point_ids[i]]; 
      for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
          {
          //cout<<i<<" mpid"<<*pos<<"  is in?"<<boundary_dofs[i][3*(*pos)]<<endl;
          if (boundary_dofs[i][3*(*pos)] == 1)
	     {
             moving_point_ids[i] = *pos;
             break;
             }
          }
      duplicates = double_nodes_set[base_point_ids[i]]; 
      for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
          {
          //cout<<i<<" bpid"<<*pos<<"  is in?"<<boundary_dofs[i][3*(*pos)]<<endl;
          if (boundary_dofs[i][3*(*pos)] == 1)
	     {
             base_point_ids[i] = *pos;
             break;
             }
          }
      }
  
  for(unsigned int i=0; i<7; ++i)
     {
     //cout<<"smoother "<<i<<endl;
     line_smoothers[i]->update_reference(base_point_ids[i],moving_point_ids[i]);
     }
     
}



void NumericalTowingTank::perform_line_smoothing(unsigned int num_smoothings)
{

//smoothing is done on a COPY of map_points
smoothing_map_points = map_points;

                              // line smoothing is performed here
for(unsigned int i=0; i<num_smoothings; ++i)
   {
   line_smoothers[i]->smooth(on_curve_option[i]);
   }


map_points = smoothing_map_points;
                              // we update the support points
update_support_points();

                             // the smoothing moved all the right keel nodes: now
                              // we loop over all the left keel nodes to move them
                              // along with their right twins
for (unsigned int k=0; k<num_smoothings;++k)
    for (unsigned int i=0; i<dh.n_dofs(); ++i)
        {
        if ( (boundary_dofs[k][3*i]) )
           {
           std::set<unsigned int> duplicates = vector_double_nodes_set[3*i];
           duplicates.erase(i); 
           for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
                for (unsigned int j=0; j<3; ++j)
	            map_points(*pos+j) +=  vector_support_points[3*i](j) - vector_support_points[*pos](j);
           }
        }


/*
                            // the smoothing moved all the boat water line nodes: now
                             // we loop over all the free surface water line nodes to
                             // move their horizontal componentsalong with their boat twins.
                             // as for the vertical components, the free surface node rules,
                             // so the vertical component of boat node is taken from that
                             // of free surface node
for (unsigned int k=2; k<6;++k)
    for (unsigned int i=0; i<dh.n_dofs(); ++i)
        {
        if ( boundary_dofs[k][3*i])
           {
           std::set<unsigned int> duplicates = vector_double_nodes_set[3*i];
           duplicates.erase(i); 
           for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
               {
               for (unsigned int j=0; j<2; ++j)
	           map_points(*pos+j) +=  vector_support_points[3*i](j) - vector_support_points[*pos](j);
               map_points(3*i+2) += vector_support_points[*pos](2) - vector_support_points[3*i](2);
               }
//           duplicates = vector_double_nodes_set[3*i];
//           for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
//               {
//               Point<3> p = ref_points[3*i];
//               for (unsigned int j=0; j<3; ++j)
//	           p(j)+= map_points(*pos+j);
//               cout<<3*i<<" ("<<*pos<<")   p("<<p<<")"<<"  ref("<<ref_points[3*i]<<")"<<endl;
//               }
           }
        }
*/

}



void NumericalTowingTank::perform_surface_projection()
{ 
                             // we first update the vector of mesh normals at nodes
compute_normals_at_nodes(map_points);
                             // we update the vector with the support points
update_support_points();

std::pair<double,double> params;
                              // we move all nodes of the right side of
                              // the boat surface using the specified projection
Point<3> proj_node;
for (unsigned int i=0; i<dh.n_dofs(); ++i)
    {
    if ( (flags[i] & boat) &&
         (flags[i] & right_side) &&
         ((flags[i] & edge)== 0) )
       {
       boat_model.boat_water_line_right->assigned_axis_projection_and_diff_forms(proj_node,
                                                                                 iges_normals[i],
                                                                                 iges_mean_curvatures[i],
                                                                                 vector_support_points[3*i],
                                                                                 node_normals[i]);  // for projection in mesh normal direction
       for (unsigned int j=0; j<3; ++j)
           map_points(3*i+j) += proj_node(j) - vector_support_points[3*i](j);
       }
    }

                              // we move all nodes of the left side of
                              // the boat surface using the specified projection
for (unsigned int i=0; i<dh.n_dofs(); ++i)
    {
    if ( (flags[i] & boat) &&
         (flags[i] & left_side) &&
         ((flags[i] & edge)== 0) )
       {
       boat_model.boat_water_line_left->assigned_axis_projection_and_diff_forms(proj_node,
                                                                                iges_normals[i],
                                                                                iges_mean_curvatures[i],
                                                                                vector_support_points[3*i],
                                                                                node_normals[i]);  // for projection in mesh normal direction
       //iges_normals[i]*=-1.0; // reflected shape has wrong orientation! we should change it instead of acting here
       //iges_mean_curvatures[i]*=-1.0; 
       for (unsigned int j=0; j<3; ++j)
           map_points(3*i+j) += proj_node(j) - vector_support_points[3*i](j);
       }
    }

/*
                              // we also need to compute the normals and curvatures on the
                              // keel nodes (right and left side)
for (unsigned int i=0; i<dh.n_dofs(); ++i)
    {
    if ( (flags[i] & keel) &&
         (flags[i] & right_side) )
       {
       boat_model.boat_surface_right->normal_projection_and_diff_forms(proj_node,
                                                                       iges_normals[i],
                                                                       iges_mean_curvatures[i],
                                                                       vector_support_points[3*i]);  // for projection in surface normal direction
       }
    }
for (unsigned int i=0; i<dh.n_dofs(); ++i)
    {
    if ( (flags[i] & keel) &&
         (flags[i] & left_side) )
       {
       boat_model.boat_surface_left->normal_projection_and_diff_forms(proj_node,
                                                                      iges_normals[i],
                                                                      iges_mean_curvatures[i],
                                                                      vector_support_points[3*i]);  // for projection in surface normal direction
       //iges_normals[i]*=-1.0; // reflected shape has wrong orientation! we should change it instead of acting here
       //iges_mean_curvatures[i]*=-1.0; 
       //cout<<"node ("<<vector_support_points[3*i]<<") proj("<<proj_node<<endl;
       //cout<<"normal ("<<iges_normals[i]<<") curvature "<<iges_mean_curvatures[i]<<endl;
       }
    }
*/
}


void NumericalTowingTank::perform_water_line_nodes_projection()
{ 
                             // we update the vector with the support points
update_support_points();

                              // we move all nodes of the right side of
                              // the water line on the boat surface
Point<3> proj_node;
for (unsigned int i=0; i<dh.n_dofs(); ++i)
    {
    if ( (flags[i] & water) &&
         (flags[i] & near_boat) &&
         (flags[i] & right_side) &&
         !(flags[i] & transom_on_water) &&
         (moving_point_ids[3] != i) &&
         (moving_point_ids[4] != i) &&
         (moving_point_ids[5] != i) &&
         (moving_point_ids[6] != i) ) // to avoid the bow and stern node
       {
       boat_model.boat_water_line_right->axis_projection_and_diff_forms(proj_node,
                                                                        iges_normals[i],
                                                                        iges_mean_curvatures[i],
                                                                        vector_support_points[3*i]);  // y axis projection
       //cout<<i<<" (rw) -->  dist "<<proj_node.distance(vector_support_points[3*i])<<" ("<<proj_node<<") Vs ("<<vector_support_points[3*i]<<")"<<endl;
       // " + n("<<iges_normals[i]<<")"<<endl;
       std::set<unsigned int> duplicates = vector_double_nodes_set[3*i];
       for (std::set<unsigned int>::iterator pos=duplicates.begin(); pos!=duplicates.end();pos++)
           for (unsigned int j=0; j<3; ++j)
               map_points(*pos+j) += proj_node(j) - vector_support_points[*pos](j);
       }
    }

                              // we move all nodes of the left side of
                              // the water line on the boat surface
for (unsigned int i=0; i<dh.n_dofs(); ++i)
    {
    if ( (flags[i] & water) &&
         (flags[i] & near_boat) &&
         (flags[i] & left_side) &&
         !(flags[i] & transom_on_water) &&
         (moving_point_ids[3] != i) &&
         (moving_point_ids[4] != i) &&
         (moving_point_ids[5] != i) &&
         (moving_point_ids[6] != i) ) // to avoid the bow and stern node
       {
       boat_model.boat_water_line_left->axis_projection_and_diff_forms(proj_node,
                                                                       iges_normals[i],
                                                                       iges_mean_curvatures[i],
                                                                       vector_support_points[3*i]);  // y axis direction projection
       //cout<<i<<" (lw) -->  dist "<<proj_node.distance(vector_support_points[3*i])<<" ("<<proj_node<<") Vs ("<<vector_support_points[3*i]<<")"<<endl;
       // " + n("<<iges_normals[i]<<")"<<endl;
       std::set<unsigned int> duplicates = vector_double_nodes_set[3*i];
       for (std::set<unsigned int>::iterator pos=duplicates.begin(); pos!=duplicates.end();pos++)
           for (unsigned int j=0; j<3; ++j)
               map_points(*pos+j) += proj_node(j) - vector_support_points[*pos](j);       }
    }

}


void NumericalTowingTank::evaluate_ref_surf_distances(Vector <double> &distances,
                                                      const bool only_surf_smoothing)
{

                             // we update the vector with the support points
update_support_points();
                             // we now update the vector of mesh normals at nodes
                             // using the current geometry
compute_normals_at_nodes(map_points);

//we work on a COPY of map_points
smoothing_map_points = map_points;


//we enforce contraint on the new geometry assigned by the DAE solver
vector_constraints.distribute(smoothing_map_points);

//////////////////////////////////////////
//this takes care of the bow and stern nodes 
if (only_surf_smoothing == false)
for (unsigned int k=3; k<7; ++k)
    { 
    unsigned int i = moving_point_ids[k];
    //cout<<k<<" "<<i<<" "<<support_points[i]<<endl;
         {
         Point <3> dP0 = support_points[i];
         Point <3> dP; 
         				   //this is the horizontal plane
         Handle(Geom_Plane) horPlane = new Geom_Plane(0.,0.,1.,-dP0(2));
         Handle(Geom_Curve) curve;
         if (boat_model.is_transom)
            {
            if (k==3 || k==4)
               curve = boat_model.equiv_keel_bspline;
            else if (k == 6)
               curve = boat_model.left_transom_bspline;
            else
               curve = boat_model.right_transom_bspline;
            }
         else
            {
            curve = boat_model.equiv_keel_bspline;
            }
         GeomAPI_IntCS Intersector(curve, horPlane);
         int npoints = Intersector.NbPoints();
         AssertThrow((npoints != 0), ExcMessage("Keel or transom curve is not intersecting with horizontal plane!"));
         //cout<<"Number of intersections: "<<npoints<<endl;
         double minDistance=1e7;
         gp_Pnt P;
         gp_Vec V1;
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
          // here temporarily for kcs hull tests
            if (minDistance > 0.5*boat_model.boatWetLength)
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
         
         std::set<unsigned int> duplicates = double_nodes_set[i];
          //duplicates.erase(i); 
          for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
              {
	      smoothing_map_points(3*(*pos)) = dP(0)-ref_points[3*(*pos)](0);
              smoothing_map_points(3*(*pos)+1) = dP(1)-ref_points[3*(*pos)](1);
              smoothing_map_points(3*(*pos)+2) = dP(2)-ref_points[3*(*pos)](2);
              edges_tangents[3*(*pos)] = V1.X();
              edges_tangents[3*(*pos)+1] = V1.Y();
              edges_tangents[3*(*pos)+2] = V1.Z();
              //cout<<*pos<<" "<<i<<" "<<smoothing_map_points(3*i)<<" vs "<<map_points(3*i)<<"  ("<<vector_support_points[3*i]<<")"<<endl;
              }
         
         }              
    }


// this cycle hooks the boat and far field double nodes
// to their water twins that have been moved
if (only_surf_smoothing == false)
for (unsigned int i=0; i<vector_dh.n_dofs(); ++i)
    {
    if ( (vector_flags[i] & water) &&
         (vector_flags[i] & edge) )
       {
       std::set<unsigned int> duplicates = vector_double_nodes_set[i];
       duplicates.erase(i); 
       for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
           {
           smoothing_map_points(*pos) = smoothing_map_points(i);
           }
       }
    }
             
//////////////////////////////////////////

                             // line smoothing on keel/transom is performed here and modifies smoothing_map_points
                              // it ONLY moves nodes on the LEFT side of the keel
for(unsigned int i=0; i<3; ++i)
   {
   line_smoothers[i]->smooth(on_curve_option[i]);
   line_smoothers[i]->get_curve_tangent_vectors_at_smoothing_dofs(edges_tangents);
   line_smoothers[i]->get_curve_length_ratios_at_smoothing_dofs(edges_length_ratios);

   for (unsigned int j=0;j<vector_dh.n_dofs();++j)
       if (boundary_dofs[i][j])
          {
          std::set<unsigned int> duplicates = vector_double_nodes_set[j];
          //duplicates.erase(j); 
          for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
              {
              //cout<<*pos<<" "<<smoothing_map_points(*pos)<<" "<<j<<" "<<smoothing_map_points(j)<<endl;
	      smoothing_map_points(*pos) = smoothing_map_points(j);
              edges_tangents(*pos) = edges_tangents(j);
              edges_length_ratios(*pos) = edges_length_ratios(j);
              }
          
          }
   }

                              // line smoothing on water_lines is performed here and modifies smoothing_map_points
                              // it ONLY moves nodes on the water that are also near_boat
if (only_surf_smoothing == false)
for(unsigned int i=3; i<7; ++i)
   {
   line_smoothers[i]->smooth(on_curve_option[i]);
   line_smoothers[i]->get_curve_tangent_vectors_at_smoothing_dofs(edges_tangents);
   line_smoothers[i]->get_curve_length_ratios_at_smoothing_dofs(edges_length_ratios);
   }


                              // we move all nodes of the right side of
                              // the water line on the boat surface
Point<3> proj_node;
if (only_surf_smoothing == false)
for (unsigned int i=0; i<dh.n_dofs(); ++i)
    {
    if ( (flags[i] & water) &&
         (flags[i] & near_boat) &&
         (flags[i] & right_side) &&
         !(flags[i] & transom_on_water) &&
         (moving_point_ids[3] != i) &&
         (moving_point_ids[4] != i) &&
         (moving_point_ids[5] != i) &&
         (moving_point_ids[6] != i)  ) // to avoid the bow and stern node
       {
       Point<3> intermadiate_point_pos(ref_points[3*i](0)+smoothing_map_points(3*i),
                                       ref_points[3*i](1)+smoothing_map_points(3*i+1),
                                       ref_points[3*i](2)+smoothing_map_points(3*i+2));

       Point<3> direction(iges_normals[i](0),iges_normals[i](1),0.0);
       if (direction.square() < 1e-3)
          {
          std::set<unsigned int> duplicates = double_nodes_set[i];
          duplicates.erase(i);
          unsigned int j = *(duplicates.begin());
          direction(0) = node_normals[j](0);
          direction(1) = node_normals[j](1);
          direction(2) = 0.0;
          }        
       //cout<<i<<"(rw) -->  ("<<intermadiate_point_pos<<") and ("<<direction<<")"<<endl;
       boat_model.boat_water_line_right->assigned_axis_projection_and_diff_forms(proj_node,
                                                                                 iges_normals[i],
                                                                                 iges_mean_curvatures[i],
                                                                                 intermadiate_point_pos,
                                                                                 direction);  // hor normal dir projection
       //cout<<i<<"(rw) -->  ("<<proj_node<<") Vs ("<<vector_support_points[3*i]<<") + n("<<iges_normals[i]<<")"<<endl;
       //cout<<i<<"(rw) -->  ("<<proj_node<<") Vs ("<<intermediate_point_pos<<")"<<endl;
       for (unsigned int j=0; j<2; ++j)
           smoothing_map_points(3*i+j) = proj_node(j) - ref_points[3*i](j);
       // we're doing this thing on the water side, but the iges_normal and iges_mean curvature belong to the boat side
       std::set<unsigned int> duplicates = double_nodes_set[i];
       for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
           {
	   iges_normals[*pos] = iges_normals[i];
           iges_mean_curvatures[*pos] = iges_mean_curvatures[i];
           }
       //iges_normals[i] = Point<3>(0.0,0.0,0.0);
       //iges_mean_curvatures[i] = 0;
       }
    }

                              // we move all nodes of the left side of
                              // the water line on the boat surface
if (only_surf_smoothing == false)
for (unsigned int i=0; i<dh.n_dofs(); ++i)
    {
    if ( (flags[i] & water) &&
         (flags[i] & near_boat) &&
         (flags[i] & left_side)  &&
         !(flags[i] & transom_on_water) &&
         (moving_point_ids[3] != i) &&
         (moving_point_ids[4] != i) &&
         (moving_point_ids[5] != i) &&
         (moving_point_ids[6] != i) ) // to avoid the bow and stern node
       {
       Point<3> intermadiate_point_pos(ref_points[3*i](0)+smoothing_map_points(3*i),
                                       ref_points[3*i](1)+smoothing_map_points(3*i+1),
                                       ref_points[3*i](2)+smoothing_map_points(3*i+2));
       Point<3> direction(iges_normals[i](0),iges_normals[i](1),0.0);
       if (direction.square() < 1e-3)
          {
          std::set<unsigned int> duplicates = double_nodes_set[i];
          duplicates.erase(i);
          unsigned int j = *(duplicates.begin());
          direction(0) = node_normals[j](0);
          direction(1) = node_normals[j](1);
          direction(2) = 0.0;
          } 
       boat_model.boat_water_line_left->assigned_axis_projection_and_diff_forms(proj_node,
                                                                                 iges_normals[i],
                                                                                 iges_mean_curvatures[i],
                                                                                 intermadiate_point_pos,
                                                                                 direction);  // hor normal dir projection
       //cout<<i<<"(lw) -->  ("<<proj_node<<") Vs ("<<vector_support_points[3*i]<<") + n("<<iges_normals[i]<<")"<<endl;
       for (unsigned int j=0; j<2; ++j)
           smoothing_map_points(3*i+j) = proj_node(j) - ref_points[3*i](j);
       // we're doing this thing on the water side, but the iges_normal and iges_mean curvature belong to the boat side
       std::set<unsigned int> duplicates = double_nodes_set[i];
       for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
           {
	   iges_normals[*pos] = iges_normals[i];
           iges_mean_curvatures[*pos] = iges_mean_curvatures[i];
           }
       //iges_normals[i] = Point<3>(0.0,0.0,0.0);
       //iges_mean_curvatures[i] = 0;
       }
    }

// this cycle hooks the boat and far field double nodes
// to their water twins that have been moved
if (only_surf_smoothing == false)
for (unsigned int i=0; i<vector_dh.n_dofs(); ++i)
    {
    if ( (vector_flags[i] & water) &&
         (vector_flags[i] & edge) )
       {
       std::set<unsigned int> duplicates = vector_double_nodes_set[i];
       duplicates.erase(i); 
       for (std::set<unsigned int>::iterator pos = duplicates.begin(); pos !=duplicates.end(); pos++)
           {
           smoothing_map_points(*pos) = smoothing_map_points(i);
           }
       }
    }
// the line treatment (smoothing and projection) part is finished here


                              // surface treatment starts here: before doing the smoothing
                              // we must project all nodes on the boat surface to obtain
                              // the normals and curvatures of the surface

                              // we move all nodes of the right side of
                              // the boat surface using the specified projection
for (unsigned int i=0; i<dh.n_dofs(); ++i)
    {
    if ( (flags[i] & boat) &&
         (flags[i] & right_side) &&
         ((flags[i] & edge)== 0))
       {
       Point<3> intermediate_point_pos = ref_points[3*i] +
                                         Point<3>(smoothing_map_points(3*i),smoothing_map_points(3*i+1),smoothing_map_points(3*i+2));
       bool succeed = 
       boat_model.boat_water_line_right->assigned_axis_projection_and_diff_forms(proj_node,
                                                                                 iges_normals[i],
                                                                                 iges_mean_curvatures[i],
                                                                                 intermediate_point_pos,
                                                                                 node_normals[i]);  // for projection in mesh normal direction
       if (succeed == false)
          boat_model.boat_surface_right->normal_projection_and_diff_forms(proj_node,
                                                                          iges_normals[i],
                                                                          iges_mean_curvatures[i],
			                                                  intermediate_point_pos);  // for projection in normal direction

       for (unsigned int j=0; j<3; ++j)
           smoothing_map_points(3*i+j) = proj_node(j) - ref_points[3*i](j);
       }
    }
                              // we move all nodes of the left side of
                              // the boat surface using the specified projection
for (unsigned int i=0; i<dh.n_dofs(); ++i)
    {
    if ( (flags[i] & boat) &&
         (flags[i] & left_side) &&
         ((flags[i] & edge)== 0))
       {
       Point<3> intermediate_point_pos = ref_points[3*i] +
                                         Point<3>(smoothing_map_points(3*i),smoothing_map_points(3*i+1),smoothing_map_points(3*i+2));
       bool succeed = 
       boat_model.boat_water_line_left->assigned_axis_projection_and_diff_forms(proj_node,
                                                                                iges_normals[i],
                                                                                iges_mean_curvatures[i],
                                                                                intermediate_point_pos,
                                                                                node_normals[i]);  // for projection in mesh normal direction
       if (succeed == false)
          boat_model.boat_surface_left->normal_projection_and_diff_forms(proj_node,
                                                                          iges_normals[i],
                                                                          iges_mean_curvatures[i],
			                                                  intermediate_point_pos);  // for projection in normal direction

       for (unsigned int j=0; j<3; ++j)
           smoothing_map_points(3*i+j) = proj_node(j) - ref_points[3*i](j);
       }
   }


                              // with all the normals and curvatures we assemble the vector
                              // used by the smoothing
for (unsigned int i=0; i<vector_dh.n_dofs()/3;++i) 
    {
    if ((boundary_dofs[0][3*i] == false) &&
	(boundary_dofs[1][3*i] == false) && 
        (boundary_dofs[2][3*i] == false) &&
	(boundary_dofs[3][3*i] == false) &&
	(boundary_dofs[4][3*i] == false) &&
	(boundary_dofs[5][3*i] == false) &&
        (boundary_dofs[6][3*i] == false))
        {   
       	for (unsigned int j=0; j<3; ++j) 
       	    {
       	    smoothing_curvature_vector(3*i+j) = -iges_normals[i][j]*(iges_mean_curvatures[i]);
            //cout<<"smooth("<<3*i+j<<") = "<<smoothing_curvature_vector(3*i+j)<<endl;
       	    }
        }
     }

// we perform the surface smoothing
surface_smoother->update_reference();
surface_smoother->smooth();

compute_normals_at_nodes(smoothing_map_points);


// this cycle is needed because surface smoothing
// on free surface must not be effective
for (unsigned int i=0; i<dh.n_dofs(); ++i)
       {
       if ( (flags[i] & water) &&
            ((flags[i] & edge)== 0) )
       {
       smoothing_map_points(3*i) = map_points(3*i);
       smoothing_map_points(3*i+1) = map_points(3*i+1);
       smoothing_map_points(3*i+2) = map_points(3*i+2);
       }
    }

// here we decide if free surface smoothing is active also on boat
// nodes
//   for (unsigned int i=0; i<dh.n_dofs(); ++i)
//       {
//       if ( (flags[i] & boat) &&
//          ((flags[i] & edge)== 0) )
//          {
//          smoothing_map_points(3*i) = map_points(3*i);
//          smoothing_map_points(3*i+1) = map_points(3*i+1);
//          smoothing_map_points(3*i+2) = map_points(3*i+2);
//          }
//       }
                              // we move all nodes of the right side of
                              // the boat surface using the specified projection
for (unsigned int i=0; i<dh.n_dofs(); ++i)
    {
    if ( (flags[i] & boat) &&
         (flags[i] & right_side) &&
         ((flags[i] & edge)== 0))
       {
       Point<3> intermediate_point_pos = ref_points[3*i] +
                                         Point<3>(smoothing_map_points(3*i),smoothing_map_points(3*i+1),smoothing_map_points(3*i+2));
       bool succeed = 
       boat_model.boat_water_line_right->assigned_axis_projection_and_diff_forms(proj_node,
                                                                                 iges_normals[i],
                                                                                 iges_mean_curvatures[i],
                                                                                 intermediate_point_pos,
                                                                                 node_normals[i]);  // for projection in mesh normal direction
       if (succeed == false)
          boat_model.boat_surface_right->normal_projection_and_diff_forms(proj_node,
                                                                          iges_normals[i],
                                                                          iges_mean_curvatures[i],
			                                                  intermediate_point_pos);  // for projection in normal direction
       for (unsigned int j=0; j<3; ++j)
           smoothing_map_points(3*i+j) = proj_node(j) - ref_points[3*i](j);
       }
    }
                              // we move all nodes of the left side of
                              // the boat surface using the specified projection
for (unsigned int i=0; i<dh.n_dofs(); ++i)
    {
    if ( (flags[i] & boat) &&
         (flags[i] & left_side) &&
         ((flags[i] & edge)== 0))
       {
       Point<3> intermediate_point_pos = ref_points[3*i] +
                                         Point<3>(smoothing_map_points(3*i),smoothing_map_points(3*i+1),smoothing_map_points(3*i+2));
       bool succeed = 
       boat_model.boat_water_line_left->assigned_axis_projection_and_diff_forms(proj_node,
                                                                                iges_normals[i],
                                                                                iges_mean_curvatures[i],
                                                                                intermediate_point_pos,
                                                                                node_normals[i]);  // for projection in mesh normal direction
       if (succeed == false)
          boat_model.boat_surface_left->normal_projection_and_diff_forms(proj_node,
                                                                          iges_normals[i],
                                                                          iges_mean_curvatures[i],
			                                                  intermediate_point_pos);  // for projection in normal direction

       for (unsigned int j=0; j<3; ++j)
           smoothing_map_points(3*i+j) = proj_node(j) - ref_points[3*i](j);
       }
   }
                           // we also need to compute the normals and curvatures on the
                           // keel nodes (right and left side)
for (unsigned int i=0; i<dh.n_dofs(); ++i)
    {
    if ( (flags[i] & keel) &&
         (flags[i] & right_side) )
       {
       Point<3> intermediate_point_pos = ref_points[3*i] +
                                         Point<3>(smoothing_map_points(3*i),smoothing_map_points(3*i+1),smoothing_map_points(3*i+2));
       boat_model.boat_surface_right->normal_projection_and_diff_forms(proj_node,
                                                                       iges_normals[i],
                                                                       iges_mean_curvatures[i],
                                                                       intermediate_point_pos);  // for projection in surface normal direction
       }
    }

for (unsigned int i=0; i<dh.n_dofs(); ++i)
    {
    if ( (flags[i] & keel) &&
         (flags[i] & left_side) )
       {
       Point<3> intermediate_point_pos = ref_points[3*i] +
                                         Point<3>(smoothing_map_points(3*i),smoothing_map_points(3*i+1),smoothing_map_points(3*i+2));
       boat_model.boat_surface_left->normal_projection_and_diff_forms(proj_node,
                                                                      iges_normals[i],
                                                                      iges_mean_curvatures[i],
                                                                      intermediate_point_pos);  // for projection in surface normal direction
       }
    }
//*/
vector_constraints.distribute(smoothing_map_points);
distances = smoothing_map_points;


distances*=-1;

distances.add(map_points);

compute_normals_at_nodes(smoothing_map_points);
old_iges_normals = iges_normals;
}


void NumericalTowingTank::perform_smoothing(bool full_treatment, const double blend_factor)
{  
      
      for(unsigned int i=0; i<vector_dh.n_dofs()/3;++i) 
       	{
          if( (boundary_dofs[0][3*i] == false) &&
	      (boundary_dofs[1][3*i] == false) && 
              (boundary_dofs[2][3*i] == false) &&
	      (boundary_dofs[3][3*i] == false) &&
	      (boundary_dofs[4][3*i] == false) &&
	      (boundary_dofs[5][3*i] == false) &&
              (boundary_dofs[6][3*i] == false) )
          {   
       	  for(unsigned int j=0; j<3; ++j) 
       	    {
       	    smoothing_curvature_vector(3*i+j) = -iges_normals[i][j]*(iges_mean_curvatures[i]);
       	    }
          }
       	}
      //smoothing is done on a COPY of map_points
      smoothing_map_points = map_points;
      surface_smoother->smooth();
      // if we are on the free surface, only horizontal components of map_points must be changed
      // according to the smoothing (the third one will be changed according to the free surface differential
      // equation). otherwise, all three components are updated 

      if (full_treatment)
         {
         for (unsigned int i=0; i<vector_dh.n_dofs()/3;++i) 
       	     {
             if ((flags[i] & water) == 0)
                {
                map_points(3*i) = smoothing_map_points(3*i);
                map_points(3*i+1) = smoothing_map_points(3*i+1);
                map_points(3*i+2) = smoothing_map_points(3*i+2);  
                }
             }
          }
        else
          {
           for (unsigned int i=0; i<vector_dh.n_dofs()/3;++i) 
       	     {
             if ( (flags[i] & water) && 
                  ((flags[i] & edge) == 0) )
                {
                map_points(3*i) = old_map_points(3*i) + blend_factor*(smoothing_map_points(3*i)-old_map_points(3*i));
                map_points(3*i+1) = old_map_points(3*i+1) + blend_factor*(smoothing_map_points(3*i+1)-old_map_points(3*i+1));           
                }
             }
          }



}


void NumericalTowingTank::compute_normals_at_nodes(Vector<double> &map_points_used)
{  
   
   vector_normals_matrix.reinit (normals_sparsity_pattern);
   vector_normals_rhs.reinit(vector_dh.n_dofs());
   vector_normals_solution.reinit(vector_dh.n_dofs());

   MappingQEulerian<2, Vector<double>, 3> mappingg(mapping_degree, map_points_used, vector_dh);

   FEValues<2,3> vector_fe_v(mappingg, vector_fe, *quadrature,
			           update_values | update_cell_normal_vectors |  
			           update_JxW_values);

   const unsigned int vector_n_q_points = vector_fe_v.n_quadrature_points;
   const unsigned int   vector_dofs_per_cell   = vector_fe.dofs_per_cell;
   std::vector<unsigned int> vector_local_dof_indices (vector_dofs_per_cell);

   FullMatrix<double>   local_normals_matrix (vector_dofs_per_cell, vector_dofs_per_cell);
   Vector<double>       local_normals_rhs (vector_dofs_per_cell);

   cell_it
   vector_cell = vector_dh.begin_active(),
   vector_endc = vector_dh.end();
            
   for (; vector_cell!=vector_endc; ++vector_cell)
     {
       vector_fe_v.reinit (vector_cell);
       local_normals_matrix = 0;
       local_normals_rhs = 0;
       const std::vector<Point<3> > &vector_node_normals = vector_fe_v.get_normal_vectors();
       unsigned int comp_i, comp_j;
       
       for (unsigned int q=0; q<vector_n_q_points; ++q)
	 for (unsigned int i=0; i<vector_dofs_per_cell; ++i)
	   {
	     comp_i = vector_fe.system_to_component_index(i).first;
	     for (unsigned int j=0; j<vector_dofs_per_cell; ++j)
	       {
		 comp_j = vector_fe.system_to_component_index(j).first;
		 if (comp_i == comp_j) 
		   {
		     local_normals_matrix(i,j) += vector_fe_v.shape_value(i,q)*
						  vector_fe_v.shape_value(j,q)*
						  vector_fe_v.JxW(q);
		   }
	       }
	   local_normals_rhs(i) += (vector_fe_v.shape_value(i, q)) *
                                    vector_node_normals[q](comp_i) * vector_fe_v.JxW(q);
	   }
        
       vector_cell->get_dof_indices (vector_local_dof_indices);
       
       vector_constraints.distribute_local_to_global
       (local_normals_matrix,
	local_normals_rhs,
	vector_local_dof_indices,
	vector_normals_matrix,
	vector_normals_rhs);
     }

   SparseDirectUMFPACK normals_inverse;
   normals_inverse.initialize(vector_normals_matrix);
   normals_inverse.vmult(vector_normals_solution, vector_normals_rhs);
   vector_constraints.distribute(vector_normals_solution);

   node_normals.clear();
   node_normals.resize(dh.n_dofs());   
   for (unsigned int i=0; i<vector_dh.n_dofs()/3; ++i)
       node_normals[i] = Point<3>(vector_normals_solution(3*i),
                                  vector_normals_solution(3*i+1),
                                  vector_normals_solution(3*i+2));


}




void NumericalTowingTank::compute_constraints(ConstraintMatrix &cc) {

  // we start clearing the constraint matrices
  cc.clear();

  // here we prepare the constraint matrices so as to account for the presence hanging
  // nodes
  DoFTools::make_hanging_node_constraints (vector_dh,cc);

  cc.close();
}


                                      // in the first layer of water cells past
                                      // the transom there can't be hanging nodes:
                                      // this method removes them
void NumericalTowingTank::remove_transom_hanging_nodes()
{
cout<<"Removing hanging nodes from transom stern..."<<endl;

cout<<"dofs before: "<<dh.n_dofs()<<endl;

    unsigned int refinedCellCounter = 0;
    unsigned int cycles_counter = 0;
    while(refinedCellCounter)
     {
     refinedCellCounter = 0;
     for (unsigned int i=0;i<dh.n_dofs();++i)
         {
         if ((flags[i] & transom_on_water) )
         //if ((flags[i] & water) && (flags[i] & near_boat))
            {
            //cout<<i<<": "<<support_points[i]<<endl;
            std::vector<cell_it>  cells = dof_to_elems[i];
            for (unsigned int k=0; k<cells.size(); ++k)
                {
                //cout<<k<<":  "<<cells[k]<<"   ("<<cells[k]->center()<<")"<<endl;
                for (unsigned int j=0; j<GeometryInfo<2>::faces_per_cell; ++j)
                    {
                    //cout<<"j: "<<j<<"  nb: "<<cells[k]->neighbor_index(j)<<"  ("<<endl;
                    if (cells[k]->neighbor_index(j) != -1)
                       if (cells[k]->neighbor(j)->at_boundary() && cells[k]->neighbor_is_coarser(j))
                          {
                          //cout<<"FOUND: "<<cells[k]->neighbor(j)<<" ("<<cells[k]->neighbor(j)->center()<<")"<<endl;
                          cells[k]->neighbor(j)->set_refine_flag();
                          refinedCellCounter++;
                          //cout<<"Cycle: "<<cycles_counter<<"  Refined Cells: "<<refinedCellCounter<<endl;
                          }
                    }
                }
            }
         cycles_counter++;
         if (cycles_counter > 20)
            {
            cout<<"Warning! Maximum number of transom stern edge uniforming cycles reached!"<<endl;
            break;
            }
         }

     //cout<<"refinedCellCounter   "<<refinedCellCounter<<endl;
     tria.execute_coarsening_and_refinement();
     dh.distribute_dofs(fe);
     vector_dh.distribute_dofs(vector_fe);  
       
     map_points.reinit(vector_dh.n_dofs());
     smoothing_map_points.reinit(vector_dh.n_dofs());
     old_map_points.reinit(vector_dh.n_dofs());
     ref_points.resize(vector_dh.n_dofs());
     DoFTools::map_dofs_to_support_points<2,3>(StaticMappingQ1<2,3>::mapping,
					       vector_dh, ref_points);
     generate_double_nodes_set();
     make_edges_conformal(true);
     make_edges_conformal(true);
     full_mesh_treatment();
     cycles_counter++;
     } 
//*/
cout<<"dofs after: "<<dh.n_dofs()<<endl;
cout<<"...Done removing hanging nodes from transom stern"<<endl;
}



                                      // this routine detects if mesh is not
                                      // conformal at edges (because of double
                                      // nodes) and makes the refinements needed
                                      // to make it conformal
void NumericalTowingTank::make_edges_conformal(bool isotropic_ref_on_opposite_side)
{
cout<<"Restoring mesh conformity..."<<endl;

cout<<"dofs before: "<<dh.n_dofs()<<endl;

double tol=1e-7;
     for (unsigned int i=0;i<dh.n_dofs();++i)
         {
         if ((flags[i] & edge) &&
             (double_nodes_set[i].size() == 1) )
            { 
            //we identify here the two vertices of the parent cell on the considered side
            //(which is the one with the non conformal node)
            vector<Point<3> > nodes(2);
            for (unsigned int kk=0; kk<2;++kk)
                {
                DoFHandler<2,3>::cell_iterator cell = dof_to_elems[i][kk];
                for (unsigned int f=0; f<GeometryInfo<2>::faces_per_cell; ++f)
                    {
	            if (cell->face(f)->at_boundary())
                       {
                       if (ref_points[3*i].distance(cell->face(f)->vertex(0)) <tol)
                          nodes[kk] = cell->face(f)->vertex(1);
                       else if (ref_points[3*i].distance(cell->face(f)->vertex(1)) <tol)
                          nodes[kk] = cell->face(f)->vertex(0);
                       }
                    }
                }
            // we can now compute the center of the parent cell face
            Point<3> parent_face_center = 0.5*(nodes[0]+nodes[1]);
            // now we look for the opposite side cell that has a face on an edge, having the same center
            DoFHandler<2,3>::cell_iterator cell1 = dof_to_elems[i][0]->parent();
            for(unsigned int f=0; f<GeometryInfo<2>::faces_per_cell; ++f)
	       if (cell1->face(f)->at_boundary())
                  {
	          for (std::set<tria_it>::iterator jt=edge_cells.begin(); jt != edge_cells.end(); ++jt)       
	               for (unsigned int d=0; d<GeometryInfo<2>::faces_per_cell; ++d)
		           if ((*jt)->face(d)->at_boundary())
		              if ( parent_face_center.distance(((*jt)->face(d)->vertex(0)+(*jt)->face(d)->vertex(1))/2) < tol)
                                 {
                                 // if we are on wall or free surf, use isotropic refinement
                                 if ( isotropic_ref_on_opposite_side )//(*jt)->material_id() == free_sur_ID1 ||
                                      //(*jt)->material_id() == free_sur_ID2 ||
                                      //(*jt)->material_id() == free_sur_ID3 ||
                                      //(*jt)->material_id() == wall_sur_ID1 ||
                                      //(*jt)->material_id() == wall_sur_ID2 ||
                                      //(*jt)->material_id() == wall_sur_ID3 )
                                    (*jt)->set_refine_flag();
                                 // otherwise, use anisotropic refinement to make edge mesh conformal
                                 else
                                    {
                                    if ((d==0) || (d==1))
                                       (*jt)->set_refine_flag(RefinementCase<2>::cut_axis(1));
                                    else
                                       (*jt)->set_refine_flag(RefinementCase<2>::cut_axis(0));
                                    }
                                 }
                 }
            }           
         }

  	//std::cout << "Refined counter: " << refinedCellCounter << std::endl;
        tria.execute_coarsening_and_refinement();
        dh.distribute_dofs(fe);
        vector_dh.distribute_dofs(vector_fe);  
        map_points.reinit(vector_dh.n_dofs());
        smoothing_map_points.reinit(vector_dh.n_dofs());
        old_map_points.reinit(vector_dh.n_dofs());
        ref_points.resize(vector_dh.n_dofs());
        DoFTools::map_dofs_to_support_points<2,3>(StaticMappingQ1<2,3>::mapping,
					          vector_dh, ref_points);
        generate_double_nodes_set();
        full_mesh_treatment();

cout<<"dofs after: "<<dh.n_dofs()<<endl;
cout<<"...Done restoring mesh conformity"<<endl;
}
                                      // this routine detects if mesh elements have
                                      // high aspect ratio and performs anisotropic
                                      // refinements until all aspect ratios are below 1.5

void NumericalTowingTank::remove_mesh_anisotropy(Triangulation<2,3> &tria)
  {
    Triangulation<2,3>::active_cell_iterator
      cell = tria.begin_active(), endc = tria.end();
    unsigned int refinedCellCounter = 1;
    unsigned int cycles_counter = 0;
    while(refinedCellCounter)
     {
	refinedCellCounter = 0;
	for (cell=tria.begin_active(); cell!= endc;++cell)
	  {
            //if (  cell->material_id() == wall_sur_ID1 ||
            //      cell->material_id() == wall_sur_ID2 ||
            //      cell->material_id() == wall_sur_ID3 )//( cell->center().distance(Point<3>(0.0,0.0,0.0)) <
            //     fmax(6.0*boat_model.boatWetLength/pow(2.0,double(cycles_counter)),boat_model.boatWetLength) )
            //   {
	       if (cell->extent_in_direction(0) > max_aspect_ratio*cell->extent_in_direction(1))
	          {
		  cell->set_refine_flag(RefinementCase<2>::cut_axis(0));
		  refinedCellCounter++;
	          }
	       else
	          {
		  if (cell->extent_in_direction(1) >max_aspect_ratio*cell->extent_in_direction(0))
		     {
		     cell->set_refine_flag(RefinementCase<2>::cut_axis(1));
		     refinedCellCounter++;
		     }
	          }
             //  }
	  }
     //cout<<refinedCellCounter<<endl;
     tria.execute_coarsening_and_refinement();
     dh.distribute_dofs(fe);
     vector_dh.distribute_dofs(vector_fe);  
       
     map_points.reinit(vector_dh.n_dofs());
     smoothing_map_points.reinit(vector_dh.n_dofs());
     old_map_points.reinit(vector_dh.n_dofs());
     ref_points.resize(vector_dh.n_dofs());
     DoFTools::map_dofs_to_support_points<2,3>(StaticMappingQ1<2,3>::mapping,
					       vector_dh, ref_points);
     generate_double_nodes_set();


     make_edges_conformal(false);
     make_edges_conformal(false);
     make_edges_conformal(false);
     full_mesh_treatment();
     cycles_counter++;

     //std::string filename = ( "meshResult_" +
	//		   Utilities::int_to_string(int(round(cycles_counter))) +
	//		   ".inp" );
     //std::ofstream logfile(filename.c_str());
     //GridOut grid_out;
     //grid_out.write_ucd(tria, logfile);

     }

/*
    cell = tria.begin_active();
    refinedCellCounter = 1;
    cycles_counter = 0;
    while(refinedCellCounter)
     {
	refinedCellCounter = 0;
	for (cell=tria.begin_active(); cell!= endc;++cell)
	  {
            if (  cell->material_id() == free_sur_ID1 ||
                  cell->material_id() == free_sur_ID2 ||
                  cell->material_id() == free_sur_ID3 )//( cell->center().distance(Point<3>(0.0,0.0,0.0)) <
            //     fmax(6.0*boat_model.boatWetLength/pow(2.0,double(cycles_counter)),boat_model.boatWetLength) )
               {
	       if (cell->extent_in_direction(0) > max_aspect_ratio*cell->extent_in_direction(1))
	          {
		  cell->set_refine_flag(RefinementCase<2>::cut_axis(0));
		  refinedCellCounter++;
	          }
	       else
	          {
		  if (cell->extent_in_direction(1) >max_aspect_ratio*cell->extent_in_direction(0))
		     {
		     cell->set_refine_flag(RefinementCase<2>::cut_axis(1));
		     refinedCellCounter++;
		     }
	          }
               }
	  }
     //cout<<refinedCellCounter<<endl;
     tria.execute_coarsening_and_refinement();
     dh.distribute_dofs(fe);
     vector_dh.distribute_dofs(vector_fe);  
       
     map_points.reinit(vector_dh.n_dofs());
     smoothing_map_points.reinit(vector_dh.n_dofs());
     old_map_points.reinit(vector_dh.n_dofs());
     ref_points.resize(vector_dh.n_dofs());
     DoFTools::map_dofs_to_support_points<2,3>(StaticMappingQ1<2,3>::mapping,
					       vector_dh, ref_points);
     generate_double_nodes_set();
     make_edges_conformal(false);
     make_edges_conformal(false);
     make_edges_conformal(false);
     full_mesh_treatment();
     cycles_counter++;
     }

*/

/*    Triangulation<2,3>::active_cell_iterator
      cell = tria.begin_active(), endc = tria.end();
    unsigned int refinedCellCounter = 1;
    double counter=0;
    while(refinedCellCounter)
      {
	refinedCellCounter = 0;
	for (cell=tria.begin_active(); cell!= endc;++cell)
	  {
	    if ( (cell->extent_in_direction(0) > 1.5*cell->extent_in_direction(1))    &&
                 (fabs(cell->center()(1)) < (Ly_domain/4+Lx_boat/2)/pow(2.0,counter)) && 
                 (cell->material_id() == free_sur_ID1 ||
                  cell->material_id() == free_sur_ID2 ||
                  cell->material_id() == free_sur_ID3 ) )
	      {
		cell->set_refine_flag(RefinementCase<2>::cut_axis(0));
		refinedCellCounter++;
	      }
	    else
	      {
		if ( (cell->extent_in_direction(1) > 1.5*cell->extent_in_direction(0))    &&
                     (fabs(cell->center()(1)) < (Ly_domain/4+Lx_boat/2)/pow(2.0,counter)) &&
                     (cell->material_id() == free_sur_ID1 ||
                      cell->material_id() == free_sur_ID2 ||
                      cell->material_id() == free_sur_ID3 ) ) 
		  {
		    cell->set_refine_flag(RefinementCase<2>::cut_axis(1));
		    refinedCellCounter++;
		  }
	      }
	  }
     counter = counter+1.0;
     tria.execute_coarsening_and_refinement();
     dh.distribute_dofs(fe);
     vector_dh.distribute_dofs(vector_fe);  
     map_points.reinit(vector_dh.n_dofs());
     smoothing_map_points.reinit(vector_dh.n_dofs());
     old_map_points.reinit(vector_dh.n_dofs());
     ref_points.resize(vector_dh.n_dofs());
     DoFTools::map_dofs_to_support_points<2,3>(StaticMappingQ1<2,3>::mapping,
					       vector_dh, ref_points);
     generate_double_nodes_set();
     make_edges_conformal(false);
     make_edges_conformal(false);
     full_mesh_treatment();
     cout<<"Current dofs number: "<<dh.n_dofs()<<endl;
     std::string filename = ( "meshResult_" +
			   Utilities::int_to_string(int(round(counter))) +
			   ".vtu" );
     std::ofstream logfile(filename.c_str());
     GridOut grid_out;
     grid_out.write_ucd(tria, logfile);
     }


    cell = tria.begin_active(), endc = tria.end();
    refinedCellCounter = 1;
    counter=0;
    while(refinedCellCounter)
      {
       std::string filename0 = ( "meshResultZeroth_" +
     Utilities::int_to_string(int(round(counter))) +
			   ".vtu" );
     std::ofstream logfile0(filename0.c_str());
     GridOut grid_out0;
     grid_out0.write_ucd(tria, logfile0);
	refinedCellCounter = 0;
	for (cell=tria.begin_active(); cell!= endc;++cell)
	  {
	    if ( (cell->extent_in_direction(0) > 2.0*cell->extent_in_direction(1))    &&
                 (cell->material_id() == wall_sur_ID1 ||
                  cell->material_id() == wall_sur_ID2 ||
                  cell->material_id() == wall_sur_ID3 ) )
	      {
		cell->set_refine_flag(RefinementCase<2>::cut_axis(0));
		refinedCellCounter++;
	      }
	    else
	      {
		if ( (cell->extent_in_direction(1) > 2.0*cell->extent_in_direction(0))    &&
                     (cell->material_id() == wall_sur_ID1 ||
                      cell->material_id() == wall_sur_ID2 ||
                      cell->material_id() == wall_sur_ID3 ) ) 
		  {
		    cell->set_refine_flag(RefinementCase<2>::cut_axis(1));
		    refinedCellCounter++;
		  }
	      }
	  }
     counter = counter+1.0;
     tria.execute_coarsening_and_refinement();
     dh.distribute_dofs(fe);
     vector_dh.distribute_dofs(vector_fe);  
     map_points.reinit(vector_dh.n_dofs());
     smoothing_map_points.reinit(vector_dh.n_dofs());
     old_map_points.reinit(vector_dh.n_dofs());
     ref_points.resize(vector_dh.n_dofs());
     DoFTools::map_dofs_to_support_points<2,3>(StaticMappingQ1<2,3>::mapping,
					       vector_dh, ref_points);
     generate_double_nodes_set();
     std::string filename1 = ( "meshResultFirst_" +
     Utilities::int_to_string(int(round(counter))) +
			   ".vtu" );
     std::ofstream logfile1(filename1.c_str());
     GridOut grid_out1;
     grid_out1.write_ucd(tria, logfile1);
     make_edges_conformal(false);
     make_edges_conformal(false);
     full_mesh_treatment();
     cout<<"Current dofs number: "<<dh.n_dofs()<<endl;
     std::string filename = ( "meshResultSecond_" +
			   Utilities::int_to_string(int(round(counter))) +
			   ".vtu" );
     std::ofstream logfile(filename.c_str());
     GridOut grid_out;
     grid_out.write_ucd(tria, logfile);
     }
//*/
     cout<<"Current dofs number: "<<dh.n_dofs()<<endl;
     std::string filename = ( "meshResultSecond.vtu" );
     std::ofstream logfile(filename.c_str());
     GridOut grid_out;
     grid_out.write_ucd(tria, logfile);

  }

void NumericalTowingTank::refine_global_on_boat(const unsigned int num_refinements)
{
for (unsigned int i=0; i<num_refinements;++i)
    {
    std::cout<<"Uniform boat refinement cycle "<<i+1<<" of "<<num_refinements<<std::endl;
    tria_it cell = tria.begin_active(), endc = tria.end();
    for (cell=tria.begin_active(); cell!= endc;++cell)
        {
        if ((cell->material_id() == wall_sur_ID1 ||
             cell->material_id() == wall_sur_ID2 ||
             cell->material_id() == wall_sur_ID3 ))
            cell->set_refine_flag();
        }
    tria.execute_coarsening_and_refinement();
    dh.distribute_dofs(fe);
    vector_dh.distribute_dofs(vector_fe);  
    map_points.reinit(vector_dh.n_dofs());
    smoothing_map_points.reinit(vector_dh.n_dofs());
    old_map_points.reinit(vector_dh.n_dofs());
    ref_points.resize(vector_dh.n_dofs());
    DoFTools::map_dofs_to_support_points<2,3>(StaticMappingQ1<2,3>::mapping,
				  	      vector_dh, ref_points);
    generate_double_nodes_set();
    full_mesh_treatment();
    make_edges_conformal(true);
    make_edges_conformal(true);
    }
/*
for (unsigned int i=0; i<4-num_refinements;++i)
    {
    tria_it cell = tria.begin_active(), endc = tria.end();
    for (cell=tria.begin_active(); cell!= endc;++cell)
        {
        if ( (cell->material_id() == wall_sur_ID1 ||
              cell->material_id() == wall_sur_ID2 ||
              cell->material_id() == wall_sur_ID3 ) &&
             ((cell->center()(0) < -1.95) && (cell->center()(2) < -0.15)) )
            cell->set_refine_flag();
        }
    tria.execute_coarsening_and_refinement();
    dh.distribute_dofs(fe);
    vector_dh.distribute_dofs(vector_fe);  
    map_points.reinit(vector_dh.n_dofs());
    smoothing_map_points.reinit(vector_dh.n_dofs());
    old_map_points.reinit(vector_dh.n_dofs());
    ref_points.resize(vector_dh.n_dofs());
    DoFTools::map_dofs_to_support_points<2,3>(StaticMappingQ1<2,3>::mapping,
					      vector_dh, ref_points);
    generate_double_nodes_set();
    full_mesh_treatment();
    make_edges_conformal(true);
    make_edges_conformal(true);
    }
//*/
/*
for (unsigned int i=0; i<4;++i)
    {
    tria_it cell = tria.begin_active(), endc = tria.end();
    for (cell=tria.begin_active(); cell!= endc;++cell)
        {
        if ( (cell->material_id() == free_sur_ID1 ||
              cell->material_id() == free_sur_ID2 ||
              cell->material_id() == free_sur_ID3 ) &&
             (cell->center()(0) > -4.10)            &&
             (cell->center()(0) < 6.10)             &&
             (cell->center()(1) > -3.5)            &&
             (cell->center()(1) < 3.5)             &&
             (cell->diameter()/2.8 > 0.174353)          )
           cell->set_refine_flag();
        }
    tria.execute_coarsening_and_refinement();
    dh.distribute_dofs(fe);
    vector_dh.distribute_dofs(vector_fe);  
    map_points.reinit(vector_dh.n_dofs());
    smoothing_map_points.reinit(vector_dh.n_dofs());
    old_map_points.reinit(vector_dh.n_dofs());
    ref_points.resize(vector_dh.n_dofs());
    DoFTools::map_dofs_to_support_points<2,3>(StaticMappingQ1<2,3>::mapping,
					      vector_dh, ref_points);
    generate_double_nodes_set();
    full_mesh_treatment();
    make_edges_conformal(true);
    make_edges_conformal(true);
    }
//*/

/*
for (unsigned int k=0; k<1; ++k)
    {
    Triangulation<2,3>::active_cell_iterator
    cell = tria.begin_active(), endc = tria.end();
    for (cell=tria.begin_active(); cell!= endc;++cell)
        {
        if ( (cell->material_id() == wall_sur_ID1 ||
             cell->material_id() == wall_sur_ID2 ||
             cell->material_id() == wall_sur_ID3 ) &&
             (cell->center()(0) > -3.50) &&
             (cell->center()(0) < 17.00) &&
             (fabs(cell->center()(1)) < 0.363970234*(cell->center()(0)+3.50)) &&
             (cell->center()(2) > -0.1) &&
             (cell->diameter() > 0.4) )
             cell->set_refine_flag();
        }
    tria.execute_coarsening_and_refinement();
    dh.distribute_dofs(fe);
    vector_dh.distribute_dofs(vector_fe);  
    map_points.reinit(vector_dh.n_dofs());
    smoothing_map_points.reinit(vector_dh.n_dofs());
    old_map_points.reinit(vector_dh.n_dofs());
    ref_points.resize(vector_dh.n_dofs());
    DoFTools::map_dofs_to_support_points<2,3>(StaticMappingQ1<2,3>::mapping,
					      vector_dh, ref_points);
    generate_double_nodes_set();
    full_mesh_treatment();
    make_edges_conformal(true);
    make_edges_conformal(true);
    }
//*/
}


void NumericalTowingTank::extract_boundary_dofs(std::vector<bool> &dofs, unsigned int id,
			                        DoFHandler<2,3> &vector_dh) 
{
  std::vector<bool> comp_sel(3, true);
  std::set<unsigned char> ids;
  ids.insert(id);
  DoFTools::extract_boundary_dofs(vector_dh, comp_sel, dofs, ids);
}

  
unsigned int NumericalTowingTank::find_point_id(const Point<3> &p, const vector<Point<3> > &ps) 
{
  //cout<<"numPoints "<<ps.size()<<endl;
  for (unsigned int i=0; i<ps.size()/3; ++i)
      {
      //cout<<"i "<<i<<" p("<<ps[3*i]<<")  d "<< p.distance(ps[3*i])<<endl;
      if(p.distance(ps[3*i]) < 1e-7)
        return(i);
      }
  return 0;
}


Handle(Geom_Curve) NumericalTowingTank::get_curve(const vector<Point<3> > &ps,
			                          const vector<bool> &id,
                                                  const Point<3> direction) 
{
  vector<Point<3> > points;
  for(unsigned int i=0; i<ps.size()/3; ++i)
    if(id[3*i] == true)
      points.push_back(ps[3*i]);

  return interpolation_curve_points_sort(points, direction);
}

class NumericalTowingTank;
