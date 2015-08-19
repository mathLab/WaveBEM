
//----------------------------  template.cc  ---------------------------
//    $Id: testsuite.html 13373 2006-07-13 13:12:08Z kanschat $
//    Version: $Name$ 
//
//    Copyright (C) 2005 by the deal.II authors 
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  template.cc  ---------------------------


// a short (a few lines) description of what the program does

//#include "../tests.h"
// Read goteborg.iges and dump its topological structure to the logfile.

#ifndef numerical_towing_tank_h
#define numerical_towing_tank_h

#include <fstream>
#include <iostream>
#include <cmath>
#include <limits>
#include <vector>

#include <deal.II/base/logstream.h>

#include "occ_utilities.h"
#include "occ_normal_projection.h"
#include "occ_arclength_projection.h"
#include "occ_axis_projection.h"

#include <TopTools.hxx>
#include <Standard_Stream.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Edge.hxx>


#include <IGESControl_Reader.hxx>
#include <IGESControl_Controller.hxx>
#include <IGESControl_Writer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Face.hxx>
#include <BRepTools.hxx>
#include <XSControl_Reader.hxx>
#include <TopTools_SequenceOfShape.hxx>
#include <Handle_Standard_Transient.hxx>
#include <TColStd_SequenceOfTransient.hxx>
#include <TColStd_HSequenceOfTransient.hxx>
#include <TopExp_Explorer.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <Standard_Real.hxx>
#include <Standard_Integer.hxx>
#include <BRep_Tool.hxx>
#include <Geom_Surface.hxx>
#include <Geom_Plane.hxx>
#include <Prs3d_ShapeTool.hxx>
#include <Bnd_Box.hxx>
#include <gp_Ax3.hxx>
#include <gp_Pln.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <GeomAPI_IntSS.hxx>
#include <TopoDS_Wire.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <Geom_Curve.hxx>
#include <Geom_BoundedCurve.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <TColGeom_Array1OfCurve.hxx>
#include <GeomAPI_IntCS.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <GeomLib_Tool.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <GeomConvert_CompCurveToBSplineCurve.hxx>

#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/base/point.h>
#include "occ_line_smoothing.h"
#include "surface_smoothing.h"
#include "computational_domain.h"


class NumericalTowingTank : public ComputationalDomain<3>
{
  public:
    
                                      // constructor: since this is the
				      // class containing all the geometry and
				      // the base instruments needed by all the
				      // other classes, it is created first and
				      // the constructor does not need
				      // arguments.
				      // For the same reason, most of the class
				      // attributes are public: we can leter
				      // make them public end introduce suitable
				      // Get and Set methods, if needed
    
    NumericalTowingTank(unsigned int, unsigned int);

    void full_mesh_treatment();
                                      // method to compute interpolated curvatures
				      // for geometrically conformal mesh refinement
    void compute_curvatures(Vector<double> &curvatures);
                                      // method to use interpolated curvatures
				      // for geometrically conformal mesh refinement
    void apply_curvatures(const Vector<double> &curvatures,  const std::vector<bool> boundary_dofs);

    void partial_mesh_treatment(const double blend_factor);

    void update_mapping(const Vector<double> &map_points);

    ~NumericalTowingTank();

    virtual void read_domain();
    virtual void refine_and_resize();
    virtual void generate_double_nodes_set();

                                      // method to declare the parameters
				      // to be read from the parameters file
    
    virtual void declare_parameters(ParameterHandler &prm);
    
                                      // method to parse the needed parameters
				      // from the parameters file
  
    virtual void parse_parameters(ParameterHandler &prm);
    
                                      // method to parse mesh from the
				      // input file
    
    void create_initial_mesh(const Point<3> PointFrontTop,
                             const Point<3> PointFrontBot,
                             const Point<3> PointMidTop,
                             const Point<3> PointMidBot,
                             const Point<3> PointBackTop,
                             const Point<3> PointBackBot,
                             const Point<3> PointLeftTransom,
                             const Point<3> PointRightTransom,
                             const Point<3> PointCenterTransom,
                             Triangulation<2,3>  &triangulation);
    
    void refine_global_on_boat(const unsigned int num_refinements);

    void compute_nodes_flags();

    void set_up_smoother();

    void initialize_smoother();

    void update_smoother();

    void perform_line_smoothing(unsigned int num_smoothings);

    void perform_surface_projection();

    void perform_water_line_nodes_projection();

    void perform_smoothing(bool full_treatment, const double blend_factor);

    void compute_normals_at_nodes(Vector<double> &map_points_used);

    void compute_constraints(ConstraintMatrix &cc);

                                      // this routine detects if mesh is not
                                      // conformal at water/boat edges (because of double
                                      // nodes) and makes the refinements needed
                                      // to make it conformal
    void make_edges_conformal(bool isotropic_ref_on_opposite_side);
                                      // this routine detects if mesh elements have
                                      // high aspect ratio and performs anisotropic
                                      // refinements until all aspect ratios are below 1.5
    void remove_mesh_anisotropy(Triangulation<2,3> &tria);

                                      // in the first layer of water cells past
                                      // the transom there can't be hanging nodes:
                                      // this method removes them
    void remove_transom_hanging_nodes();

                                      // this routine finds the index of the vector point p
    unsigned int find_point_id(const Point<3> &p, const std::vector<Point<3> > &ps);
                                      // method to get the ids of all dofs lying on a boundary
    void extract_boundary_dofs(std::vector<bool> &dofs, unsigned int id,
			       DoFHandler<2,3> &vector_dh);
                                      // method to obtain an approximated geom_curve  
    Handle(Geom_Curve) get_curve(const std::vector<Point<3> > &ps,
			         const std::vector<bool> &id,
                                 const Point<3> direction);

                                      // function that computes (vectorial)
                                      // distances of each node from the reference
                                      // (CAD) surface. Retunrn zero if point is
                                      // on the surface or if no reference surface
                                      // is specified for the node considered
    void evaluate_ref_surf_distances(Vector <double> &distances, const bool only_surf_smoothing=false);


 				   // number of refining cycles to be performed on the boat 				      
unsigned int n_cycles;
 				   // name of the iges file to be used 				      
std::string iges_file_name;
				   // we first need a boat model, which
                                   // will provide deatails about the
                                   // boat geometry
BoatModel boat_model;
				   // surface smoother for mesh refinements: will use mapping
SurfaceSmoothing *restart_surface_smoother;
				   // surface smoother for mesh smoothing: will use static mapping
SurfaceSmoothing *surface_smoother;
                                   // vector containing all the line smoothers in the
                                   // following order
				   // 0. line smoothing class for front part of the keel
				   // 1. line smoothing class for left transom edge/rear part of the keel
                                   // 2. line smoothing class for right transom edge/rear part of the keel
				   // 3. line smoothing class for right front water line
				   // 4. line smoothing class for left front water line
				   // 5. line smoothing class for right rear water line
				   // 6. line smoothing class for left rear water line
std::vector<OpenCascade::LineSmoothing *> line_smoothers;

				   // vector of vectors of booleans needed to decide which
                                   // nodes can be displaced by each smoothing
                                   // and which ones are fixed. order is the same as line_smoothers
std::vector< std::vector<bool> > boundary_dofs;
                                   // contains the ids of the boundaries involved in each
                                   // line smoothing 
std::vector<unsigned int> boundary_ids;
                                   // contains the base points of each line smoothing 
std::vector<Point<3> > base_points;
                                   // contains the moving points of each line smoothing 
std::vector<Point<3> > moving_points;
                                   // contains the base points ids of each line smoothing 
std::vector<unsigned int> base_point_ids;
                                   // contains the moving points ids of each line smoothing 
std::vector<unsigned int> moving_point_ids;
                                   // contains the curves of each line smoothing 
std::vector<Handle(Geom_Curve)> curves;
                                   // contains booleans if curve must be rebuilt at each
                                   // call of smoothing routine. false for keel smoothings
                                   // true for water line smoothings 
std::vector<bool> on_curve_option;
                                   // locations are needed when reference configuration
                                   // is roto-translated
std::vector< TopLoc_Location *> smoothers_locations;
 				   // nodes flags on the scalar dof_handler
std::vector<GeometryFlags> flags;
 				   // nodes flags on the vector dof_handler
std::vector<GeometryFlags> vector_flags;
 				   // vector containing the normals on boat nodes
                                   // zero vectors on other nodes
std::vector< Point<3> > iges_normals;
 				   // vector containing the normals on boat nodes
                                   // zero vectors on other nodes, but referred to
                                   // grid after remesh
std::vector< Point<3> > old_iges_normals;
 				   // vector containing the mean curvatures on boat
                                   // nodes, and zeros on other nodes
std::vector<double> iges_mean_curvatures;
                                   // sparsity pattern for the normals problem
SparsityPattern      normals_sparsity_pattern;
                                   // constraint matrix for normals problem
ConstraintMatrix     vector_constraints;
                                   // matrix for the problem for node normals computation
SparseMatrix<double> vector_normals_matrix;
                                   // vector for the rhs of the problem for node normals computation
Vector<double>       vector_normals_rhs;
				   // solution vector to problem for node normals computation
Vector<double>       vector_normals_solution;
				   // node normals at mesh nodes
std::vector<Point<3> > node_normals;
                                   // set containing cells on edges
std::set<tria_it> edge_cells;
                                   // set containing cells on boat bordering with free surface
std::set<tria_it> boat_edge_cells;
                                   // set containing cells on water bordering with boat
std::set<tria_it> water_edge_cells;
                                   // map containing cells on boat and bordering cell on free surface
std::map<tria_it, tria_it> boat_to_water_edge_cells;
                                   // map containing cells on boat and bordering cell on free surface
std::map<tria_it, tria_it> water_to_boat_edge_cells;
                                   // vector containing right hand side for smoothing beltrami problem
Vector<double> smoothing_curvature_vector;
                                   // vector containing the euler vector for smoothing (always initialized
                                   // as a copy of map_points)
Vector<double> smoothing_map_points;
                                   // the euler vector at the simulation start
Vector<double> initial_map_points;
                                   // the euler vector obtained right after each restart
Vector<double> old_map_points;
                                   // x domain dimension
double Lx_domain;
                                   // y domain dimension
double Ly_domain;
                                   // z domain dimension
double Lz_domain;
                                   // boat wet length
double Lx_boat;

const unsigned int mapping_degree;

                                   // vector with edges tangents at dofs
Vector<double > edges_tangents;
                                   // vector with edges dofs length ratios
Vector<double > edges_length_ratios;
                                   // boat displacement (in Kg) to be passed to boat model class
double boat_displacement;
                                   // this is the sink (positive downwards) to be applied to the hydrostatic
                                   // equilibrium position
double assigned_sink;
                                   // this is the trim (positive when bow goes up) computed from the
                                   // hydrostatic equilibrium position (which is assumed to be the angular
                                   // position assigned in the CAD file)
double assigned_trim;
                                   // here we have the maximum aspect ratio allowed for the cells
double max_aspect_ratio;
                                   // the next is a factor determining the inclination of the mesh longitudinal
                                   // lines on the front
double front_mesh_inclination_coeff;
                                   // the next is a factor determining the inclination of the mesh longitudinal
                                   // lines on the back
double back_mesh_inclination_coeff;
                                   // this determines the location of the stern bottom point of the first quadrilaterals
                                   // built on the hull (in % of the total immersed length of the keel, starting form
                                   // the aft perpendicular)
double back_keel_length;
                                   // this determines the location of the bow bottom point of the first quadrilaterals
                                   // built on the hull (in % of the total immersed length of the keel, starting form
                                   // the fore perpendicular)
double front_keel_length;
                                   // this determines the location of the central bottom point of the first quadrilaterals
                                   // built on the hull (in % of the total immersed length of the keel, starting form
                                   // the fore perpendicular)
double middle_keel_length;
                                   // number of edges composing the transom stern (temporary)
unsigned int number_of_transom_edges;
                                   // number of uniform refinements on the boat surface requested for the initial mesh
unsigned int init_global_boat_refs;
                                   // number of non uniform (curvature based) refinements on the boat surface requested
                                   // for the initial mesh
unsigned int init_adaptive_boat_refs;
                                   // fraction of boat cells to be refined per each cycle in the initial curvature
                                   // based refinement
double init_adaptive_boat_refs_fraction;
                                   // a flag that determines if a boat surface has to be used or the numerical wave tank
                                   // has to be prepared for the simulation of free waves
bool no_boat;

};

#endif
