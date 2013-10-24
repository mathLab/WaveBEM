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
				 

#include "../include/boat_model.h"
#include <gp_Trsf.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepBuilderAPI_GTransform.hxx>

#include <TopTools.hxx>
#include <Standard_Stream.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Edge.hxx>
#include <GC_MakeSegment.hxx>

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
#include <GeomLProp_SLProps.hxx>
#include <BRepTools_ReShape.hxx>

#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/grid_tools.h>

#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/grid_out.h>
#include <grid/tria_iterator.h>
#include <grid/grid_generator.h>
#include <grid/tria_boundary_lib.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <dofs/dof_renumbering.h>
#include <fe/fe_q.h>
#include <fe/fe_values.h>
#include <fe/fe_system.h>
#include <fe/mapping_q_eulerian.h>
#include <fe/mapping_q.h>
#include <numerics/vector_tools.h>


using namespace dealii;
using namespace OpenCascade;

BoatModel::BoatModel() :
		boat_surface_right(NULL),
		boat_surface_left(NULL),
		boat_water_line_right(NULL),
		boat_water_line_left(NULL),
		boat_keel(NULL),
		boat_keel_norm(NULL),
                boat_transom_left(NULL),
                boat_transom_right(NULL),
                wake_line_left(NULL),
                wake_line_right(NULL),
		water_line_right(NULL),
		water_line_left(NULL)
{
}

BoatModel::~BoatModel()
{
  if(boat_surface_right)
    delete boat_surface_right;
  if(boat_surface_left)
    delete boat_surface_left;
  if(boat_water_line_right)
    delete boat_water_line_right;
  if(boat_water_line_left)
    delete boat_water_line_left;
  if(boat_keel)
    delete boat_keel;
  if(boat_keel_norm)
    delete boat_keel_norm;
  if (boat_transom_left)
    delete boat_transom_left;
  if (boat_transom_right)
    delete boat_transom_right;
  if(water_line_right)
    delete water_line_right;
  if(water_line_left)
    delete water_line_left;
}


TopoDS_Shape BoatModel::ReverseFaceOrientation(const TopoDS_Shape& shape,
                                               const TopoDS_Face& face)
{
     Handle(BRepTools_ReShape) rebuild = new BRepTools_ReShape();
     rebuild->ModeConsiderOrientation() = Standard_True;
     TopoDS_Shape newface = face.Complemented();
     rebuild->Replace(face, newface, Standard_True );
     TopoDS_Shape newshape = rebuild->Apply(shape, TopAbs_FACE );
     return newshape;
}


void BoatModel::start_iges_model(std::string igesFileName,
                                 double scale,
                                 double displacement,
                                 double assigned_sink,
                                 double assigned_trim,                                 
				 double back_keel_length,
				 double front_keel_length,
				 double middle_keel_length,
                                 unsigned int number_of_transom_edges)
{
				   //let's read the (iges) hull shape from file
  sh =read_IGES(igesFileName,scale);


    TopExp_Explorer faceExplorer(sh, TopAbs_FACE);
    TopoDS_Face face;
    gp_Pnt Pint;
    std::vector<bool> to_be_changed;
    unsigned int face_count = 0;
    while(faceExplorer.More())
      {
	face = TopoDS::Face(faceExplorer.Current());
        Standard_Real umin, umax, vmin, vmax;
        BRepTools::UVBounds(face, umin, umax, vmin, vmax);          // create surface
        Handle(Geom_Surface) surf=BRep_Tool::Surface(face);          // get surface properties
        GeomLProp_SLProps props(surf, (umin+umax)/2.0, (vmin+vmax)/2.0, 1, 0.01);          // get surface normal
        gp_Dir norm=props.Normal();   
        if (norm.Y() > 0)
           {
           to_be_changed.push_back(true);
           //cout<<"i "<<face_count<<endl;
           }
        else if (norm.Y() == 0)
           {
           if (norm.Z() < 0)
              {
              to_be_changed.push_back(true);
              }
           else
              {
              to_be_changed.push_back(false);
              }              
           }
        else
           to_be_changed.push_back(false);
	faceExplorer.Next();
        ++face_count;
      }
    //cout<<"size "<<to_be_changed.size()<<endl;
    TopoDS_Shape newShape;
    for (unsigned int i=0; i<face_count; ++i)
        {
        //cout<<"*i "<<i<<"  of "<<face_count<<endl;
        if (to_be_changed.at(i))
           {
           //cout<<"**i "<<i<<"  of "<<face_count<<endl;
           faceExplorer.Init(sh,TopAbs_FACE);
           for (unsigned int j=0; j<i; ++j)
               faceExplorer.Next();
           face = TopoDS::Face(faceExplorer.Current());
           newShape = ReverseFaceOrientation(sh,face);
           sh = newShape;
           //cout<<"***i "<<i<<"  of "<<face_count<<endl;
           }
        //cout<<"****i "<<i<<"  of "<<face_count<<endl;
        }
 //Standard::Purge();
				   //let's create the reflected side of the hull
				   //this is the y axis definition for planar mirroring
  gp_Ax2 yAx(gp_Pnt(0.0,0.0,0.0), gp_Dir(0.0,1.0,0.0), gp_Dir(1.0,0.0,0.0));
				   //here we define the y mirroring transformation
  gp_Trsf y_mirroring;
  y_mirroring.SetMirror(yAx);
				   //here we apply to the hull the y mirroring transformation
  BRepBuilderAPI_Transform boat_surface_transformation(sh, y_mirroring);
				   //we use the transformation to define the mirrored shape
  refl_sh = boat_surface_transformation.Shape();


                                   // if displacement is set to 0 we keep the boat position
                                   // in the iges file
  if (displacement > 0) 
     {
                                   //here we perform hydrostatic computations to evaluate the
                                   //correct hull sinkage
     double sink;
     double weight = displacement*9.81;
     compute_hydrostatic_sink(sink,weight);
                                   //now we have to use the computed sink to diplace vertically 
                                   //the left and right sides of the hull

     TopoDS_Shape sh_copy(sh);
     TopoDS_Shape refl_sh_copy(refl_sh);
     gp_Vec vert_displ(0.0,0.0,-sink); 
     gp_Trsf vert_translation;
     vert_translation.SetTranslation(vert_displ);
     BRepBuilderAPI_GTransform right_vert_transl(sh_copy, vert_translation);
     sh = right_vert_transl.Shape();
     BRepBuilderAPI_GTransform left_vert_transl(refl_sh_copy, vert_translation);
     refl_sh = left_vert_transl.Shape();
     }
/*
     Point<3> hs_force = compute_hydrostatic_force(0.0);
     Point<3> hs_moment = compute_hydrostatic_moment(0.0);

     double rot_axis_x_coor = -hs_moment(1)/hs_force(2);

                                  //here we prepare the rotation of the boat of the requested trim angle     
     gp_Pnt rot_center(0.0,rot_axis_x_coor,0.0);
     gp_Dir rot_dir(0.0,1.0,0.0);
     gp_Ax1 rot_axis(rot_center, rot_dir);
     gp_Trsf rotation;
     rotation.SetRotation(rot_axis,assigned_trim);
                                  //here we prepare the translation of the boat of the requested sink
     gp_Trsf translation;
     gp_Vec vrt_displ(0.0,0.0,-assigned_sink);
     translation.SetTranslation(vrt_displ);
                                  //the rotation and translation are combined in a single transformation
     gp_Trsf Tcomp = translation*rotation;
                                  //the transformation is applied to the two sides of the boat
     BRepBuilderAPI_GTransform right_transf(sh, Tcomp);
     sh = right_transf.Shape();
     BRepBuilderAPI_GTransform left_transf(refl_sh, Tcomp);
     refl_sh = left_transf.Shape();
*/   


     cout<<"The hull has been placed in the correct position"<<endl;
     Point<3> hs_force = compute_hydrostatic_force(0.0);
                                   //now the boat is in the correct position 
// These lines can be used to dump the keel edge (or other shapes) on an .igs file

IGESControl_Controller::Init();
IGESControl_Writer ICW ("MM", 0);
Standard_Boolean ok = ICW.AddShape (sh);
ICW.ComputeModel();
Standard_Boolean OK = ICW.Write ("shape.igs");
//*/



				   //here we extract the keel edge from the hull shape
  keel_edge = extract_xz_edges(sh,1e-4,30);
				   //here we extract the transom edge from the right hull shape
  right_transom_edge = extract_transom_edges(sh,number_of_transom_edges,1e-2);
                                   //here we extract the transom edge from the right hull shape
                                   //by applying
                                   //to the right water line the y mirroring transformation
  BRepBuilderAPI_Transform transom_transformation(right_transom_edge, y_mirroring);
                                   //we use the transformation to define the mirrored shape
  left_transom_edge = transom_transformation.Shape();


// These lines can be used to dump the keel edge (or other shapes) on an .igs file
/*
IGESControl_Controller::Init();
IGESControl_Writer ICW ("MM", 0);
Standard_Boolean ok = ICW.AddShape (keel_edge);
ICW.ComputeModel();
Standard_Boolean OK = ICW.Write ("keel.igs");
*/
				   //here we extract the undisturbed right water line from
                                   //the hull shape
  intersect_plane(sh,right_undist_water_line,0.0,0.0,1.0,0.0,1e-2); // 1e-2 tolerance for comacina
				   //here we extract the undisturbed left water line from
                                   //the hull shape by applying
                                   //to the right water line the y mirroring transformation
  BRepBuilderAPI_Transform waterline_transformation(right_undist_water_line, y_mirroring);
                                   //we use the transformation to define the mirrored shape
  left_undist_water_line = waterline_transformation.Shape();


  TopLoc_Location L;
  Standard_Real First;
  Standard_Real Last;
				   //we get the curve corresponding to right undisturbed
                                   // waterline
  TopExp_Explorer edge1Explorer(right_undist_water_line, TopAbs_EDGE);
  TopoDS_Edge edge1 =  TopoDS::Edge(edge1Explorer.Current());
  right_undisturbed_waterline_curve = BRep_Tool::Curve(edge1,L,First,Last);
				   //we get the curve corresponding to left undisturbed
                                   // waterline
  TopLoc_Location L2;
  Standard_Real First2;
  Standard_Real Last2;  
  TopExp_Explorer edge2Explorer(left_undist_water_line, TopAbs_EDGE);
  TopoDS_Edge edge2 =  TopoDS::Edge(edge2Explorer.Current());
  left_undisturbed_waterline_curve = BRep_Tool::Curve(edge2,L2,First2,Last2);

  gp_Vec V1;
  gp_Pnt P1;
  gp_Vec V2;
  gp_Pnt P2;  
  left_undisturbed_waterline_curve->D1(First,P1,V1);
  left_undisturbed_waterline_curve->D1(Last,P2,V2);
  double slope = 0;
  if (P1.X() > P2.X())
     slope = fabs(V1.Y()/V1.X());
  else
     slope = fabs(V2.Y()/V2.X());
  
				   //we now find the keel intersections with the xy plane
  std::vector<gp_Pnt> intPoints(2);
				   //this is the xy plane
  Handle(Geom_Plane) xyPlane = new Geom_Plane(0.,0.,1.,0.);
				   // here we intersect the keel with
				   // the xy plane

  TopExp_Explorer edgeExplorer(keel_edge, TopAbs_EDGE);
  TopoDS_Edge edge =  TopoDS::Edge(edgeExplorer.Current());
  equiv_keel_bspline = BRep_Tool::Curve(edge,L,First,Last);
  
  TopExp_Explorer edge3Explorer(left_transom_edge, TopAbs_EDGE);
  TopoDS_Edge edge3 =  TopoDS::Edge(edge3Explorer.Current());
  left_transom_bspline = BRep_Tool::Curve(edge3,L,First,Last);
  //cout<<First<<" "<<Last<<endl;
  if (left_transom_bspline->Value(First).Z() > left_transom_bspline->Value(Last).Z())
     PointCenterTransom = Pnt(left_transom_bspline->Value(Last));
  else
     PointCenterTransom = Pnt(left_transom_bspline->Value(First));

  TopExp_Explorer edge4Explorer(right_transom_edge, TopAbs_EDGE);
  TopoDS_Edge edge4 =  TopoDS::Edge(edge4Explorer.Current());
  right_transom_bspline = BRep_Tool::Curve(edge4,L,First,Last);

  GeomAPI_IntCS Intersector(equiv_keel_bspline, xyPlane);
  int npoints = Intersector.NbPoints();

  gp_Pnt gpFront(0.0,0.0,0.0);
  gp_Pnt gpBack(0.0,0.0,0.0);
  gp_Pnt transomLeft(0.0,0.0,0.0);
  gp_Pnt transomRight(0.0,0.0,0.0);

      				   //this is the xy plane
    Handle(Geom_Plane) xyPlane1 = new Geom_Plane(0.,0.,1.,0.);
    GeomAPI_IntCS IntersectorLeft(left_transom_bspline, xyPlane1);
    				   //this is the xy plane
    Handle(Geom_Plane) xyPlane2 = new Geom_Plane(0.,0.,1.,0.);    
    GeomAPI_IntCS IntersectorRight(right_transom_bspline, xyPlane2);


  cout<<"npoints "<<npoints<<endl;
  if (npoints == 2)
    {
    is_transom = false;
    if (Intersector.Point(1).X() < Intersector.Point(2).X())
       {
       gpFront = Intersector.Point(1);
       gpBack = Intersector.Point(2);
       }
    else
       {
       gpFront = Intersector.Point(2);
       gpBack = Intersector.Point(1);
       }
    boatWetLength = fabs(gpBack.X()-gpFront.X());
    Point<3> far(50*boatWetLength,0.0,0.0);
    GC_MakeSegment first_wake(Pnt(far), gpBack);
    left_wake_bspline = first_wake.Value();
    right_wake_bspline = first_wake.Value();
    }     
  else if (npoints == 1)
    {
    
    is_transom = true;
    cout<<"Number of keel intersections with xy plane: "<<npoints<<endl;

    if( (IntersectorLeft.NbPoints() != 1) || (IntersectorRight.NbPoints() != 1))
        AssertThrow((IntersectorLeft.NbPoints() == 1) && (IntersectorRight.NbPoints() == 1),
                    ExcMessage("Transom edges don't possess a single intersection with horizontal plane!"));
    cout<<"Transom Point Left: "<<Pnt(IntersectorLeft.Point(1))<<endl;
    cout<<"Transom Point Right: "<<Pnt(IntersectorRight.Point(1))<<endl;
    cout<<"Transom Point Center: "<<PointCenterTransom<<endl;
    transomLeft = IntersectorLeft.Point(1);
    transomRight = IntersectorRight.Point(1);
    gpFront = Intersector.Point(1);
    gpBack = Pnt(PointCenterTransom);
    boatWetLength = fabs((transomLeft.X()+transomRight.X())/2.0-gpFront.X());
    Point<3> junction = (Pnt(transomLeft) + Pnt(transomRight))/2.0;
    junction(1) = 0.0;
    junction(0) += Pnt(transomLeft).distance(Pnt(transomRight))/slope;
    cout<<junction<<endl;
    Point<3> far(50*boatWetLength,0.0,0.0);
    GC_MakeSegment first_wake(Pnt(far), Pnt(junction));
    GC_MakeSegment second_wake_left(Pnt(junction), transomLeft);
    GC_MakeSegment second_wake_right(Pnt(junction), transomRight);
    Handle(Geom_Curve) first_wake_curve = first_wake.Value();
    Handle(Geom_Curve) second_wake_left_curve = second_wake_left.Value();
    Handle(Geom_Curve) second_wake_right_curve = second_wake_right.Value();
    Handle(Geom_BoundedCurve) first_wake_bcurve = Handle(Geom_BoundedCurve)::DownCast(first_wake_curve);
    Handle(Geom_BoundedCurve) second_wake_left_bcurve = Handle(Geom_BoundedCurve)::DownCast(second_wake_left_curve);
    Handle(Geom_BoundedCurve) second_wake_right_bcurve = Handle(Geom_BoundedCurve)::DownCast(second_wake_right_curve);
    bool check = false;
    GeomConvert_CompCurveToBSplineCurve convert_left_wake_bspline(first_wake_bcurve,Convert_TgtThetaOver2);
    GeomConvert_CompCurveToBSplineCurve convert_right_wake_bspline(first_wake_bcurve,Convert_TgtThetaOver2);
    check = convert_left_wake_bspline.Add(second_wake_left_bcurve,1e-7,0,1,0);
    if (check == false)
       cout<<"Failed joining left wake line"<<endl;
    left_wake_bspline = convert_left_wake_bspline.BSplineCurve();
    check = convert_right_wake_bspline.Add(second_wake_right_bcurve,1e-7,0,1,0);
    if (check == false)
       cout<<"Failed joining left wake line"<<endl;
    right_wake_bspline = convert_right_wake_bspline.BSplineCurve();
    }
  else
    {
    AssertThrow((npoints != 1) && (npoints != 2),
                 ExcMessage("Keel and has no intersection or too many intersections with horizontal plane!"));
    }
				   //we define the keel arclength and normal projection 
  wake_line_left = new ArclengthProjection(BRepBuilderAPI_MakeEdge(left_wake_bspline));
  wake_line_right = new ArclengthProjection(BRepBuilderAPI_MakeEdge(right_wake_bspline));

  cout<<"gpFront: "<<Pnt(gpFront)<<endl;
  cout<<"gpBack:"<<Pnt(gpBack)<<endl;
  cout<<"Boat Wet Lenght Lpp: "<<boatWetLength<<endl;

				   //we define the hull surface normal projections 
  boat_surface_right = new NormalProjection<2>(sh);
  boat_surface_left = new NormalProjection<2>(refl_sh);
				   //we define the hull surface y direction projections 
  boat_water_line_right = new AxisProjection(sh, Point<3>(0,1,0),1e-7,1e-2*boatWetLength);
  boat_water_line_left = new AxisProjection(refl_sh, Point<3>(0,-1,0),1e-7,1e-2*boatWetLength);
				   //we define the corresponding waterline arclength and projection 
  water_line_right = new ArclengthProjection(right_undist_water_line,1e-5*boatWetLength);
  water_line_left = new ArclengthProjection(left_undist_water_line,1e-5*boatWetLength);
				   //we define the keel arclength and normal projection 
  boat_keel = new ArclengthProjection(keel_edge,1e-5*boatWetLength);
  boat_keel_norm = new NormalProjection<1>(keel_edge);
				   //we define the left and right transom arclength projection 
  boat_transom_left = new ArclengthProjection(left_transom_edge,1e-5*boatWetLength);
  boat_transom_right = new ArclengthProjection(right_transom_edge,1e-5*boatWetLength);


  PointFrontTop = Pnt(gpFront);

  PointBackTop= Pnt(gpBack);				  

  PointBackBot = boat_keel->arclength_projection(PointBackTop, PointFrontTop,
  						 back_keel_length); 
  	  
  PointFrontBot = boat_keel->arclength_projection(PointBackTop, PointFrontTop,
  						  1.0-front_keel_length);

  PointLeftTransom = Pnt(transomLeft);

  PointRightTransom = Pnt(transomRight);
  
  PointMidBot = boat_keel->arclength_projection(PointBackTop, PointFrontTop,
  						middle_keel_length);
  

  boat_water_line_right->axis_projection(PointMidTop, Point<3>(0.0,0.0,0.0));



}

Point<3> BoatModel::compute_hydrostatic_force(const double &sink)
{
  double rho = 998.1;
  double g = 9.81;  

  double z_zero = sink;
  Point<3> hydrostatic_force(0.0,0.0,0.0);
  double wet_surface = 0.0;
  // we will need a quadrature
  QGauss<2> quad(300);
  // we now loop over all the CAD faces of sh
  TopExp_Explorer faceExplorer(sh, TopAbs_FACE);
  TopoDS_Face face;
  gp_Pnt Pint;
  std::vector<bool> to_be_changed;
  unsigned int face_count = 0;
  while(faceExplorer.More())
       {
       face = TopoDS::Face(faceExplorer.Current());
       Standard_Real umin, umax, vmin, vmax;
       BRepTools::UVBounds(face, umin, umax, vmin, vmax);          // create surface
       Handle(Geom_Surface) surf=BRep_Tool::Surface(face);          // get surface associated with face
       // creating a 2d triangulation here with a single cell here with points located on the umin, umax, vmin, vmax boundaries
       Triangulation<2,2> ref_triangulation;

       std::vector<Point<2> > ref_vertices;
       std::vector<CellData<2> > ref_cells;
       SubCellData ref_subcelldata;

       ref_vertices.resize(4);
       ref_vertices[0](0)=umin; ref_vertices[0](1)=vmin;
       ref_vertices[1](0)=umax; ref_vertices[1](1)=vmin;
       ref_vertices[2](0)=umin; ref_vertices[2](1)=vmax;
       ref_vertices[3](0)=umax; ref_vertices[3](1)=vmax;

       ref_cells.resize(1);

       ref_cells[0].vertices[0]=0; ref_cells[0].vertices[1]=1; ref_cells[0].vertices[2]=3; ref_cells[0].vertices[3]=2;
       ref_cells[0].material_id = 1;

       GridTools::delete_unused_vertices (ref_vertices, ref_cells, ref_subcelldata);
       GridReordering<2,2>::reorder_cells (ref_cells);

       ref_triangulation.create_triangulation_compatibility(ref_vertices, ref_cells, ref_subcelldata );

       // with this triangulation we create a and a FE and a FEValues (the jacobian will account for
       // transformation from [0,1]x[0,1] to [umin,umax]x[vmin,vmax], we'll have to add the other part)

       FE_Q<2,2> fe(1);
       DoFHandler<2,2> ref_dh(ref_triangulation);
       ref_dh.distribute_dofs(fe);

       FEValues<2,2> ref_fe_v(StaticMappingQ1<2,2>::mapping, fe,quad,
   		              update_values |
		              update_quadrature_points |
		              update_JxW_values);
       
       ref_fe_v.reinit(ref_dh.begin_active());
       const unsigned int n_q_points = ref_fe_v.n_quadrature_points;
       
       // we now perform the actual pressure integral over the surface patch
       const std::vector<Point<2> > &quad_nodes = ref_fe_v.get_quadrature_points();
       for (unsigned int i=0; i<n_q_points; ++i)
           {
           gp_Pnt q;
           gp_Vec Du, Dv, normal;
           surf->D1(quad_nodes[i](0),quad_nodes[i](1), q, Du, Dv);
           double jacobian = Du.CrossMagnitude(Dv);
           normal = Du^ Dv/jacobian;
           Point<3> Normal(normal.X(),normal.Y(),normal.Z());
           // adjusting normal orientation
           if (face.Orientation()==TopAbs_REVERSED)
              {
              Normal*=-1.0;
              }
           hydrostatic_force += (rho*g*fmax(z_zero-q.Z(),0.0))*Normal*jacobian*ref_fe_v.JxW(i);
           if (q.Z() <= z_zero)
              wet_surface += jacobian*ref_fe_v.JxW(i);
           //cout<<"q("<<Pnt(q)<<")  p:"<<(rho*g*fmax(z_zero-q.Z(),0.0))<<"  n("<<Normal<<")"<<endl;
           }

       //GeomLProp_SLProps props(surf, (umin+umax)/2.0, (vmin+vmax)/2.0, 1, 0.01);          // get surface normal
       //gp_Dir norm=props.Normal();   

       faceExplorer.Next();
       ++face_count;
      //cout<<"Face count: "<<face_count<<endl;
      }

  // we now instead loop over all the CAD faces of refl_sh
  TopExp_Explorer reflFaceExplorer(refl_sh, TopAbs_FACE);

  face_count = 0;
  while(reflFaceExplorer.More())
       {
       //Point<3> face_force(0.0,0.0,0.0);
       face = TopoDS::Face(reflFaceExplorer.Current());
       Standard_Real umin, umax, vmin, vmax;
       BRepTools::UVBounds(face, umin, umax, vmin, vmax);          // create surface
       Handle(Geom_Surface) surf=BRep_Tool::Surface(face);          // get surface associated with face
       // creating a 2d triangulation here with a single cell here with points located on the umin, umax, vmin, vmax boundaries
       Triangulation<2,2> ref_triangulation;

       std::vector<Point<2> > ref_vertices;
       std::vector<CellData<2> > ref_cells;
       SubCellData ref_subcelldata;

       ref_vertices.resize(4);
       ref_vertices[0](0)=umin; ref_vertices[0](1)=vmin;
       ref_vertices[1](0)=umax; ref_vertices[1](1)=vmin;
       ref_vertices[2](0)=umin; ref_vertices[2](1)=vmax;
       ref_vertices[3](0)=umax; ref_vertices[3](1)=vmax;

       ref_cells.resize(1);

       ref_cells[0].vertices[0]=0; ref_cells[0].vertices[1]=1; ref_cells[0].vertices[2]=3; ref_cells[0].vertices[3]=2;
       ref_cells[0].material_id = 1;

       GridTools::delete_unused_vertices (ref_vertices, ref_cells, ref_subcelldata);
       GridReordering<2,2>::reorder_cells (ref_cells);

       ref_triangulation.create_triangulation_compatibility(ref_vertices, ref_cells, ref_subcelldata );

       // with this triangulation we create a and a FE and a FEValues (the jacobian will account for
       // transformation from [0,1]x[0,1] to [umin,umax]x[vmin,vmax], we'll have to add the other part)

       FE_Q<2,2> fe(1);
       DoFHandler<2,2> ref_dh(ref_triangulation);
       ref_dh.distribute_dofs(fe);

       FEValues<2,2> ref_fe_v(StaticMappingQ1<2,2>::mapping, fe,quad,
   		              update_values |
		              update_quadrature_points |
		              update_JxW_values);
       
       ref_fe_v.reinit(ref_dh.begin_active());
       const unsigned int n_q_points = ref_fe_v.n_quadrature_points;
       
       // we now perform the actual pressure integral over the surface patch
       const std::vector<Point<2> > &quad_nodes = ref_fe_v.get_quadrature_points();
       for (unsigned int i=0; i<n_q_points; ++i)
           {
           gp_Pnt q;
           gp_Vec Du, Dv, normal;
           surf->D1(quad_nodes[i](0),quad_nodes[i](1), q, Du, Dv);
           double jacobian = Du.CrossMagnitude(Dv);
           normal = (Du^Dv)/jacobian;
           Point<3> Normal(normal.X(),normal.Y(),normal.Z());
           // adjusting normal orientation
           if (face.Orientation()==TopAbs_REVERSED)
              {
              Normal*=-1.0;
              }
           hydrostatic_force += (rho*g*fmax(z_zero-q.Z(),0.0))*Normal*jacobian*ref_fe_v.JxW(i);
           if (q.Z() <= z_zero)
              wet_surface += jacobian*ref_fe_v.JxW(i);
           //face_force += (rho*g*fmax(z_zero-q.Z(),0.0))*Normal*jacobian*ref_fe_v.JxW(i); 
           //cout<<"q("<<Pnt(q)<<")  p:"<<(rho*g*fmax(z_zero-q.Z(),0.0))<<"  n("<<Normal<<")"<<endl;
           }

       //GeomLProp_SLProps props(surf, (umin+umax)/2.0, (vmin+vmax)/2.0, 1, 0.01);          // get surface normal
       //gp_Dir norm=props.Normal();   

       reflFaceExplorer.Next();
       ++face_count;
      //cout<<"Face count: "<<face_count<<"  Face force: "<<face_force<<endl;
      }


      cout<<"Current hydrostatic force: "<<hydrostatic_force<<"   Current wet surface: "<<wet_surface<<endl;
      return hydrostatic_force;
}


Point<3> BoatModel::compute_hydrostatic_moment(const double &sink)
{
  double rho = 998.1;
  double g = 9.81;  

  double z_zero = sink;
  Point<3> hydrostatic_moment(0.0,0.0,0.0);
  double wet_surface = 0.0;
  // we will need a quadrature
  QGauss<2> quad(300);
  // we now loop over all the CAD faces of sh
  TopExp_Explorer faceExplorer(sh, TopAbs_FACE);
  TopoDS_Face face;
  gp_Pnt Pint;
  std::vector<bool> to_be_changed;
  unsigned int face_count = 0;
  while(faceExplorer.More())
       {
       face = TopoDS::Face(faceExplorer.Current());
       Standard_Real umin, umax, vmin, vmax;
       BRepTools::UVBounds(face, umin, umax, vmin, vmax);          // create surface
       Handle(Geom_Surface) surf=BRep_Tool::Surface(face);          // get surface associated with face
       // creating a 2d triangulation here with a single cell here with points located on the umin, umax, vmin, vmax boundaries
       Triangulation<2,2> ref_triangulation;

       std::vector<Point<2> > ref_vertices;
       std::vector<CellData<2> > ref_cells;
       SubCellData ref_subcelldata;

       ref_vertices.resize(4);
       ref_vertices[0](0)=umin; ref_vertices[0](1)=vmin;
       ref_vertices[1](0)=umax; ref_vertices[1](1)=vmin;
       ref_vertices[2](0)=umin; ref_vertices[2](1)=vmax;
       ref_vertices[3](0)=umax; ref_vertices[3](1)=vmax;

       ref_cells.resize(1);

       ref_cells[0].vertices[0]=0; ref_cells[0].vertices[1]=1; ref_cells[0].vertices[2]=3; ref_cells[0].vertices[3]=2;
       ref_cells[0].material_id = 1;

       GridTools::delete_unused_vertices (ref_vertices, ref_cells, ref_subcelldata);
       GridReordering<2,2>::reorder_cells (ref_cells);

       ref_triangulation.create_triangulation_compatibility(ref_vertices, ref_cells, ref_subcelldata );

       // with this triangulation we create a and a FE and a FEValues (the jacobian will account for
       // transformation from [0,1]x[0,1] to [umin,umax]x[vmin,vmax], we'll have to add the other part)

       FE_Q<2,2> fe(1);
       DoFHandler<2,2> ref_dh(ref_triangulation);
       ref_dh.distribute_dofs(fe);

       FEValues<2,2> ref_fe_v(StaticMappingQ1<2,2>::mapping, fe,quad,
   		              update_values |
		              update_quadrature_points |
		              update_JxW_values);
       
       ref_fe_v.reinit(ref_dh.begin_active());
       const unsigned int n_q_points = ref_fe_v.n_quadrature_points;
       
       // we now perform the actual pressure integral over the surface patch
       const std::vector<Point<2> > &quad_nodes = ref_fe_v.get_quadrature_points();
       for (unsigned int i=0; i<n_q_points; ++i)
           {
           gp_Pnt q;
           gp_Vec Du, Dv, normal;
           surf->D1(quad_nodes[i](0),quad_nodes[i](1), q, Du, Dv);
           double jacobian = Du.CrossMagnitude(Dv);
           normal = Du^ Dv/jacobian;
           Point<3> Normal(normal.X(),normal.Y(),normal.Z());
           // adjusting normal orientation
           if (face.Orientation()==TopAbs_REVERSED)
              {
              Normal*=-1.0;
              }
           Point<3> Q(q.X(),q.Y(),q.Z());
           Point<3> Kross(Q(1)*Normal(2)-Q(2)*Normal(1),Q(2)*Normal(0)-Q(0)*Normal(2),Q(0)*Normal(1)-Q(1)*Normal(0));
           hydrostatic_moment += (rho*g*fmax(z_zero-q.Z(),0.0))*Kross*jacobian*ref_fe_v.JxW(i);
           //cout<<"q("<<Pnt(q)<<")  p:"<<(rho*g*fmax(z_zero-q.Z(),0.0))<<"  n("<<Normal<<")"<<endl;
           }

       //GeomLProp_SLProps props(surf, (umin+umax)/2.0, (vmin+vmax)/2.0, 1, 0.01);          // get surface normal
       //gp_Dir norm=props.Normal();   

       faceExplorer.Next();
       ++face_count;
      //cout<<"Face count: "<<face_count<<endl;
      }

  // we now instead loop over all the CAD faces of refl_sh
  TopExp_Explorer reflFaceExplorer(refl_sh, TopAbs_FACE);

  face_count = 0;
  while(reflFaceExplorer.More())
       {
       //Point<3> face_force(0.0,0.0,0.0);
       face = TopoDS::Face(reflFaceExplorer.Current());
       Standard_Real umin, umax, vmin, vmax;
       BRepTools::UVBounds(face, umin, umax, vmin, vmax);          // create surface
       Handle(Geom_Surface) surf=BRep_Tool::Surface(face);          // get surface associated with face
       // creating a 2d triangulation here with a single cell here with points located on the umin, umax, vmin, vmax boundaries
       Triangulation<2,2> ref_triangulation;

       std::vector<Point<2> > ref_vertices;
       std::vector<CellData<2> > ref_cells;
       SubCellData ref_subcelldata;

       ref_vertices.resize(4);
       ref_vertices[0](0)=umin; ref_vertices[0](1)=vmin;
       ref_vertices[1](0)=umax; ref_vertices[1](1)=vmin;
       ref_vertices[2](0)=umin; ref_vertices[2](1)=vmax;
       ref_vertices[3](0)=umax; ref_vertices[3](1)=vmax;

       ref_cells.resize(1);

       ref_cells[0].vertices[0]=0; ref_cells[0].vertices[1]=1; ref_cells[0].vertices[2]=3; ref_cells[0].vertices[3]=2;
       ref_cells[0].material_id = 1;

       GridTools::delete_unused_vertices (ref_vertices, ref_cells, ref_subcelldata);
       GridReordering<2,2>::reorder_cells (ref_cells);

       ref_triangulation.create_triangulation_compatibility(ref_vertices, ref_cells, ref_subcelldata );

       // with this triangulation we create a and a FE and a FEValues (the jacobian will account for
       // transformation from [0,1]x[0,1] to [umin,umax]x[vmin,vmax], we'll have to add the other part)

       FE_Q<2,2> fe(1);
       DoFHandler<2,2> ref_dh(ref_triangulation);
       ref_dh.distribute_dofs(fe);

       FEValues<2,2> ref_fe_v(StaticMappingQ1<2,2>::mapping, fe,quad,
   		              update_values |
		              update_quadrature_points |
		              update_JxW_values);
       
       ref_fe_v.reinit(ref_dh.begin_active());
       const unsigned int n_q_points = ref_fe_v.n_quadrature_points;
       
       // we now perform the actual pressure integral over the surface patch
       const std::vector<Point<2> > &quad_nodes = ref_fe_v.get_quadrature_points();
       for (unsigned int i=0; i<n_q_points; ++i)
           {
           gp_Pnt q;
           gp_Vec Du, Dv, normal;
           surf->D1(quad_nodes[i](0),quad_nodes[i](1), q, Du, Dv);
           double jacobian = Du.CrossMagnitude(Dv);
           normal = (Du^Dv)/jacobian;
           Point<3> Normal(normal.X(),normal.Y(),normal.Z());

           // adjusting normal orientation
           if (face.Orientation()==TopAbs_REVERSED)
              {
              Normal*=-1.0;
              }
           Point<3> Q(q.X(),q.Y(),q.Z());
           Point<3> Kross(Q(1)*Normal(2)-Q(2)*Normal(1),Q(2)*Normal(0)-Q(0)*Normal(2),Q(0)*Normal(1)-Q(1)*Normal(0));
           hydrostatic_moment += (rho*g*fmax(z_zero-q.Z(),0.0))*Kross*jacobian*ref_fe_v.JxW(i);

           //face_force += (rho*g*fmax(z_zero-q.Z(),0.0))*Normal*jacobian*ref_fe_v.JxW(i); 
           //cout<<"q("<<Pnt(q)<<")  p:"<<(rho*g*fmax(z_zero-q.Z(),0.0))<<"  n("<<Normal<<")"<<endl;
           }

       //GeomLProp_SLProps props(surf, (umin+umax)/2.0, (vmin+vmax)/2.0, 1, 0.01);          // get surface normal
       //gp_Dir norm=props.Normal();   

       reflFaceExplorer.Next();
       ++face_count;
      //cout<<"Face count: "<<face_count<<"  Face force: "<<face_force<<endl;
      }


      cout<<"Current hydrostatic moment: "<<hydrostatic_moment<<endl;
      return hydrostatic_moment;
}


void BoatModel::compute_hydrostatic_sink(double &sink, const double &weight)
{
double delta_sink = 0.001;
double delta_z;
Point<3> hydrostatic_force;
Point<3> hydrostatic_force_n_minus_2 = compute_hydrostatic_force(0.0);
Point<3> hydrostatic_force_n_minus_1 = compute_hydrostatic_force(delta_sink);
double z;
double z_n_minus_1 = delta_sink;
double z_n_minus_2 = 0.0;

while (fabs(hydrostatic_force(2)-weight)/weight > 1e-7)
{
double slope = (hydrostatic_force_n_minus_1(2)-hydrostatic_force_n_minus_2(2))/(z_n_minus_1-z_n_minus_2);
z = z_n_minus_1 - (hydrostatic_force_n_minus_1(2)-weight)/slope; 
z_n_minus_2 = z_n_minus_1;
z_n_minus_1 = z;
hydrostatic_force = compute_hydrostatic_force(z);
hydrostatic_force_n_minus_2 = hydrostatic_force_n_minus_1;
hydrostatic_force_n_minus_1 = hydrostatic_force;
//cout<<"z: "<<z<<"   hydrostatic force: "<<hydrostatic_force<<"   check: "<<fabs(hydrostatic_force(2)-weight)/weight<<endl;
}

sink = z;



cout<<"z: "<<z<<"  hydrostatic force (out): "<<hydrostatic_force<<"  (weight: "<<weight<<")"<<endl;
}

class BoatModel;
