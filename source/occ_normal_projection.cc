#include "occ_normal_projection.h"
#include "occ_utilities.h"

#include <iostream>
#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Face.hxx>
#include <BRepTools.hxx>
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
#include <BRepBuilderAPI_MakeEdge.hxx>
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
#include <ShapeAnalysis_Surface.hxx>
#include <GeomLProp_SLProps.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <TColgp_Array1OfPnt2d.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <vector>
#define FALSE 0
#define TRUE 1


#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>

namespace OpenCascade 
{
  
  template <int dim>
  NormalProjection<dim>::NormalProjection(const TopoDS_Shape &sh) :
  sh(sh)
  {
    if(dim == 1) 
      {
					 // Check that we have at least
					 // one edge.
	TopExp_Explorer edgeExplorer(sh , TopAbs_EDGE);
	unsigned int n_edges = 0;
	while(edgeExplorer.More())
	  {
	    n_edges++;
	    edgeExplorer.Next();
	  }
	AssertThrow(n_edges > 0, ExcMessage("We have no edges to process"));
      }
    else if(dim == 2)
      {
					 // Check that we have at least
					 // one face.
	TopExp_Explorer faceExplorer(sh , TopAbs_FACE);
	unsigned int n_faces = 0;
	while(faceExplorer.More())
	  {
	    n_faces++;
	    faceExplorer.Next();
	  }
	AssertThrow(n_faces > 0, ExcMessage("We have no faces to process"));
      }
    else
      AssertThrow(false, ExcNotImplemented());
  }

  
  template <>
  void NormalProjection<1>::normal_projection
  (Point<3> &projection, const Point<3> &origin) const 
  {
    TopLoc_Location L = sh.Location();
    TopLoc_Location L_inv = L.Inverted();	
				     // translating original
				     // Point<dim> to gp point
    gp_Pnt P0(origin(0),origin(1),origin(2));
    P0.Transform(L_inv);
				     // destination point
    gp_Pnt Pproj(0.0,0.0,0.0);
				     // we prepare now the surface
				     // for the projection we get
				     // the whole shape from the
				     // iges model

// and here we loop on the faces of the shape
    TopExp_Explorer edgeExplorer(sh , TopAbs_EDGE);
    unsigned int tot_proj_points = 0;
    double minDistance = 1e7;
    TopoDS_Edge edge;
    while(edgeExplorer.More())
      {
	edge = TopoDS::Edge(edgeExplorer.Current());
        edge.Location(L);
        BRepAdaptor_Curve AC(edge);
	Standard_Real First;
	Standard_Real Last;
					 // the projection
					 // function needs a
					 // surface, so we obtain
					 // the surface upon which
					 // the face is defined
	Handle(Geom_Curve) CurveToProj = BRep_Tool::Curve(edge,L,First,Last);

	TopoDS_Shape shnew = BRepBuilderAPI_MakeEdge(CurveToProj);
	    
	GeomAPI_ProjectPointOnCurve Proj(P0,CurveToProj);
	unsigned int num_proj_points = Proj.NbPoints();
	tot_proj_points+=num_proj_points;
	if ((num_proj_points > 0) && (Proj.LowerDistance() < minDistance))
	  {
	    Pproj = Proj.NearestPoint();
            Pproj.Transform(L);
            minDistance = Proj.LowerDistance();
	  }
	edgeExplorer.Next();
      }

    Assert(tot_proj_points > 0,
	   ExcMessage("Point projection on curve in normal direction does not exist"));

// translating destination point
    projection(0) = Pproj.X();
    projection(1) = Pproj.Y();
    projection(2) = Pproj.Z();

  }

  template <>
  void NormalProjection<2>::normal_projection
  (Point<3> &projection, const Point<3> &origin) const 
  {
  TopLoc_Location L = sh.Location();
  TopLoc_Location L_inv = L.Inverted();	
				     // translating original
				     // Point<dim> to gp point
  gp_Pnt P0(origin(0),origin(1),origin(2));
  P0.Transform(L_inv);
				     // destination point
  gp_Pnt Pproj(0.0,0.0,0.0);
				     // we prepare now the surface
				     // for the projection we get
				     // the whole shape from the
				     // iges model


  // and here we loop on the faces of the shape
  unsigned int face_count=0;
  TopExp_Explorer faceExplorer(sh, TopAbs_FACE);
  double minDistance = 1e7;
  TopoDS_Face face;
  gp_Pnt face_proj(0.0,0.0,0.0);

  while(faceExplorer.More())
      {
      face = TopoDS::Face(faceExplorer.Current());
					 // the projection
					 // function needs a
					 // surface, so we obtain
					 // the surface upon which
					 // the face is defined
      Handle(Geom_Surface) SurfToProj = BRep_Tool::Surface(face);


      ShapeAnalysis_Surface projector(SurfToProj);
      gp_Pnt2d proj_params = projector.ValueOfUV(P0, 1e-7);
      SurfToProj->D0(proj_params.X(),proj_params.Y(),face_proj);
      if (Pnt(face_proj).distance(origin) < minDistance)
         {
         minDistance = Pnt(face_proj).distance(origin);
         Pproj = face_proj;
         }
      faceExplorer.Next();
      ++face_count;
      }
  Pproj.Transform(L);


// translating destination point
    projection = Pnt(Pproj);


}



  template <>
  void NormalProjection<1>::normal_projection_and_diff_forms(Point<3> &/*projection*/,
                                                             Point<3> &/*normal*/,
                                                             double &/*mean_curvature*/,
			                                     const Point<3> &/*origin*/) const
  {
      Assert(0 > 1,
	   ExcMessage("Method is not implemented for dim = 1"));
  }


  template <>
  void NormalProjection<2>::normal_projection_and_diff_forms(Point<3> &projection,
                                                             Point<3> &normal,
                                                             double &mean_curvature,
			                                     const Point<3> &origin) const
{
TopLoc_Location L = sh.Location();
TopLoc_Location L_inv = L.Inverted();	
				     // translating original
				     // Point<dim> to gp point
gp_Pnt P0 = Pnt(origin);
//P0.Transform(L_inv);
				     // destination point
gp_Pnt Pproj(0.0,0.0,0.0);

				     // we prepare now the surface
				     // for the projection we get
				     // the whole shape from the
				     // iges model


// and here we loop on the faces of the shape
unsigned int face_count=0;
TopExp_Explorer faceExplorer(sh, TopAbs_FACE);
double minDistance = 1e7;
TopoDS_Face face;
gp_Pnt face_proj(0.0,0.0,0.0);
gp_Dir Normal;
Standard_Real Mean_Curvature = 0;

while(faceExplorer.More())
    {
    face = TopoDS::Face(faceExplorer.Current());
    face.Location(L);
    BRepAdaptor_Surface AF(face);
    //cout<<"** "<<Pnt(AF.Value(AF.FirstUParameter(),AF.FirstVParameter()))<<endl;
    //cout<<"** "<<Pnt(AF.Value(AF.LastUParameter(),AF.FirstVParameter()))<<endl;
    //cout<<"** "<<Pnt(AF.Value(AF.LastUParameter(),AF.LastVParameter()))<<endl;
    //cout<<"** "<<Pnt(AF.Value(AF.FirstUParameter(),AF.LastVParameter()))<<endl;
					 // the projection
					 // function needs a
					 // surface, so we obtain
					 // the surface upon which
					 // the face is defined
    
    
    Handle(Geom_Surface) SurfToProj = BRep_Tool::Surface(face);
    Standard_Real uufirst,uulast,vvfirst,vvlast;
    SurfToProj->Bounds (uufirst, uulast, vvfirst, vvlast);
    //cout<<"*#* "<<Pnt(SurfToProj->Value(uufirst,vvfirst))<<endl;
    //cout<<"*#* "<<Pnt(SurfToProj->Value(uulast,vvfirst))<<endl;
    //cout<<"*#* "<<Pnt(SurfToProj->Value(uulast,vvlast))<<endl;
    //cout<<"*#* "<<Pnt(SurfToProj->Value(uufirst,vvlast))<<endl;
    //cout<<Pnt(P0)<<endl;
    ShapeAnalysis_Surface projector(SurfToProj);
    if (projector.IsDegenerated(P0,1e-7) )
       {
       //cout<<"Gotcha!"<<endl;
       Standard_Real ufirst,ulast,vfirst,vlast;
       projector.Bounds (ufirst, ulast, vfirst, vlast);
       //cout<<"Ubound "<<ufirst<<" "<<ulast<<endl;
       //cout<<"Vbound "<<vfirst<<" "<<vlast<<endl;
       gp_Pnt Pmid = SurfToProj->Value((ufirst+ulast)/2.0, (vfirst+vlast)/2.0);
        
       gp_Pnt P0Mod(P0.X()+(Pmid.X()-P0.X())/1000.0,
                    P0.Y()+(Pmid.Y()-P0.Y())/1000.0,
                    P0.Z()+(Pmid.Z()-P0.Z())/1000.0);
       P0=P0Mod;
       //cout<<"Point "<<Pnt(P0)<<endl;
       }
    gp_Pnt2d proj_params = projector.ValueOfUV(P0, 1e-7);
    SurfToProj->D0(proj_params.X(),proj_params.Y(),face_proj);
    if (Pnt(face_proj).distance(origin) < minDistance)
       {
       minDistance = Pnt(face_proj).distance(origin);
       Pproj = face_proj;

       GeomLProp_SLProps props(SurfToProj, proj_params.X(), proj_params.Y(),1, 1e-7);
       Normal = props.Normal();
       Mean_Curvature = props.MeanCurvature();

       // adjusting normal orientation
       if (face.Orientation()==TopAbs_REVERSED)
          {
          Normal.Reverse();
          Mean_Curvature*=-1;
          }
       }
    faceExplorer.Next();
    ++face_count;
    }
    //Pproj.Transform(L);
    
// translating destination point
    projection = Pnt(Pproj);

// translating normal vector
    //Normal.Transform(L);
    normal(0) = Normal.X();
    normal(1) = Normal.Y();
    normal(2) = Normal.Z();
    // translating mean curvature
    //cout<<"????? "<<normal<<endl;

    mean_curvature = double(Mean_Curvature);

  }

  template <int dim>
  Point<3> NormalProjection<dim>::get_new_point_on_line
  (const Triangulation< 2,3 >::line_iterator &line) const
  {
    Point<3> projected_point;
    Point<3> source_point = StraightBoundary<2,3>::get_new_point_on_line(line);
    normal_projection(projected_point, source_point);
    return projected_point;
  }


  template <int dim>
  Point<3> NormalProjection<dim>::get_new_point_on_quad
  (const Triangulation< 2,3 >::quad_iterator &quad) const
  {
    AssertThrow(dim == 2, ExcMessage("Cannot project on quad on this"
				     " geometry type."));
    Point<3> proj_point;
    Point<3> init_point =
      StraightBoundary<2,3>::get_new_point_on_quad(quad);
    normal_projection(proj_point, init_point);
    return proj_point;
  }


  
  template class NormalProjection<1>;
  template class NormalProjection<2>;

}
