#include "occ_arclength_projection.h"


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
#include <ShapeAnalysis_Curve.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <vector>
#define FALSE 0
#define TRUE 1


#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>


namespace OpenCascade
{

  ArclengthProjection::ArclengthProjection(const TopoDS_Shape &sh,
                                           double tolerance) :
    sh(sh),
    tolerance(tolerance)
  {
    // Check that we have at most
    // one edge.
    TopExp_Explorer edgeExplorer(sh , TopAbs_EDGE);
    unsigned int n_edges = 0;
    while (edgeExplorer.More())
      {
        n_edges++;
        edgeExplorer.Next();
      }
    AssertThrow(n_edges == 1, ExcMessage("We can do this only "
                                         "on a single edge: "
                                         "split your geometry."));

  }

  Point<3> ArclengthProjection::arclength_projection
  (const Point<3> &p0, const Point<3> &p1, const double distance) const
  {

    cout<<p0<<"  --"<<distance<<"--  "<<p1<<endl;
    Point<3> projected_point;

    Assert((0. <= distance) &&
           (distance <= 1.),
           ExcMessage("Distance should be between 0. and 1."));

    gp_Pnt P0(p0(0),p0(1),p0(2));
    gp_Pnt P1(p1(0),p1(1),p1(2));

    TopExp_Explorer edgeExplorer(sh , TopAbs_EDGE);
    TopoDS_Edge edge;
    edgeExplorer.More();
    edge = TopoDS::Edge(edgeExplorer.Current());
    TopLoc_Location L = sh.Location();
    edge.Location(L);
    Standard_Real First;
    Standard_Real Last;


    // the projection function
    // needs a curve, so we
    // obtain the curve upon
    // which the edge is defined
    BRepAdaptor_Curve gg_curve(edge);
    First = gg_curve.FirstParameter();
    Last = gg_curve.LastParameter();
    //Handle(Geom_Curve) g_curve = BRep_Tool::Curve(edge,L,First,Last);


    //gp_Trsf L_transformation = this->sh.Location().Transformation();
    //TopLoc_Location L_inv = L.Inverted();
    //gp_Trsf L_inv_transformation = L_inv.Transformation();
    //P0.Transform(L_inv_transformation);
    //P1.Transform(L_inv_transformation);
    Standard_Real t0;
    Standard_Real t1;
    Standard_Real t2;


    gp_Pnt proj;
    ShapeAnalysis_Curve curve_analysis;

    double off = curve_analysis.Project(gg_curve, P0, tolerance, proj, t0, Standard_True);
    AssertThrow( (off < 1000*tolerance), ExcMessage("Point of the edge to be refined is not on curve."));
    off = curve_analysis.Project(gg_curve, P1, tolerance, proj, t1, Standard_True);
    AssertThrow( (off < 1000*tolerance), ExcMessage("Point of the edge to be refined is not on curve."));

    //cout<<"First: "<<First<<"  Last: "<<Last<<endl;
    //cout<<"P0: "<<p0<<" ("<<t0<<")  vs  P1: "<<p1<<" ("<<t1<<")"<<endl;

    //GeomLib_Tool tool;
    //tool.Parameter(gg_curve, P0, tolerance, t0);
    //tool.Parameter(gg_curve, P1, tolerance, t1);

    //Check that we are in the right
    //range.
    AssertThrow((First-1000*tolerance <= t0) &&
                (t0 <= Last+1000*tolerance),
                ExcMessage("Parameter 1 is out of range!"));
    AssertThrow((First-1000*tolerance <= t1) &&
                (t1 <= Last+1000*tolerance),
                ExcMessage("Parameter 2 is out of range!"));


    // we now get the mean
    // curvilinear distance point
    //cout<<endl;
    //cout<<"Points "<<p0<<" vs "<<p1<<endl;
    //cout<<"Params "<<t0<<" vs "<<t1<<endl;
    double direction = t1 > t0 ? 1.0 : -1.0;
    //cout<<"Direction "<<direction<<endl;
    //GeomAdaptor_Curve AC(g_curve);
    Standard_Real arcLength = GCPnts_AbscissaPoint::Length(gg_curve,t0,t1);
    //cout<<"arcLength "<<arcLength<<endl;
    GCPnts_AbscissaPoint AP(gg_curve, direction*arcLength*distance, t0);
    //cout<<"direction*arcLength*distance "<<direction*arcLength*distance<<endl;
    t2 = AP.Parameter();

    AssertThrow((First <= t2) &&
                (t2 <= Last),
                ExcMessage("Parameter 3 is out of range!"));


    gp_Pnt P2 = gg_curve.Value(t2);
    cout<<"new_point "<<P2.X()<<" "<<P2.Y()<<" "<<P2.Z()<<endl;
    //P2.Transform(L_transformation);
    //cout<<"new_point re-located "<<P2.X()<<" "<<P2.Y()<<" "<<P2.Z()<<endl;

    projected_point(0) = P2.X();
    projected_point(1) = P2.Y();
    projected_point(2) = P2.Z();
    return projected_point;
  }

  Point<3> ArclengthProjection::get_new_point_on_line
  (const Triangulation< 2,3 >::line_iterator &line) const
  {
    return arclength_projection(line->vertex(0), line->vertex(1));
  }

  Point<3> ArclengthProjection::project_to_manifold(const std::vector<Point<3> > &points,
                               const Point<3 > & candidate) const
  {
    return arclength_projection(points[0], points[1]);
  }
}
