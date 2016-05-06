#include "occ_utilities.h"


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <set>
#include <IGESControl_Reader.hxx>
#include <IGESControl_Controller.hxx>
#include <IGESControl_Writer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
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
#include <GeomAPI_IntSS.hxx>
#include <Bnd_Box.hxx>
#include <gp_Trsf.hxx>
#include <gp_Ax3.hxx>
#include <gp_Pln.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <GeomConvert_CompCurveToBSplineCurve.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <TColGeom_Array1OfCurve.hxx>
#include <TColgp_HArray1OfPnt.hxx>
#include <Geom_Curve.hxx>
#include <Geom_BoundedCurve.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Geom_BSplineCurve.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <BRepAlgo_Section.hxx>
#include <boost/bind.hpp>
#include <GeomLib_Tool.hxx>
#include <TColGeom_Array2OfBezierSurface.hxx>
#include <ProjLib_ProjectOnPlane.hxx>
#include <Adaptor3d_HCurve.hxx>
#include <GeomAdaptor_HCurve.hxx>
#include <ShapeAnalysis_Surface.hxx>
#include <GeomLProp_SLProps.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <GeomLProp_CLProps.hxx>
#include <BRep_Builder.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>
#include <BRepMesh.hxx>
#include <GeomPlate_BuildPlateSurface.hxx>
#include <GeomPlate_MakeApprox.hxx>
#include <GeomAdaptor_HCurve.hxx>
#include <Adaptor3d_HCurveOnSurface.hxx>
#include <GeomPlate_CurveConstraint.hxx>
#include <TColgp_SequenceOfXY.hxx>
#include <TColgp_SequenceOfXYZ.hxx>
#include <GeomPlate_PlateG0Criterion.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <Poly_Triangulation.hxx>
#include <BRepMesh_FastDiscret.hxx>

#define FALSE 0
#define TRUE 1

#include <deal.II/base/utilities.h>
#include <deal.II/base/exceptions.h>
#include <vector>
#include <algorithm>

using namespace dealii;
using namespace std;
using namespace numbers;

namespace OpenCascade
{

  TopoDS_Shape read_IGES(string filename, double scale_factor)
  {
    IGESControl_Reader reader;
    IFSelect_ReturnStatus stat;
    stat = reader.ReadFile(filename.c_str());
    Standard_Boolean failsonly = Standard_False;
    IFSelect_PrintCount mode = IFSelect_ItemsByEntity;
    reader.PrintCheckLoad (failsonly, mode);
    Standard_Integer nIgesFaces,nTransFaces;

    Handle(TColStd_HSequenceOfTransient) myList = reader.GiveList("iges-faces");
    //selects all IGES faces in the
    //file and puts them into a list
    //called MyList,
    nIgesFaces = myList->Length();
    nTransFaces = reader.TransferList(myList);

    AssertThrow(nTransFaces > 0, ExcMessage("Read nothing from file."));

    // Handle IGES Scale here.
    gp_Pnt Origin;
    gp_Trsf scale;
    scale.SetScale (Origin, scale_factor);

    TopoDS_Shape sh = reader.OneShape();
    BRepBuilderAPI_Transform trans(sh, scale);

    return trans.Shape();   // this is the actual translation
  }

  TopoDS_Shape extract_xz_edges(const TopoDS_Shape &in_shape,
                                const double tolerance,
                                const unsigned int max_num_edges)
  {

    //// this is to loop on the edges
    TopExp_Explorer edgeExplorer(in_shape , TopAbs_EDGE);
    unsigned int numKeelEdges = 0;
    unsigned int numEdges = 0;

    TopoDS_Edge edge = TopoDS::Edge(edgeExplorer.Current());
    TopLoc_Location L;
    Standard_Real First;
    Standard_Real Last;
    gp_Pnt PIn(0.0,0.0,0.0);
    gp_Pnt PFin(0.0,0.0,0.0);
    gp_Pnt PMid(0.0,0.0,0.0);

    Handle(Geom_Curve) curveOne = BRep_Tool::Curve(edge,L,First,Last);
    std::vector< Handle(Geom_Curve) > please;

    while (edgeExplorer.More())
      {
        edge = TopoDS::Edge(edgeExplorer.Current());
        if (BRep_Tool::Degenerated(edge))
          {
          }
        else
          {
            Handle(Geom_Curve) curve = BRep_Tool::Curve(edge,L,First,Last);
            curve->D0(First,PIn);
            curve->D0(Last,PFin);
            curve->D0((First+Last)/2.0,PMid);
            //cout<<Pnt(PIn)<<" | "<<Pnt(PMid)<<" | "<<Pnt(PFin)<<" ---> ";
            //cout<<fabs(double(PIn.Y())+double(PFin.Y())+double(PMid.Y()))/3.0<<endl;
            if ((fabs(double(PIn.Y())+double(PFin.Y())+double(PMid.Y()))/3.0 < tolerance) &&
                (PIn.Distance(PMid)+PFin.Distance(PMid) > tolerance) )
              {
                cout<<numKeelEdges<<":  "<<Pnt(PIn)<<" | "<<Pnt(PMid)<<" | "<<Pnt(PFin)<<endl;
                please.push_back(curve);
                ++numKeelEdges;
              }
          }
        ++numEdges;
        edgeExplorer.Next();
      }

    numKeelEdges = please.size();
    cout<<"numKeelEdges "<<numKeelEdges<<endl;
    AssertThrow(numKeelEdges>0,
                ExcMessage("No edges on xz plane"));

    for (unsigned int i=0; i<numKeelEdges; ++i)
      {
        Handle(Geom_Curve) curvez = please[i];
        First = curvez->FirstParameter();
        Last = curvez->LastParameter();
        curvez->D0(First,PIn);
        curvez->D0(Last,PFin);
        curvez->D0((First+Last)/2.0,PMid);
        cout<<i<<": "<<Pnt(PIn)<<" | "<<Pnt(PMid)<<" | "<<Pnt(PFin)<<endl;
      }
//*/
    Handle(Geom_Curve) curve1 = please[0];
    Handle(Geom_BoundedCurve) bcurve1 = Handle(Geom_BoundedCurve)::DownCast(curve1);
    GeomConvert_CompCurveToBSplineCurve convert_keel_bspline(bcurve1,Convert_TgtThetaOver2);
    bool check = false, one_added = true, one_failed=true;
    vector<bool> added(numKeelEdges, false);
    added[0] = true;
    unsigned int added_count = 1;
    while (one_added == true)
      {

        one_added = false;
        one_failed = false;
        for (unsigned int i=1; i<numKeelEdges; ++i)
          if (added[i] == false)
            {
              cout<<"Curve "<<i<<" of "<<numKeelEdges<<" ("<<added_count<<" already added"<<")"<<endl;
              Handle(Geom_Curve) curve = please[i];
              First = curve->FirstParameter();
              Last = curve->LastParameter();
              curve->D0(First,PIn);
              curve->D0(Last,PFin);
              curve->D0((First+Last)/2.0,PMid);
              cout<<"To be added "<<Pnt(PIn)<<" | "<<Pnt(PMid)<<" | "<<Pnt(PFin)<<endl;
              Handle(Geom_Curve) temp = convert_keel_bspline.BSplineCurve();
              First = temp->FirstParameter();
              Last = temp->LastParameter();
              temp->D0(First,PIn);
              temp->D0(Last,PFin);
              temp->D0((First+Last)/2.0,PMid);
              cout<<"Keel "<<Pnt(PIn)<<" | "<<Pnt(PMid)<<" | "<<Pnt(PFin)<<endl;
              Handle(Geom_BoundedCurve) bcurve = Handle(Geom_BoundedCurve)::DownCast(curve);
              check = convert_keel_bspline.Add(bcurve,tolerance,0,1,0);
              one_failed = one_failed || (check == false);
              one_added = one_added || (check == true);
              added[i] = check;
              if (added[i] == true)
                {
                  added_count++;
                  cout<<"One Added!"<<endl;
                }
            }
      }

    AssertThrow(one_failed == false,
                ExcMessage("Bspline conversion has failed."));


    Handle(Geom_BSplineCurve) equiv_keel_bspline = convert_keel_bspline.BSplineCurve();

    Handle(Geom_Curve) bspline = convert_keel_bspline.BSplineCurve();

    GeomAdaptor_Curve AC = GeomAdaptor_Curve(bspline);
    Handle(GeomAdaptor_HCurve) ACH = new GeomAdaptor_HCurve(AC);
    gp_Ax3 Ax(gp_Pnt(0.0,0.0,0.0),gp_Dir(0.0,1.0,0.0),gp_Dir(1.0,0.0,0.0));
    ProjLib_ProjectOnPlane projOnPlane(Ax);

    projOnPlane.Load(ACH, 1e-7, Standard_True);

    Handle(Geom_BSplineCurve) equiv_keel_bspline_xy = projOnPlane.BSpline();

    Handle(Geom_Curve) bspline_xy = projOnPlane.BSpline();

    if (bspline_xy->IsCN(1))
      cout<<"XY edges curve is at least C1"<<endl;
    else
      cout<<"XY edges curve is not C1"<<endl;

    TopoDS_Edge new_edge = BRepBuilderAPI_MakeEdge(bspline_xy);
    new_edge.Location(L);


    return new_edge;
  }


  TopoDS_Shape extract_transom_edges(const TopoDS_Shape &in_shape,
                                     const unsigned int num_transom_edges,
                                     const double tolerance)
  {
    cout<<"Requested transom edges: "<<num_transom_edges<<endl;
    //// this is to loop on the edges
    TopExp_Explorer edgeExplorer(in_shape , TopAbs_EDGE);
    unsigned int numEdges = 0;

    TopoDS_Edge edge;
    TopLoc_Location L;
    Standard_Real First;
    Standard_Real Last;
    gp_Pnt PIn(0.0,0.0,0.0);
    gp_Pnt PFin(0.0,0.0,0.0);
    gp_Pnt PMid(0.0,0.0,0.0);
    TColGeom_Array1OfCurve curve_array(0,num_transom_edges);
    std::multiset< pair<unsigned int,double> , edge_ordering_rule > edges_set;
    typedef std::multiset< std::pair<unsigned int,double>, edge_ordering_rule >::iterator iterator;


    while (edgeExplorer.More())
      {
        edge = TopoDS::Edge(edgeExplorer.Current());
        if (BRep_Tool::Degenerated(edge))
          {
          }
        else
          {
            Handle(Geom_Curve) curve = BRep_Tool::Curve(edge,L,First,Last);
            curve->D0((First+Last)/2.0,PMid);
            //cout<<"** "<<numEdges<<" "<<Pnt(PMid)<<endl;
            edges_set.insert(std::pair<unsigned int, double >(numEdges, PMid.X()));
          }
        ++numEdges;
        edgeExplorer.Next();
      }

    //for (iterator pos=edges_set.begin(); pos!=edges_set.end();++pos)
    //    cout<<"*** "<<pos->first<<" "<<pos->second<<endl;


    iterator it = edges_set.end();
    for (unsigned int i=0; i<num_transom_edges; ++i)
      {
        it--;
        numEdges = 0;
        edgeExplorer.ReInit();
        while (edgeExplorer.More())
          {
            edge = TopoDS::Edge(edgeExplorer.Current());
            if (BRep_Tool::Degenerated(edge))
              {
              }
            else
              {
                Handle(Geom_Curve) curve = BRep_Tool::Curve(edge,L,First,Last);
                if (numEdges == it->first)
                  {
                    curve->D0((First+Last)/2.0,PMid);
                    //cout<<"* "<<i<<" "<<Pnt(PMid)<<endl;
                    curve_array(i) = curve;
                  }
              }
            ++numEdges;
            edgeExplorer.Next();
          }
      }


    Handle(Geom_Curve) curve1 = curve_array(0);
    Handle(Geom_BoundedCurve) bcurve1 = Handle(Geom_BoundedCurve)::DownCast(curve1);
    GeomConvert_CompCurveToBSplineCurve convert_transom_bspline(bcurve1,Convert_TgtThetaOver2);
    bool check = false, one_added = true, one_failed=true;
    vector<bool> added(num_transom_edges, false);
    added[0] = true;
    while (one_added == true)
      {
        one_added = false;
        one_failed = false;
        for (unsigned int i=1; i<num_transom_edges; ++i)
          {
            if (added[i] == false)
              {
                Handle(Geom_Curve) curve = curve_array(i);
                Handle(Geom_BoundedCurve) bcurve = Handle(Geom_BoundedCurve)::DownCast(curve);
                check = convert_transom_bspline.Add(bcurve,tolerance,0,1,0);
                one_failed = one_failed || (check == false);
                one_added = one_added || (check == true);
                added[i] = check;
              }
          }
      }

    AssertThrow(one_failed == false,
                ExcMessage("Bspline conversion of transom edges has failed."));


    Handle(Geom_BSplineCurve) equiv_transom_bspline = convert_transom_bspline.BSplineCurve();

    Handle(Geom_Curve) bspline = convert_transom_bspline.BSplineCurve();


    if (bspline->IsCN(1))
      cout<<"Transom curve is at least C1"<<endl;
    else
      cout<<"Transom curve is not C1"<<endl;


    TopoDS_Edge new_transom_edge = BRepBuilderAPI_MakeEdge(bspline);
    new_transom_edge.Location(L);

// These lines can be used to dump the keel bslpine on an .igs file
    /*
    IGESControl_Controller::Init();
    IGESControl_Writer ICW ("MM", 0);
    Standard_Boolean ok = ICW.AddShape (new_transom_edge);
    ICW.ComputeModel();
    Standard_Boolean OK = ICW.Write ("MyTransom.igs");
    */


    return new_transom_edge;



  }

  TopoDS_Shape merge_surfaces(const TopoDS_Shape &/*in_shape*/,
                              const double /*tolerance*/)
  {

    TopoDS_Shape shape;
    /*
       TopExp_Explorer faceExplorer(sh, TopAbs_FACE);
       TopoDS_Face face;
       unsigned int face_count = 0;
       while(faceExplorer.More())
         {
    face = TopoDS::Face(faceExplorer.Current());
           //Standard_Real umin, umax, vmin, vmax;
           //BRepTools::UVBounds(face, umin, umax, vmin, vmax);          // create surface
           //Handle(Geom_Surface) surf=BRep_Tool::Surface(face);          // get surface properties
    faceExplorer.Next();
           ++face_count;
         }

       TColGeom_Array2OfBezierSurface surface_array(0, face_count, 0, 0);
       faceExplorer.Reinit();
       face_count = 0;
       while(faceExplorer.More())
         {
    face = TopoDS::Face(faceExplorer.Current());
           //Standard_Real umin, umax, vmin, vmax;
           //BRepTools::UVBounds(face, umin, umax, vmin, vmax);
           Handle(Geom_Surface) surf=BRep_Tool::Surface(face);          // create surface
           surface_array.SetValue(face_count,0,surf);
    faceExplorer.Next();
           ++face_count;
         }

       AssertThrow(numKeelEdges>0,
     ExcMessage("No edges on xz plane"));

       Handle(Geom_Curve) curve1 = curve_array(0);
       Handle(Geom_BoundedCurve) bcurve1 = Handle(Geom_BoundedCurve)::DownCast(curve1);
       GeomConvert_CompCurveToBSplineCurve convert_keel_bspline(bcurve1,Convert_TgtThetaOver2);
       bool check = false, one_added = true, one_failed=true;
       vector<bool> added(numKeelEdges, false);
       added[0] = true;
       while(one_added == true)
         {
    one_added = false;
    one_failed = false;
    for (unsigned int i=1; i<numKeelEdges; ++i)
      if(added[i] == false)
        {
          Handle(Geom_Curve) curve = curve_array(i);
          Handle(Geom_BoundedCurve) bcurve = Handle(Geom_BoundedCurve)::DownCast(curve);
          check = convert_keel_bspline.Add(bcurve,1e-5,0,1,0);
          one_failed = one_failed || (check == false);
          one_added = one_added || (check == true);
          added[i] = check;
        }
         }

       AssertThrow(one_failed == false,
     ExcMessage("Bspline convertion has failed."));


       Handle(Geom_BSplineCurve) equiv_keel_bspline = convert_keel_bspline.BSplineCurve();

       Handle(Geom_Curve) bspline = convert_keel_bspline.BSplineCurve();
       return BRepBuilderAPI_MakeEdge(bspline);
    */
    return shape;

  }

  void intersect_plane(const TopoDS_Shape &in_shape,
                       TopoDS_Shape &out_shape,
                       const double c_x,
                       const double c_y,
                       const double c_z,
                       const double c,
                       const double tolerance)
  {
    Handle(Geom_Plane) plane = new Geom_Plane(c_x,c_y,c_z,c);

    // This is to loop on the faces,
    // that extracts all
    // intersections.

    //TopExp_Explorer faceExplorer(in_shape , TopAbs_FACE);

    std::vector<Handle_Geom_BoundedCurve> intersections;

    BRepAlgo_Section section(in_shape, plane);
    TopoDS_Shape edges = section.Shape();

    TopoDS_Edge edge;
    TopLoc_Location L;
    Standard_Real First;
    Standard_Real Last;
    gp_Pnt PIn(0.0,0.0,0.0);
    gp_Pnt PFin(0.0,0.0,0.0);
    gp_Pnt PMid(0.0,0.0,0.0);
    TopExp_Explorer edgeExplorer(edges , TopAbs_EDGE);
    unsigned int count = 0;
    while (edgeExplorer.More())
      {
        edge = TopoDS::Edge(edgeExplorer.Current());
        Handle(Geom_Curve) curve = BRep_Tool::Curve(edge,L,First,Last);
        intersections.push_back(Handle(Geom_BoundedCurve)::DownCast(curve));

        curve->D0(First,PIn);
        curve->D0(Last,PFin);
        curve->D0((First+Last)/2.0,PMid);
        //cout<<Pnt(PIn)<<" | "<<Pnt(PMid)<<" | "<<Pnt(PFin)<<" ---> ";
        //cout<<fabs(double(PIn.Y())+double(PFin.Y())+double(PMid.Y()))/3.0<<endl;

        cout<<count<<":  "<<Pnt(PIn)<<" | "<<Pnt(PMid)<<" | "<<Pnt(PFin)<<endl;

        //IGESControl_Controller::Init();
        //IGESControl_Writer ICW ("MM", 0);
        //Standard_Boolean ok = ICW.AddShape (edge);
        //ICW.ComputeModel();
        //Standard_Boolean OK = ICW.Write ("water_line.igs");
        //count++;

        edgeExplorer.Next();

      }

    // Now we build a single bspline
    // out of all the geometrical
    // curves
    unsigned int numIntersEdges = intersections.size();
    for (unsigned int i = 0; i<intersections.size(); ++i)
      {
        if (intersections[i]->Value(intersections[i]->FirstParameter()).X() >
            intersections[i]->Value(intersections[i]->LastParameter()).X()  )
          intersections[i]->Reverse();
      }

    GeomConvert_CompCurveToBSplineCurve
    convert_bspline(intersections[0], Convert_TgtThetaOver2);
    bool check = false, one_added = true, one_failed=true;
    vector<bool> added(numIntersEdges, false);
    added[0] = true;
    while (one_added == true)
      {
        one_added = false;
        one_failed = false;
        for (unsigned int i=1; i<numIntersEdges; ++i)
          if (added[i] == false)
            {
              Handle(Geom_Curve) curve = intersections[i];
              Handle(Geom_BoundedCurve) bcurve = Handle(Geom_BoundedCurve)::DownCast(curve);
              check = convert_bspline.Add(bcurve,tolerance,0,1,0);
              one_failed = one_failed || (check == false);
              one_added = one_added || (check == true);
              added[i] = check;
              //cout<<i<<" -->  "<<added[i]<<"  "<<false<<endl;
            }
      }

    AssertThrow(one_failed == false,
                ExcMessage("Bspline convertion of intersection with plane has failed."));

    Handle(Geom_Curve) bspline = convert_bspline.BSplineCurve();

    out_shape = BRepBuilderAPI_MakeEdge(bspline);

    if (bspline->IsCN(1))
      cout<<"Intersection with plane is at least a C1 curve"<<endl;
    else
      cout<<"Intersection with plane is not a C1 curve "<<endl;


  }



  Handle_Geom_Curve  interpolation_curve_points_sort(std::vector<Point<3> > &curve_points,
                                                     Point<3> direction)
  {

    unsigned int n_vertices = curve_points.size();

    if (direction*direction > 0)
      {
        std::sort(curve_points.begin(), curve_points.end(),
                  boost::bind(&point_compare, _1, _2, direction));
      }

    // set up array of vertices
    Handle(TColgp_HArray1OfPnt) vertices = new TColgp_HArray1OfPnt(1,n_vertices);
    for (unsigned int vertex=0; vertex<n_vertices; ++vertex)
      {
        vertices->SetValue(vertex+1,OpenCascade::Pnt(curve_points[vertex]));
      }


    GeomAPI_Interpolate bspline_generator(vertices, Standard_False, 1.0e-7);
    bspline_generator.Perform();
    AssertThrow( (bspline_generator.IsDone()), ExcMessage("Interpolated bspline generation failed"));
    //Handle(Geom_BSplineCurve) bspline = new Geom_BSplineCurve(*(bspline_generator.Curve().operator->()));
    Handle(Geom_BSplineCurve) bspline = bspline_generator.Curve();

    return bspline;
  }


  Handle_Geom_Curve  interpolation_curve(std::vector<Point<3> > &curve_points)
  {

    unsigned int n_vertices = curve_points.size();

    // set up array of vertices
    Handle(TColgp_HArray1OfPnt) vertices = new TColgp_HArray1OfPnt(1,n_vertices);
    for (unsigned int vertex=0; vertex<n_vertices; ++vertex)
      {
        vertices->SetValue(vertex+1,OpenCascade::Pnt(curve_points[vertex]));
      }


    GeomAPI_Interpolate bspline_generator(vertices, Standard_False, 1.0e-7);
    bspline_generator.Perform();

    //Handle(Geom_Curve) bspline = new Geom_BSplineCurve(*(bspline_generator.Curve().operator->()));

    Handle(Geom_BSplineCurve) bspline = bspline_generator.Curve();
    return bspline;
  }



  void feature_edges_detection()
  {
    //TopoDS_Shape sh = read_IGES("/home/amola/workspace/FRANCO_TETGEN/tower_bridge_cf_IGES.igs",0.001);
    //TopoDS_Shape sh = read_IGES("/home/amola/workspace/FRANCO_TETGEN/VIOLA_BRIDGE/viola-bridge-1/viola_bridge.iges",0.0254);
    TopoDS_Shape sh = read_IGES("/home/amola/workspace/openship/trunk/WaveBEM/utilities/kcs.iges",0.001);
    //TopoDS_Shape sh = read_IGES("/home/amola/workspace/openship/trunk/WaveBEM/utilities/goteborg.iges",0.001);
    //TopoDS_Shape sh = read_IGES("/home/amola/workspace/FRANCO_TETGEN/cognit.iges",0.001);
    //TopoDS_Shape sh = read_IGES("/home/amola/workspace/FRANCO_TETGEN/COFANO/cofanoSpostatoTanto.igs",0.001);
    //TopoDS_Shape sh = read_IGES("/home/amola/workspace/FRANCO_TETGEN/CANOTTAGGIO/doppio_placido.iges",0.001);
    //TopoDS_Shape sh = read_IGES("/home/amola/workspace/openship/trunk/WaveBEM/utilities/fc_ship.iges",0.001);
    //TopoDS_Shape sh = read_IGES("/home/amola/workspace/FRANCO_TETGEN/ESEMPIO_HANG/22imr_retro.igs",0.001);

    using dealii::numbers::PI;
    Standard_Real AngTol = 10.0/180.0;
    Standard_Real LinTol = 2e-3;

    std::vector<TopoDS_Edge> feature_edges;

    TopExp_Explorer faceExplorer(sh, TopAbs_FACE);
    TopoDS_Face face;
    unsigned int face_count = 0;
    std::vector<Bnd_Box> facesBndBoxes;

    while (faceExplorer.More())
      {
        face = TopoDS::Face(faceExplorer.Current());
        //BRepMesh::Mesh(face, 0.01f);
        Bnd_Box bndBox;
        BRepBndLib::Add(face, bndBox);
        bndBox.Enlarge(100*LinTol);
        facesBndBoxes.push_back(bndBox);
        ++face_count;
        faceExplorer.Next();
      }

    unsigned int num_faces = face_count;
    std::vector< std::set<unsigned int> > faces_neighbors_sets(num_faces);

    face_count = 0;
    faceExplorer.Init(sh, TopAbs_FACE);

    while (faceExplorer.More())
      {
        face = TopoDS::Face(faceExplorer.Current());
        Handle(Geom_Surface) ssurface = BRep_Tool::Surface(face);
        double u1,u2,v1,v2;
        ssurface->Bounds(u1,u2,v1,v2);
        cout<<"Face "<<face_count<<"   Center ("<<Pnt(ssurface->Value((u1+u2)/2.0,(v1+v2)/2.0))<<")"<<endl;
        faces_neighbors_sets[face_count].insert(face_count);
        TopExp_Explorer edgeExplorer(face, TopAbs_EDGE);
        TopoDS_Edge edge;
        gp_Pnt PMid;
        gp_Vec TangMid;
        unsigned int edge_count = 0;
        while (edgeExplorer.More())
          {
            edge = TopoDS::Edge(edgeExplorer.Current());
            if (BRep_Tool::Degenerated(edge))
              {
              }
            else
              {
                TopLoc_Location L;
                Standard_Real First;
                Standard_Real Last;
                Handle(Geom_Curve) curve = BRep_Tool::Curve(edge,L,First,Last);
                curve->D1((First+Last)/2.0,PMid,TangMid);
                TangMid.Normalize();
                Handle(Geom_Surface) this_surface = BRep_Tool::Surface(face);
                ShapeAnalysis_Surface self_projector(this_surface);
                gp_Pnt2d self_proj_params = self_projector.ValueOfUV(PMid, 1e-7);
                GeomLProp_SLProps self_props(this_surface, self_proj_params.X(), self_proj_params.Y(),1, 1e-7);
                gp_Dir self_Normal = self_props.Normal();
                //cout<<"Face "<<face_count<<"  Point "<<edge_count<<": "<<Pnt(PMid)<<endl;

                TopExp_Explorer innerFaceExplorer(sh, TopAbs_FACE);
                TopoDS_Face innerFace;
                unsigned int inner_face_count = 0;
                double minDistance = 1e7;
                gp_Pnt nearest_projection;
                unsigned int projFaceID;
                gp_Dir Normal;
                while (innerFaceExplorer.More())
                  {
                    if ((inner_face_count != face_count) && (!facesBndBoxes[inner_face_count].IsOut(PMid)) )
                      {
                        innerFace = TopoDS::Face(innerFaceExplorer.Current());
                        BRepExtrema_DistShapeShape distSS(innerFace,BRepBuilderAPI_MakeVertex(PMid));
                        distSS.Perform();
                        if (distSS.IsDone())
                          {
                            //cout<<"Distance Eval: "<<distSS.Value()<<" Point: "<<Pnt(distSS.PointOnShape1(1))<<endl;
                          }
                        else
                          {
                            cout<<"WARNING: Distance Eval Computation FAILED: should stop everything here "<<endl;
                          }
                        double U1,U2,V1,V2;
                        Handle(Geom_Surface) SurfToProj = BRep_Tool::Surface(innerFace);
                        SurfToProj->Bounds(U1,U2,V1,V2);
                        //cout<<"Inner Face "<<inner_face_count<<"   Center ("<<Pnt(SurfToProj->Value((U1+U2)/2.0,(V1+V2)/2.0))<<")"<<endl;
                        ShapeAnalysis_Surface projector(SurfToProj);
                        gp_Pnt2d proj_params = projector.ValueOfUV(distSS.PointOnShape1(1), 1e-7);
                        gp_Pnt projection;
                        SurfToProj->D0(proj_params.X(),proj_params.Y(),projection);
                        if (distSS.Value() < minDistance)
                          {
                            minDistance = distSS.Value();
                            nearest_projection = distSS.PointOnShape1(1);
                            projFaceID=inner_face_count;
                            GeomLProp_SLProps props(SurfToProj, proj_params.X(), proj_params.Y(),1, 1e-7);
                            Normal = props.Normal();
                          }

                      }
                    innerFaceExplorer.Next();
                    ++inner_face_count;
                  }

                //cout<<"Nearest Projection "<<Pnt(nearest_projection)<<" on face "<<projFaceID<<endl;
                if (minDistance < LinTol)
                  {
                    //cout<<"Fake edge? "<<endl;
                    //cout<<"n_a("<<self_Normal.X()<<","<<self_Normal.Y()<<","<<self_Normal.Z()<<")"<<endl;
                    //cout<<"n_b("<<Normal.X()<<","<<Normal.Y()<<","<<Normal.Z()<<")"<<endl;
                    //cout<<"Angle "<<180/PI*self_Normal.Angle(Normal)<<"  or   ";
                    //cout<<180.0/PI*self_Normal.Angle(-Normal)<<endl;

                    if ( fmin(fabs(180/PI*self_Normal.Angle(Normal)),
                              fabs(180/PI*self_Normal.Angle(-Normal))) < 5.0)
                      {
                        //cout<<"From closest point and normal could be the edge is fake "<<face_count<<"----->"<<projFaceID<<endl;

                        // we try here to understand if the normals are "aligned" because the faces have smooth
                        // junction or instead they form a cusp
                        innerFaceExplorer.ReInit();
                        for (unsigned int kk=0; kk<projFaceID; ++kk)
                          innerFaceExplorer.Next();
                        innerFace = TopoDS::Face(innerFaceExplorer.Current());
                        Handle(Geom_Surface) SurfToProj = BRep_Tool::Surface(innerFace);
                        gp_Vec normal_vec(Normal);
                        gp_Vec cross_vec = TangMid.Crossed(normal_vec);
                        gp_Pnt P_1(PMid.X()+50*LinTol*cross_vec.X(),PMid.Y()+50*LinTol*cross_vec.Y(),PMid.Z()+50*LinTol*cross_vec.Z());
                        gp_Pnt P_2(PMid.X()-50*LinTol*cross_vec.X(),PMid.Y()-50*LinTol*cross_vec.Y(),PMid.Z()-50*LinTol*cross_vec.Z());
                        //cout<<"Point: "<<Pnt(PMid)<<endl;
                        //cout<<"Normal: "<<normal_vec.X()<<" "<<normal_vec.Y()<<" "<<normal_vec.Z()<<endl;
                        //cout<<"Edge Tangent: "<<TangMid.X()<<" "<<TangMid.Y()<<" "<<TangMid.Z()<<endl;
                        //cout<<"P1: "<<Pnt(P_1)<<"  P2: "<<Pnt(P_2)<<endl;
                        BRepExtrema_DistShapeShape distSS1(innerFace,BRepBuilderAPI_MakeVertex(P_1));
                        distSS1.Perform();
                        double distance1 = distSS1.Value();
                        gp_Pnt proj_point1 = distSS1.PointOnShape1(1);
                        distSS1.LoadS2(BRepBuilderAPI_MakeVertex(P_2));
                        distSS1.Perform();
                        if (distSS1.Value() < distance1)
                          proj_point1 = distSS1.PointOnShape1(1);
                        //cout<<"F1 P1 Dist: "<<distance1<<"  F1 P2 Dist: "<<distSS1.Value()<<endl;
                        BRepExtrema_DistShapeShape distSS2(face,BRepBuilderAPI_MakeVertex(P_1));
                        distSS2.Perform();
                        double distance2 = distSS2.Value();
                        gp_Pnt proj_point2 = distSS2.PointOnShape1(1);
                        distSS2.LoadS2(BRepBuilderAPI_MakeVertex(P_2));
                        distSS2.Perform();
                        if (distSS2.Value() < distance2)
                          proj_point2 = distSS2.PointOnShape1(1);
                        //cout<<"F2 P1 Dist: "<<distance2<<"  F2 P2 Dist: "<<distSS2.Value()<<endl;

                        gp_Vec v1(PMid,proj_point1);
                        gp_Vec v2(PMid,proj_point2);
                        //cout<<"V1: "<<v1.X()<<" "<<v1.Y()<<" "<<v1.Z()<<endl;
                        //cout<<"V2: "<<v2.X()<<" "<<v2.Y()<<" "<<v2.Z()<<endl;
                        if ( fabs(180/PI*v1.Angle(v2)) > 90)
                          {
                            //cout<<"Yes, the edge is fake "<<face_count<<"----->"<<projFaceID<<endl;
                            faces_neighbors_sets[face_count].insert(projFaceID);
                          }
                        else
                          {
                            //cout<<"No, the edge is sharp "<<endl;
                            feature_edges.push_back(edge);
                          }
                      }
                    else
                      {
                        //cout<<"No, the edge is sharp "<<endl;
                        feature_edges.push_back(edge);
                      }
                  }
                else
                  {
                    //cout<<"Outer edge "<<endl;
                    feature_edges.push_back(edge);
                  }

              }
            edgeExplorer.Next();
            ++edge_count;
          }
        faceExplorer.Next();
        ++face_count;
      }

    cout<<endl;
    cout<<"Let's group faces into colors "<<endl;
    cout<<endl;
    for (unsigned int i=0; i<num_faces; ++i)
      {
        cout<<"Face "<<i<<": ";
        for (std::set<unsigned int>::iterator pos=faces_neighbors_sets[i].begin(); pos!=faces_neighbors_sets[i].end(); ++pos)
          cout<<" "<<*pos;
        cout<<endl;
      }
    std::set < std::set<unsigned int> > colors;
    for (unsigned int k=0; k<num_faces; ++k)
      {
        std::set<unsigned int> color=faces_neighbors_sets[k];
        unsigned int prev_size = color.size();
        unsigned int new_size = color.size()+1;
        while (prev_size != new_size)
          {
            prev_size=color.size();
            for (unsigned int i=0; i<num_faces; ++i)
              {
                for (std::set<unsigned int>::iterator pos=color.begin(); pos!=color.end(); ++pos)
                  if (faces_neighbors_sets[i].count(*pos) > 0)
                    {
                      for (std::set<unsigned int>::iterator pos=faces_neighbors_sets[i].begin(); pos!=faces_neighbors_sets[i].end(); ++pos)
                        color.insert(*pos);
                      break;
                    }
              }
            new_size=color.size();
          }
        colors.insert(color);
      }


    /*
    IGESControl_Controller::Init();
    IGESControl_Writer ICW ("MM", 0);
    for (unsigned int i=0; i<feature_edges.size(); ++i)
        Standard_Boolean ok = ICW.AddShape (feature_edges[i]);
    ICW.ComputeModel();
    Standard_Boolean OK = ICW.Write ("feature_edges.igs");
    */

    // trying to organize the edges in feature_edges vector into a limited number of
    // curves
    std::vector<bool> active_flag(feature_edges.size(),true);
    std::vector<gp_Pnt> start_points(feature_edges.size());
    std::vector<gp_Pnt> mid_points(feature_edges.size());
    std::vector<gp_Pnt> end_points(feature_edges.size());
    std::vector<gp_Dir> start_tangents(feature_edges.size());
    std::vector<gp_Dir> end_tangents(feature_edges.size());
    for (unsigned int i=0; i<feature_edges.size(); ++i)
      {
        TopLoc_Location L;
        Standard_Real First;
        Standard_Real Last;
        Handle(Geom_Curve) curve = BRep_Tool::Curve(feature_edges[i],L,First,Last);
        gp_Pnt PMid;
        curve->D0(First,start_points[i]);
        curve->D0(Last,end_points[i]);
        GeomLProp_CLProps curve_props_First(curve, First,1, LinTol);
        GeomLProp_CLProps curve_props_Last(curve, Last, 1, LinTol);
        if (curve_props_First.IsTangentDefined() == false )
          std::cout<<"Warning! Tangent vector at start of edge "<<i<<" is undefined!"<<std::endl;
        if (curve_props_Last.IsTangentDefined() == false )
          std::cout<<"Warning! Tangent vector at end of edge "<<i<<" is undefined!"<<std::endl;
        curve_props_First.Tangent(start_tangents[i]);
        curve_props_Last.Tangent(end_tangents[i]);
        GeomAdaptor_Curve AC(curve);
        Standard_Real arcLength = GCPnts_AbscissaPoint::Length(AC,First,Last);
        if (arcLength < LinTol)
          active_flag[i] = false;
        //cout<<"arcLength "<<arcLength<<endl;
        GCPnts_AbscissaPoint AP(AC, arcLength/2.0, First);
        //cout<<"direction*arcLength*distance "<<direction*arcLength*distance<<endl;
        Standard_Real t2 = AP.Parameter();
        AssertThrow((First <= t2) &&
                    (t2 <= Last),
                    ExcMessage("Parameter 3 is out of range!"));
        mid_points[i] = AC.Value(t2);

        /*
        std::cout<<std::endl;
        std::cout<<i<<std::endl;
        std::cout<<"Start Point: "<<Pnt(start_points[i])<<std::endl;
        std::cout<<"Mid Point: "<<Pnt(mid_points[i])<<std::endl;
        std::cout<<"End Point: "<<Pnt(end_points[i])<<std::endl;
        std::cout<<"Start Tangent: "<<start_tangents[i].X()<<" "<<start_tangents[i].Y()<<" "<<start_tangents[i].Z()<<std::endl;
        std::cout<<"End Tangent: "<<end_tangents[i].X()<<" "<<end_tangents[i].Y()<<" "<<end_tangents[i].Z()<<std::endl;
        std::cout<<"F: "<<First<<" "<<"  L: "<<Last<<"  Or: "<<feature_edges[i].Orientation()<<std::endl;
        std::cout<<"F?: "<<curve_props_First.IsTangentDefined()<<" "<<"  L?: "<<curve_props_Last.IsTangentDefined()<<std::endl;
        */
      }

    // we first take out all possible copies of same edge: here detected as
    // edges with same end, mid and last point, the procedude could be changed if needed
    for (unsigned int i=0; i<feature_edges.size(); ++i)
      if (active_flag[i])
        {
          for (unsigned int j=i+1; j<feature_edges.size(); ++j)
            if (mid_points[i].IsEqual(mid_points[j],LinTol))
              if (start_points[i].IsEqual(start_points[j],LinTol) || start_points[i].IsEqual(end_points[j],LinTol))
                if (end_points[i].IsEqual(end_points[j],LinTol) || end_points[i].IsEqual(start_points[j],LinTol))
                  active_flag[j] = false;
        }

    std::vector<TopoDS_Edge> unified_edge_vector;
    std::vector<Handle(Geom_Curve)> unified_curve_handles_vector;
    unsigned int unified_edges_count=0;

    for (unsigned int k=0; k<feature_edges.size(); ++k)
      if (active_flag[k])
        {
          //std::cout<<"Unified edge "<<unified_edges_count<<" is to be grown from edge: "<<k<<std::endl;
          std::set<unsigned int> joined_edges;
          joined_edges.insert(k);
          active_flag[k] = false;
          gp_Pnt curve_start_point = start_points[k];
          gp_Pnt curve_end_point = end_points[k];
          gp_Vec curve_start_tangent = start_tangents[k];
          gp_Vec curve_end_tangent = end_tangents[k];
          unsigned int added_count = 1;

          while (added_count > 0)
            {
              added_count = 0;
              for (unsigned int i=0; i<feature_edges.size(); ++i)
                if (active_flag[i])
                  {
                    if (curve_end_point.IsEqual(start_points[i],LinTol) && (curve_end_tangent.Angle(start_tangents[i]) < AngTol))
                      {
                        //cout<<"aFOUND!  "<<i<<endl;
                        joined_edges.insert(i);
                        active_flag[i] = false;
                        curve_end_point = end_points[i];
                        curve_end_tangent = end_tangents[i];
                        added_count++;
                      }
                    else if (curve_end_point.IsEqual(end_points[i],LinTol) && (PI-curve_end_tangent.Angle(end_tangents[i]) < AngTol))
                      {
                        //cout<<"bFOUND!  "<<i<<endl;
                        joined_edges.insert(i);
                        active_flag[i] = false;
                        curve_end_point = start_points[i];
                        curve_end_tangent = -start_tangents[i];
                        added_count++;
                      }
                    else if (curve_start_point.IsEqual(start_points[i],LinTol) && (PI-curve_start_tangent.Angle(start_tangents[i]) < AngTol))
                      {
                        //cout<<"cFOUND!  "<<i<<endl;
                        joined_edges.insert(i);
                        active_flag[i] = false;
                        curve_start_point = end_points[i];
                        curve_start_tangent = -end_tangents[i];
                        added_count++;
                      }
                    else if (curve_start_point.IsEqual(end_points[i],LinTol) && (curve_start_tangent.Angle(end_tangents[i]) < AngTol))
                      {
                        //cout<<"dFOUND!  "<<i<<endl;
                        joined_edges.insert(i);
                        active_flag[i] = false;
                        curve_start_point = start_points[i];
                        curve_start_tangent = start_tangents[i];
                        added_count++;
                      }
                  }
              //std::cout<<"Added: "<<added_count<<endl;
              if (curve_start_point.IsEqual(curve_end_point,LinTol))
                added_count = 0;
            }// end while

          std::vector<TopoDS_Edge> edges_to_be_joined;
          std::cout<<"Unified edge "<<unified_edges_count<<" is grown from edge: "<<k<<std::endl;
          for (std::set<unsigned int>::iterator pos=joined_edges.begin(); pos!=joined_edges.end(); ++pos)
            {
              std::cout<<*pos<<" ";
              edges_to_be_joined.push_back(feature_edges[*pos]);
            }
          std::cout<<std::endl;

          // after we found all the "joinable" edges, we join them
          TopLoc_Location L;
          Standard_Real First;
          Standard_Real Last;
          unsigned int numKeelEdges = edges_to_be_joined.size();
          Handle(Geom_Curve) curve1 = BRep_Tool::Curve(edges_to_be_joined[0],L,First,Last);
          //Handle(Geom_BoundedCurve) bcurve1 = Handle(Geom_BoundedCurve)::DownCast(curve1);
          Handle(Geom_TrimmedCurve) bcurve1 = new Geom_TrimmedCurve(curve1, First, Last);

          GeomConvert_CompCurveToBSplineCurve convert_keel_bspline(bcurve1,Convert_TgtThetaOver2);

          bool check = false, one_added = true, one_failed=true;
          vector<bool> added(numKeelEdges, false);
          added[0] = true;
          while (one_added == true)
            {
              one_added = false;
              one_failed = false;
              for (unsigned int ii=1; ii<numKeelEdges; ++ii)
                if (added[ii] == false)
                  {
                    //cout<<"Curve "<<ii<<" of "<<numKeelEdges<<endl;
                    Handle(Geom_Curve) curve = BRep_Tool::Curve(edges_to_be_joined[ii],L,First,Last);
                    Handle(Geom_TrimmedCurve) bcurve = new Geom_TrimmedCurve(curve, First, Last);
                    gp_Pnt PIn, PMid, PFin;
                    //First = curve->FirstParameter();
                    //Last = curve->LastParameter();
                    curve->D0(First,PIn);
                    curve->D0(Last,PFin);
                    curve->D0((First+Last)/2.0,PMid);
                    //cout<<"To be added "<<Pnt(PIn)<<" | "<<Pnt(PMid)<<" | "<<Pnt(PFin)<<endl;
                    Handle(Geom_Curve) temp = convert_keel_bspline.BSplineCurve();
                    First = temp->FirstParameter();
                    Last = temp->LastParameter();
                    temp->D0(First,PIn);
                    temp->D0(Last,PFin);
                    temp->D0((First+Last)/2.0,PMid);
                    //cout<<"Keel "<<Pnt(PIn)<<" | "<<Pnt(PMid)<<" | "<<Pnt(PFin)<<endl;
                    //Handle(Geom_BoundedCurve) bcurve = Handle(Geom_BoundedCurve)::DownCast(curve);
                    check = convert_keel_bspline.Add(bcurve,LinTol,0,1,0);
                    one_failed = one_failed || (check == false);
                    one_added = one_added || (check == true);
                    added[ii] = check;
                  }
            }

          AssertThrow(one_failed == false,
                      ExcMessage("Bspline conversion has failed."));

          Handle(Geom_BSplineCurve) equiv_keel_bspline = convert_keel_bspline.BSplineCurve();

          Handle(Geom_Curve) bspline = convert_keel_bspline.BSplineCurve();
          cout<<Pnt(bspline->Value(bspline->FirstParameter()))<<" and "<<Pnt(bspline->Value(bspline->LastParameter()))<<endl;
          unified_curve_handles_vector.push_back(bspline);


          //if (bspline->IsCN(1))
          //   cout<<"Transom curve is at least C1"<<endl;
          //else
          //   cout<<"Transom curve is not C1"<<endl;

          TopoDS_Edge unified_edge = BRepBuilderAPI_MakeEdge(bspline);
          unified_edge_vector.push_back(unified_edge);

          IGESControl_Controller::Init();
          IGESControl_Writer ICW ("MM", 0);
          Standard_Boolean ok = ICW.AddShape(unified_edge);
          ICW.ComputeModel();
          std::string filename = ( "unified_edges_" +
                                   Utilities::int_to_string(unified_edges_count+1) +
                                   ".igs" );
          std::ofstream file(filename.c_str());
          Standard_Boolean OK = ICW.Write (file);
          unified_edges_count++;
        } // end for



    cout<<"Number of unified edges generated: "<<unified_edge_vector.size()<<endl;

    BRep_Builder builder;

    std::vector<TopoDS_Compound> color_compounds;

    unsigned int color_count=0;
    cout<<"Color "<<color_count<<"?"<<endl;
    for (std::set< std::set<unsigned int> >::iterator posExt=colors.begin(); posExt!=colors.end(); ++posExt)
      {
        cout<<"COLOR "<<color_count<<":";
        for (std::set<unsigned int>::iterator pos= (*posExt).begin(); pos!= (*posExt).end(); ++pos)
          {
            cout<<" "<<*pos;
          }
        cout<<endl;

        faceExplorer.Init(sh, TopAbs_FACE);
        TopoDS_Compound color_compound;
        builder.MakeCompound(color_compound);

        unsigned int face_count = 0;
        cout<<"Added Faces :";
        while (faceExplorer.More())
          {
            face = TopoDS::Face(faceExplorer.Current());
            if (posExt->count(face_count) == 1)
              {
                builder.Add(color_compound,face);
                cout<<" "<<face_count;
              }
            faceExplorer.Next();
            ++face_count;
          }
        cout<<endl;
        color_compounds.push_back(color_compound);

        IGESControl_Controller::Init();
        IGESControl_Writer ICW ("MM", 0);
        Standard_Boolean ok = ICW.AddShape(color_compounds[color_count]);
        ICW.ComputeModel();
        std::string filename = ( "color_compound_" +
                                 Utilities::int_to_string(color_count+1) +
                                 ".igs" );
        std::ofstream file(filename.c_str());
        Standard_Boolean OK = ICW.Write (file);

        ++color_count;


        // let's try and mesh them here
        Standard_Integer deg = 3;
        Standard_Integer n_points = 30;
        Standard_Integer n_iter = 6;
        Standard_Real tol_3d = 1e-4;
        Standard_Real tol_2d = 1e-5;
        Standard_Real tol_ang = 0.001;
        Standard_Real tol_curv = 0.001;
        GeomPlate_BuildPlateSurface plate_surface_builder(deg, n_points, n_iter, tol_2d, tol_3d, tol_ang, tol_curv, Standard_False);
        plate_surface_builder.SetNbBounds(unified_edge_vector.size());
        for (unsigned int i=0; i<unified_edge_vector.size(); ++i)
          {
            TopoDS_Edge edge = unified_edge_vector[i];
            TopLoc_Location L;
            Standard_Real First;
            Standard_Real Last;
            Handle(Geom_Curve) curve = unified_curve_handles_vector[i]; //BRep_Tool::Curve(feature_edges[i],L,First,Last);
            cout<<Pnt(curve->Value(curve->FirstParameter()))<<" and "<<Pnt(curve->Value(curve->LastParameter()))<<endl;
            //cout<<Pnt(curve->Value(First))<<" and "<<Pnt(curve->Value(Last))<<endl;
            GeomAdaptor_HCurve geom_adaptor_hcurve(curve);
            const Handle(GeomAdaptor_HCurve)& aHC = new GeomAdaptor_HCurve(geom_adaptor_hcurve);
            Handle (GeomPlate_CurveConstraint) aConst = new GeomPlate_CurveConstraint (aHC, 0, n_points, tol_3d);
            plate_surface_builder.Add(aConst);
            gp_Pnt first_point;
            gp_Pnt last_point;
            aConst->D0(aConst->FirstParameter(), first_point);
            aConst->D0(aConst->LastParameter(), last_point);
            cout<<"## "<<Pnt(first_point)<<" and "<<Pnt(last_point)<<endl;
          }
        cout<<"Here"<<endl;
        plate_surface_builder.Perform();


        if (plate_surface_builder.IsDone())
          {
            cout<<"Made it! Error: "<<plate_surface_builder.G0Error()<<endl;
            Handle (GeomPlate_Surface) plate_surf = plate_surface_builder.Surface();
            Standard_Real dmax = plate_surface_builder.G0Error();
            TColgp_SequenceOfXY aS2d;
            TColgp_SequenceOfXYZ aS3d;
            plate_surface_builder.Disc2dContour (40, aS2d);
            plate_surface_builder.Disc3dContour (40, 0, aS3d);
            for (unsigned i=0; i<40; ++i)
              {
                gp_XY uv_point = aS2d.Value(i+1);
                gp_XYZ xyz_point = aS3d.Value(i+1);
                cout<<i<<endl;
                cout<<"("<<uv_point.X()<<","<<uv_point.Y()<<")"<<endl;
                cout<<"("<<xyz_point.X()<<","<<xyz_point.Y()<<","<<xyz_point.Z()<<")"<<endl;
              }
            //Standard_Real aMax = Max (aTol3d, 10. * aDMax);
            GeomPlate_PlateG0Criterion criterion (aS2d, aS3d, tol_3d);
            Standard_Integer max_bezier_pieces = 50;
            GeomPlate_MakeApprox make_approx(plate_surf, criterion, tol_3d, max_bezier_pieces, deg, GeomAbs_C0, 1.1);
            //GeomPlate_MakeApprox make_approx(plate_surf, tol_3d, max_bezier_pieces, deg, dmax, 0, GeomAbs_C0, 3.0);
            Handle(Geom_BSplineSurface) bspline_surf = make_approx.Surface();
            IGESControl_Controller::Init();
            IGESControl_Writer ICW2 ("MM", 0);
            Standard_Boolean ok = ICW2.AddGeom(bspline_surf);
            ICW2.ComputeModel();
            std::string new_filename = ( "approx_surface" +
                                         Utilities::int_to_string(1) +
                                         ".igs" );
            std::ofstream file(new_filename.c_str());
            Standard_Boolean OK = ICW2.Write (file);


            cout<<"Re-Made it!"<<endl;


            BRepBuilderAPI_MakeWire wireMaker;

            for (unsigned int i = 0; i< unified_curve_handles_vector.size(); ++i)
              {
                Handle(Geom_Curve) curve = unified_curve_handles_vector[i];
                TopoDS_Edge myEdge = BRepBuilderAPI_MakeEdge(curve);
                wireMaker.Add(myEdge);
              }

            TopoDS_Wire myWire = wireMaker.Wire();
            BRepBuilderAPI_MakeFace faceMaker(bspline_surf, myWire, Standard_True);
            TopoDS_Face trimmed_approx_face = faceMaker.Face();

            IGESControl_Controller::Init();
            IGESControl_Writer ICW3 ("MM", 0);
            Standard_Boolean ook = ICW3.AddShape(trimmed_approx_face);
            ICW3.ComputeModel();
            std::string neww_filename = ( "trimmed_approx_surface" +
                                          Utilities::int_to_string(1) +
                                          ".igs" );
            std::ofstream file2(neww_filename.c_str());
            Standard_Boolean OOK = ICW3.Write (file2);
            BRepMesh_IncrementalMesh (trimmed_approx_face, 0.7f);

            //BRepMesh::Mesh(trimmed_approx_face, 0.7f);

            BRepBuilderAPI_MakeFace faceMaker2(bspline_surf, 1e-6);
            TopoDS_Face aFace = faceMaker2.Face();
            TopLoc_Location loc;
            loc.Identity();
            Handle(Poly_Triangulation) tri = BRep_Tool::Triangulation(aFace, loc);




            if (tri.IsNull())
              {
                cout<<"Tringulation didn't work"<<endl;
                Bnd_Box bndBox;
                BRepBndLib::Add(trimmed_approx_face, bndBox);
                cout<<"Extent: "<<sqrt(bndBox.SquareExtent())<<endl;
                //bndBox.Enlarge(100*LinTol);
                //BRepMesh_FastDiscret mesh_fast_discret(0.0000001,trimmed_approx_face,bndBox, 0.1, Standard_True, Standard_True, Standard_True, Standard_True);
                //cout<<"Vertices: "<<mesh_fast_discret.NbVertices()<<endl;
                BRepMesh_IncrementalMesh aMesh(aFace, 0.25);
                aMesh.Perform();
                cout<<"** "<<aMesh.IsDone()<<endl;
                cout<<aMesh.GetStatusFlags()<<endl;

                //BRepMesh::Mesh(trimmed_approx_face, 0.25);
              }
            tri = BRep_Tool::Triangulation(aFace, loc);
            if  (tri.IsNull())
              {
                cout<<"Tringulation didn't work AGAIN"<<endl;
              }
            const TColgp_Array1OfPnt &nodes = tri->Nodes();
            cout<<nodes.Upper()<<endl;
          }



      }





  } // end method

} // end namespace


