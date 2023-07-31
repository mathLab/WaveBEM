#include "occ_line_smoothing.h"

#include <fstream>
#include <map>
#include <deal.II/base/logstream.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>

#include <IGESControl_Reader.hxx>
#include <IGESControl_Controller.hxx>
#include <IGESControl_Writer.hxx>
#include <gp_Pnt.hxx>
#include <gp_Dir.hxx>
#include <gp_Ax2.hxx>
#include <GC_MakeCircle.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <GeomConvert_CompCurveToBSplineCurve.hxx>
#include <Geom_BoundedCurve.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Geom_BSplineCurve.hxx>
#include <GeomLib_Tool.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <GC_MakeArcOfCircle.hxx>
#include <ShapeAnalysis_Curve.hxx>
#include <GeomAPI_IntCS.hxx>
#include <BRepAdaptor_Curve.hxx>

// all include files you need here

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/numerics/data_out.h>

#include "occ_utilities.h"

using namespace dealii;
using namespace std;
using namespace OpenCascade;


namespace OpenCascade
{
  LineSmoothing::LineSmoothing(Vector<double> &euler_vector,
                               Handle(Geom_Curve) ref_curve,
                               TopLoc_Location *curr_loc,
                               const DoFHandler<2,3> &dh,
                               const vector<bool> &smoothing_dofs,
                               const unsigned int b_point_id,
                               const unsigned int d_point_id,
                               const double tolerance) :
    euler_vector(euler_vector),
    ref_curve(ref_curve),
    curr_location(curr_loc),
    curve(ref_curve),
    projection(BRepBuilderAPI_MakeEdge(ref_curve)),
    dh(dh),
    smoothing_dofs(smoothing_dofs),
    tolerance(tolerance)
  {
    if (curr_location)
      {
        ref_location = *curr_location;
        used_location = *curr_location;
      }
    else
      {
        ref_location.Identity();
        used_location.Identity();
      }
    update_reference(b_point_id, d_point_id);
  }

  void LineSmoothing::update_reference(unsigned int bid, unsigned int did)
  {
    base_point_id = bid;
    driving_point_id = did;

    smoothing_list.clear();

    TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(ref_curve);
    edge.Location(ref_location);
    BRepAdaptor_Curve AC(edge);
    //GeomAdaptor_Curve AC(ref_curve);


    // Compute the reference support
    // points.
    support_points.resize(dh.n_dofs());
    //cout<<"dimension "<<dh.n_dofs()<<endl;
    DoFTools::map_dofs_to_support_points<2,3>(StaticMappingQ1<2,3>::mapping,
                                              dh, support_points);

    // Store the base point
    // parameter.
    gp_Pnt proj;
    ShapeAnalysis_Curve curve_analysis;
    double off = curve_analysis.Project(AC, Pnt(support_points[3*base_point_id]), tolerance, proj, occ_base_t, Standard_True);
    //cout<<"base_point_id "<<base_point_id<<" (of "<<dh.n_dofs()<<") "<<support_points[3*base_point_id]<<" vs "<<Pnt(proj)<<"  off: "<<off<<endl;
    AssertThrow( (off < tolerance), ExcMessage("Base point is not on curve."));

    // Store all parameters and
    // distances from the base point.
//IGESControl_Controller::Init();
//IGESControl_Writer ICW ("MM", 0);
//Standard_Boolean ok = ICW.AddGeom (AC);
//ICW.ComputeModel();
//Standard_Boolean OK = ICW.Write ("MyFile2.igs");

//cout<<"*"<<endl;
    for (unsigned int i=0; i<dh.n_dofs()/3; ++i)
      if (smoothing_dofs[3*i] == true)
        { 
          Point<3> p = support_points[3*i];
          //cout<<"Point: "<<p<<endl;
          gp_Pnt P = Pnt(p);
          double t;
          off = curve_analysis.Project(AC, P, tolerance, proj, t, Standard_True);
          //cout<<"Proj: "<<Pnt(proj)<<endl;
          //cout<<"in "<<3*i<<"-> p("<<support_points[3*i]<<")   proj="<<Pnt(proj)<<"  off "<<off<<endl;
          AssertThrow( (off < tolerance), ExcMessage("Point is not on ref curve!"));
          double dist = GCPnts_AbscissaPoint::Length(AC,t,occ_base_t);
          smoothing_list[std::pair<double, double >(t, dist)]=i;
          //cout<<"in "<<3*i<<"-> p("<<support_points[3*i]<<")   t="<<t<<"  dist="<<dist<<"  off "<<off<<endl;
        }

    // Store the new Reference length
    iterator endmap = smoothing_list.end();
    endmap--;
    ref_L = endmap->first.second;

    lengths_before_smoothing.reinit(smoothing_list.size());
    lengths_after_smoothing.reinit(smoothing_list.size());
    node_indices.resize(smoothing_list.size());
    fixed_length_ratios.reinit(smoothing_list.size());

    unsigned int count = 0;
    for (iterator it = smoothing_list.begin(); it != smoothing_list.end(); ++it)
      {
        fixed_length_ratios(count) = (it->first.second)/ref_L;
        node_indices[count] = it->second;
        //cout<<count<<"  "<<node_indices[count]<<"  "<<fixed_length_ratios(count)<<endl;
        count++;
      }


    // Sanity check: the base point
    // has to have zero length
    AssertThrow( smoothing_list.begin()->first.second < tolerance, ExcInternalError());

  }

  void LineSmoothing::smooth(bool maintain_on_original_curve)
  {
    Point<3> dP0 = support_points[driving_point_id*3], dP;
    for (unsigned int i=0; i<3; ++i)
      dP0(i) += euler_vector(3*driving_point_id+i);
    //cout<<driving_point_id<<endl;

    if (maintain_on_original_curve == true)
      {
        curve = ref_curve;

        ////////////////////////////////////////////////////////////////
        // general procedure
        //projection.normal_projection(dP, dP0);
        ////////////////////////////////////////////////////////////////

        // procedure only needed for waveBem
        used_location = *curr_location;
        //this is the horizontal plane
        TopLoc_Location L_inv = used_location.Inverted();

        Handle(Geom_Plane) horPlane = new Geom_Plane(0.,0.,1.,-dP0(2));
        horPlane->Transform(L_inv.Transformation());
        GeomAPI_IntCS Intersector(ref_curve, horPlane);
        int npoints = Intersector.NbPoints();
        AssertThrow((npoints != 0), ExcMessage("Reference curve is not intersecting with horizontal plane!"));
        double minDistance=1e7;
        for (int i=0; i<npoints; ++i)
          {
            gp_Pnt inters_point = Intersector.Point(i+1);
            inters_point.Transform(used_location.Transformation());
            Point<3> inters = Pnt(inters_point);
            if (dP0.distance(inters) < minDistance)
              {
                minDistance = dP0.distance(inters);
                dP = inters;
              }
          }
        //*/
        ////////////////////////////////////////////////////////////////
      }
    else
      {
        dP = dP0;

        std::vector<Point<3> > current_points;
        Point<3> direction(1.0,0.0,0.0);
        for (iterator it = smoothing_list.begin(); it != smoothing_list.end(); ++it)
          {
            unsigned int id=it->second;
            Point<3> p = support_points[3*id];
            for (unsigned int d=0; d<3; ++d)
              p(d) += euler_vector(3*id+d);
            current_points.push_back(p);
          }
        // Get the current curve
        curve = interpolation_curve(current_points);
        used_location.Identity();

        Point<3> base_point;
        base_point = support_points[base_point_id*3];
        for (unsigned int i=0; i<3; ++i)
          base_point(i) += euler_vector(3*base_point_id+i);
        gp_Pnt occ_driving_point = Pnt(dP);

        gp_Pnt proj;
        ShapeAnalysis_Curve curve_analysis;
        double off;
        off = curve_analysis.Project(curve, Pnt(base_point), tolerance, proj, occ_base_t, Standard_True);
        //cout<<"Point: "<<support_points[base_point_id*3]<<"    Proj: "<<Pnt(proj)<<"  Off: "<<off<<endl;
        AssertThrow( (off < tolerance), ExcMessage("Base point is not on curve."));


        // These lines can be used to dump the bslpine on an .igs file
        //IGESControl_Controller::Init();
        //IGESControl_Writer ICW ("MM", 0);
        //Standard_Boolean ok = ICW.AddGeom (curve);
        //ICW.ComputeModel();
        //Standard_Boolean OK = ICW.Write ("MyCurve.igs");

      }
    // Create a geometry adaptor on
    // the current curve
    TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(curve);
    edge.Location(used_location);
    BRepAdaptor_Curve AC(edge);
    Point<3> driving_point =dP;

    double off;
    gp_Pnt proj;
    ShapeAnalysis_Curve curve_analysis;
    off = curve_analysis.Project(AC, Pnt(driving_point), tolerance, proj, occ_driving_t, Standard_True);
    AssertThrow((off < tolerance), ExcMessage("Driving point is not on curve!"));

    // Total length
    double L = GCPnts_AbscissaPoint::Length(AC,occ_driving_t,occ_base_t);

    double direction = occ_driving_t > occ_base_t ? 1.0 : -1.0;

    // Now perform the smoothing
    unsigned int counter = 0;
    for (iterator it = smoothing_list.begin(); it != smoothing_list.end(); ++it)
      {
        unsigned int i=it->second;
        double dist=it->first.second;
        Point<3> p = support_points[3*i];
        for (unsigned int d=0; d<3; ++d)
          p(d) += euler_vector(3*i+d);
        // we now must compute where the old p lies on the curve (its parameter)
        // to be able to record the old length, useful for possible
        // interpolation of functions which might reconstruct the
        // new function values at the nodes after they have been moved by smoothing
        double old_t;
        off = curve_analysis.Project(AC, Pnt(p), tolerance, proj, old_t, Standard_True);
        //cout<<"Smoothing point distance from curve: "<<off<<"  Point: "<<p<<" vs "<<Pnt(proj)<<endl;
        AssertThrow( (off <L), ExcMessage("Smoothing point waaay off the curve."));
        lengths_before_smoothing(counter) = GCPnts_AbscissaPoint::Length(AC,old_t,occ_base_t);
        node_indices[counter] = i;
        //cout<<"out "<<3*i<<"--> p("<<p<<") t="<<t<<"  dist="<<dist<<endl;
        double new_dist = direction*fixed_length_ratios(counter)*L;
        lengths_after_smoothing(counter) = new_dist;
        GCPnts_AbscissaPoint AP(AC, new_dist, occ_base_t);
        double new_t = AP.Parameter();
        Point<3> pnew = Pnt(AC.Value(new_t));
        //cout<<"pnew("<<pnew<<") tnew="<<new_t<<"  distnew="<<new_dist<<endl;
        // Euler distance
        pnew -= support_points[3*i];
        for (unsigned int j=0; j<3; ++j)
          {
            //cout<<3*i+j<<" "<<pnew(j)<<endl;
            euler_vector(3*i+j) = pnew(j);
          }
        counter++;
      }
  }

}


void LineSmoothing::get_curve_tangent_vectors_at_smoothing_dofs(Vector<double> &tangents)
{
  AssertThrow(tangents.size()==dh.n_dofs(),
              ExcMessage("Tangent vector has wrong size"));

  TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(curve);
  edge.Location(used_location);
  BRepAdaptor_Curve AC(edge);

  for (iterator it = smoothing_list.begin(); it != smoothing_list.end(); ++it)
    {
      unsigned int i=it->second;
      double t=it->first.first;
      gp_Pnt P;
      gp_Vec V1;
      AC.D1(t,P,V1);

      tangents[3*i] = V1.X();
      tangents[3*i+1] = V1.Y();
      tangents[3*i+2] = V1.Z();
    }

}

void LineSmoothing::get_curve_length_ratios_at_smoothing_dofs(Vector<double> &length_ratios)
{
  AssertThrow(length_ratios.size()==dh.n_dofs(),
              ExcMessage("Length ratios vector has wrong size"));
  // Create a geometry adaptor on
  // the current curve
  TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(curve);
  edge.Location(used_location);
  BRepAdaptor_Curve AC(edge);
  // Total length
  double L = GCPnts_AbscissaPoint::Length(AC,occ_driving_t,occ_base_t);

  for (iterator it = smoothing_list.begin(); it != smoothing_list.end(); ++it)
    {
      unsigned int i=it->second;
      double dist = it->first.second;

      for (unsigned int j=0; j<3; ++j)
        length_ratios(3*i+j) = dist/L;
    }
}
