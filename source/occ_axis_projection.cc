#include "occ_axis_projection.h"
#include "occ_utilities.h"

#include <iostream>

#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <Standard_Real.hxx>
#include <Standard_Integer.hxx>
#include <BRep_Tool.hxx>
#include <Geom_Surface.hxx>
#include <Prs3d_ShapeTool.hxx>
#include <gp_Ax3.hxx>
#include <gp_Lin.hxx>
#include <GeomLProp_SLProps.hxx>
#include <IntCurvesFace_ShapeIntersector.hxx>
#include <vector>
#define FALSE 0
#define TRUE 1


#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>


namespace OpenCascade
{


  AxisProjection::AxisProjection(const TopoDS_Shape &sh,
                                 Tensor<1,3> direction,
                                 double tolerance,
                                 double recovery_tolerance) :
    sh(sh),
    direction(direction[0],
              direction[1],
              direction[2]),
    Direction(direction),
    tolerance(tolerance),
    recovery_tolerance(recovery_tolerance)
  {
    // Check that we have at least
    // one face.
    TopExp_Explorer faceExplorer(sh , TopAbs_FACE);
    unsigned int n_faces = 0;
    while (faceExplorer.More())
      {
        n_faces++;
        faceExplorer.Next();
      }
    AssertThrow(n_faces > 0, ExcMessage("We have no faces to process"));


  }


  bool AxisProjection::axis_projection_and_diff_forms(Point<3> &projection,
                                                      Tensor<1,3> &normal,
                                                      double &mean_curvature,
                                                      const Point<3> &origin) const
  {
  
  cout<<"******* "<<origin<<endl;
    // translating original
    // Point<dim> to gp point
    gp_Pnt P0 = Pnt(origin);
    gp_Ax1 gpaxis(P0, direction);
    gp_Lin line(gpaxis);

    // destination point, normal
    // and mean curvature
    gp_Pnt Pproj(0.0,0.0,0.0);
    gp_Dir Normal;
    Standard_Real Mean_Curvature;

    // we prepare now the surface
    // for the projection we get
    // the whole shape from the
    // iges model
    IntCurvesFace_ShapeIntersector Inters;
    Inters.Load(sh,tolerance);
    Inters.Perform(line,-RealLast(),+RealLast());

    Point<3> average(0.0,0.0,0.0);
    Point<3> av_normal(0.0,0.0,0.0);
    double av_curvature = 0.0;
    if (Inters.NbPnt() == 0)
      {
        unsigned int succeeded = 0;
        cout<<"Axis A("<<Direction<<") direction projection of point P("<<origin<<")  on shape FAILED!"<<endl;
        cout<<"Trying to fix this"<<endl;

        gp_Pnt P1 = Pnt(origin+recovery_tolerance*Point<3>(1.0,0.0,0.0));
        gp_Ax1 gpaxis1(P1, direction);
        gp_Lin line1(gpaxis1);
        IntCurvesFace_ShapeIntersector Inters1;
        Inters1.Load(sh,tolerance);
        Inters1.Perform(line1,-RealLast(),+RealLast());
        //AssertThrow(Inters1.NbPnt() > 0,
        //   ExcMessage("Recovery point 1 projection on surface in given direction does not exist"));
        if (Inters1.NbPnt() > 0)
          {
            succeeded++;
            Standard_Real min_distance = 1e15;
            Standard_Real distance;
            int lowest_dist_int=0;
            for (int i=0; i<Inters1.NbPnt(); ++i)
              {
                distance = Pnt(origin).Distance(Inters1.Pnt(i+1));
                //cout<<"Point "<<i<<": "<<Pnt(Inters1.Pnt(i+1))<<"  distance: "<<distance<<endl;
                if (distance < min_distance)
                  {
                    min_distance = distance;
                    Pproj = Inters1.Pnt(i+1);
                    lowest_dist_int = i+1;
                  }
              }
            average += Pnt(Inters1.Pnt(lowest_dist_int));
            TopoDS_Face face1 =  Inters1.Face(lowest_dist_int);
            Handle(Geom_Surface) SurfToProj1 = BRep_Tool::Surface(face1);
            GeomLProp_SLProps props1(SurfToProj1, Inters1.UParameter(lowest_dist_int), Inters1.VParameter(lowest_dist_int), 1, tolerance);
            gp_Dir Normal1 = props1.Normal();
            Standard_Real Mean_Curvature1 = props1.MeanCurvature();
            // adjusting normal orientation
            if (face1.Orientation()==TopAbs_REVERSED)
              {
                Normal1.SetCoord(-1.0*Normal1.X(),-1.0*Normal1.Y(),-1.0*Normal1.Z());
                Mean_Curvature1*=-1.0;
              }
            av_normal += Point<3>(Normal1.X(),Normal1.Y(),Normal1.Z());
            av_curvature += Mean_Curvature1;
          }

        gp_Pnt P2 = Pnt(origin+recovery_tolerance*Point<3>(-1.0,0.0,0.0));
        gp_Ax1 gpaxis2(P2, direction);
        gp_Lin line2(gpaxis2);
        IntCurvesFace_ShapeIntersector Inters2;
        Inters2.Load(sh,tolerance);
        Inters2.Perform(line2,-RealLast(),+RealLast());
        //AssertThrow(Inters2.NbPnt() > 0,
        //  ExcMessage("Recovery point 2 projection on surface in given direction does not exist"));
        if (Inters2.NbPnt() > 0)
          {
            succeeded++;
            Standard_Real min_distance = 1e15;
            Standard_Real distance;
            int lowest_dist_int=0;
            for (int i=0; i<Inters2.NbPnt(); ++i)
              {
                distance = Pnt(origin).Distance(Inters2.Pnt(i+1));
                //cout<<"Point "<<i<<": "<<Pnt(Inters2.Pnt(i+1))<<"  distance: "<<distance<<endl;
                if (distance < min_distance)
                  {
                    min_distance = distance;
                    Pproj = Inters2.Pnt(i+1);
                    lowest_dist_int = i+1;
                  }
              }
            average += Pnt(Inters2.Pnt(lowest_dist_int));
            TopoDS_Face face2 =  Inters2.Face(lowest_dist_int);
            Handle(Geom_Surface) SurfToProj2 = BRep_Tool::Surface(face2);
            GeomLProp_SLProps props2(SurfToProj2, Inters2.UParameter(lowest_dist_int), Inters2.VParameter(lowest_dist_int), 1, tolerance);
            gp_Dir Normal2 = props2.Normal();
            Standard_Real Mean_Curvature2 = props2.MeanCurvature();
            // adjusting normal orientation
            if (face2.Orientation()==TopAbs_REVERSED)
              {
                Normal2.SetCoord(-1.0*Normal2.X(),-1.0*Normal2.Y(),-1.0*Normal2.Z());
                Mean_Curvature2*=-1.0;
              }
            av_normal += Point<3>(Normal2.X(),Normal2.Y(),Normal2.Z());
            av_curvature += Mean_Curvature2;
          }

        gp_Pnt P3 = Pnt(origin+recovery_tolerance*Point<3>(0.0,0.0,1.0));
        gp_Ax1 gpaxis3(P3, direction);
        gp_Lin line3(gpaxis3);
        IntCurvesFace_ShapeIntersector Inters3;
        Inters3.Load(sh,tolerance);
        Inters3.Perform(line3,-RealLast(),+RealLast());
        //AssertThrow(Inters3.NbPnt() > 0,
        //   ExcMessage("Recovery point 3 projection on surface in given direction does not exist"));
        if (Inters3.NbPnt() > 0)
          {
            succeeded++;
            Standard_Real min_distance = 1e15;
            Standard_Real distance;
            int lowest_dist_int=0;
            for (int i=0; i<Inters3.NbPnt(); ++i)
              {
                distance = Pnt(origin).Distance(Inters3.Pnt(i+1));
                //cout<<"Point "<<i<<": "<<Pnt(Inters3.Pnt(i+1))<<"  distance: "<<distance<<endl;
                if (distance < min_distance)
                  {
                    min_distance = distance;
                    Pproj = Inters3.Pnt(i+1);
                    lowest_dist_int = i+1;
                  }
              }
            average += Pnt(Inters3.Pnt(lowest_dist_int));
            TopoDS_Face face3 =  Inters3.Face(lowest_dist_int);
            Handle(Geom_Surface) SurfToProj3 = BRep_Tool::Surface(face3);
            GeomLProp_SLProps props3(SurfToProj3, Inters3.UParameter(lowest_dist_int), Inters3.VParameter(lowest_dist_int), 1, tolerance);
            gp_Dir Normal3 = props3.Normal();
            Standard_Real Mean_Curvature3 = props3.MeanCurvature();
            // adjusting normal orientation
            if (face3.Orientation()==TopAbs_REVERSED)
              {
                Normal3.SetCoord(-1.0*Normal3.X(),-1.0*Normal3.Y(),-1.0*Normal3.Z());
                Mean_Curvature3*=-1.0;
              }
            av_normal += Point<3>(Normal3.X(),Normal3.Y(),Normal3.Z());
            av_curvature += Mean_Curvature3;
          }

        gp_Pnt P4 = Pnt(origin+recovery_tolerance*Point<3>(0.0,0.0,-1.0));
        gp_Ax1 gpaxis4(P4, direction);
        gp_Lin line4(gpaxis4);
        IntCurvesFace_ShapeIntersector Inters4;
        Inters4.Load(sh,tolerance);
        Inters4.Perform(line4,-RealLast(),+RealLast());
        //AssertThrow(Inters4.NbPnt() > 0,
        //  ExcMessage("Recovery point 4 projection on surface in given direction does not exist"));
        if (Inters4.NbPnt() > 0)
          {
            succeeded++;
            Standard_Real min_distance = 1e15;
            Standard_Real distance;
            int lowest_dist_int=0;
            for (int i=0; i<Inters4.NbPnt(); ++i)
              {
                distance = Pnt(origin).Distance(Inters4.Pnt(i+1));
                //cout<<"Point "<<i<<": "<<Pnt(Inters4.Pnt(i+1))<<"  distance: "<<distance<<endl;
                if (distance < min_distance)
                  {
                    min_distance = distance;
                    Pproj = Inters4.Pnt(i+1);
                    lowest_dist_int = i+1;
                  }
              }
            average += Pnt(Inters4.Pnt(lowest_dist_int));
            TopoDS_Face face4 =  Inters4.Face(lowest_dist_int);
            Handle(Geom_Surface) SurfToProj4 = BRep_Tool::Surface(face4);
            GeomLProp_SLProps props4(SurfToProj4, Inters4.UParameter(lowest_dist_int), Inters4.VParameter(lowest_dist_int), 1, tolerance);
            gp_Dir Normal4 = props4.Normal();
            Standard_Real Mean_Curvature4 = props4.MeanCurvature();
            // adjusting normal orientation
            if (face4.Orientation()==TopAbs_REVERSED)
              {
                Normal4.SetCoord(-1.0*Normal4.X(),-1.0*Normal4.Y(),-1.0*Normal4.Z());
                Mean_Curvature4*=-1.0;
              }
            av_normal += Point<3>(Normal4.X(),Normal4.Y(),Normal4.Z());
            av_curvature += Mean_Curvature4;
          }
        if (succeeded > 0)
          {
            cout<<"Recovery attempt of point projection on surface in given direction FAILED"<<endl;
            return false;
          }
        //AssertThrow(succeeded > 0,
        //   ExcMessage("Recovery attempt of point projection on surface in given direction FAILED"));

        Pproj = Pnt(average/succeeded);
        av_normal/=succeeded;
        Normal.SetCoord(av_normal(0),av_normal(1),av_normal(2));
        Mean_Curvature = av_curvature/succeeded;
      }
    else
      {
        //cout<<"Intersections found: "<<Inters.NbPnt()<<endl;
        Standard_Real min_distance = 1e15;
        Standard_Real distance;
        int lowest_dist_int=0;
        for (int i=0; i<Inters.NbPnt(); ++i)
          {
            distance = Pnt(origin).Distance(Inters.Pnt(i+1));
            //cout<<"Point "<<i<<": "<<Pnt(Inters.Pnt(i+1))<<"  distance: "<<distance<<endl;
            if (distance < min_distance)
              {
                min_distance = distance;
                Pproj = Inters.Pnt(i+1);
                lowest_dist_int = i+1;
              }
          }
        TopoDS_Face face =  Inters.Face(lowest_dist_int);
        Handle(Geom_Surface) SurfToProj = BRep_Tool::Surface(face);
        GeomLProp_SLProps props(SurfToProj, Inters.UParameter(lowest_dist_int), Inters.VParameter(lowest_dist_int), 1, tolerance);
        Normal = props.Normal();
        Mean_Curvature = props.MeanCurvature();
        // adjusting normal orientation
        if (face.Orientation()==TopAbs_REVERSED)
          {
            Normal.SetCoord(-1.0*Normal.X(),-1.0*Normal.Y(),-1.0*Normal.Z());
            Mean_Curvature*=-1.0;
          }
      }


// translating destination point
    projection = Pnt(Pproj);

// translating normal vector
    normal[0] = Normal.X();
    normal[1] = Normal.Y();
    normal[2] = Normal.Z();


// translating mean curvature
    mean_curvature = double(Mean_Curvature);

    return true;
  }

  bool AxisProjection::axis_projection(Point<3> &projection,
                                       const Point<3> &origin) const
  {
    // translating original
    // Point<dim> to gp point


    gp_Pnt P0 = Pnt(origin);
    gp_Ax1 gpaxis(P0, direction);
    gp_Lin line(gpaxis);

    // destination point
    gp_Pnt Pproj(0.0,0.0,0.0);

    // we prepare now the surface
    // for the projection we get
    // the whole shape from the
    // iges model
    IntCurvesFace_ShapeIntersector Inters;
    Inters.Load(sh,tolerance);

    Inters.Perform(line,-RealLast(),+RealLast());

    Point<3> average(0.0,0.0,0.0);
    if (Inters.NbPnt() == 0)
      {
        unsigned int succeeded = 0;
        cout<<"Axis A("<<Direction<<") direction projection of point P("<<origin<<")  on shape FAILED!"<<endl;
        cout<<"Trying to fix this"<<endl;

        gp_Pnt P1 = Pnt(origin+recovery_tolerance*Point<3>(1.0,0.0,0.0));
        gp_Ax1 gpaxis1(P1, direction);
        gp_Lin line1(gpaxis1);
        IntCurvesFace_ShapeIntersector Inters1;
        Inters1.Load(sh,tolerance);
        Inters1.Perform(line1,-RealLast(),+RealLast());
        //AssertThrow(Inters1.NbPnt() > 0,
        //  ExcMessage("Recovery point 1 projection on surface in given direction does not exist"));
        if (Inters1.NbPnt() > 0)
          {
            Standard_Real min_distance = 1e15;
            Standard_Real distance;
            int lowest_dist_int=0;
            for (int i=0; i<Inters1.NbPnt(); ++i)
              {
                distance = Pnt(origin).Distance(Inters1.Pnt(i+1));
                //cout<<"Point "<<i<<": "<<Pnt(Inters1.Pnt(i+1))<<"  distance: "<<distance<<endl;
                if (distance < min_distance)
                  {
                    min_distance = distance;
                    Pproj = Inters1.Pnt(i+1);
                    lowest_dist_int = i+1;
                  }
              }
            average += Pnt(Inters1.Pnt(lowest_dist_int));
            succeeded++;
          }

        gp_Pnt P2 = Pnt(origin+recovery_tolerance*Point<3>(-1.0,0.0,0.0));
        gp_Ax1 gpaxis2(P2, direction);
        gp_Lin line2(gpaxis2);
        IntCurvesFace_ShapeIntersector Inters2;
        Inters2.Load(sh,tolerance);
        Inters2.Perform(line2,-RealLast(),+RealLast());
        //AssertThrow(Inters2.NbPnt() > 0,
        //  ExcMessage("Recovery point 2 projection on surface in given direction does not exist"));
        if (Inters2.NbPnt() > 0)
          {
            Standard_Real min_distance = 1e15;
            Standard_Real distance;
            int lowest_dist_int=0;
            for (int i=0; i<Inters2.NbPnt(); ++i)
              {
                distance = Pnt(origin).Distance(Inters2.Pnt(i+1));
                //cout<<"Point "<<i<<": "<<Pnt(Inters2.Pnt(i+1))<<"  distance: "<<distance<<endl;
                if (distance < min_distance)
                  {
                    min_distance = distance;
                    Pproj = Inters2.Pnt(i+1);
                    lowest_dist_int = i+1;
                  }
              }
            average += Pnt(Inters2.Pnt(lowest_dist_int));
            succeeded++;
          }

        gp_Pnt P3 = Pnt(origin+recovery_tolerance*Point<3>(0.0,0.0,1.0));
        gp_Ax1 gpaxis3(P3, direction);
        gp_Lin line3(gpaxis3);
        IntCurvesFace_ShapeIntersector Inters3;
        Inters3.Load(sh,tolerance);
        Inters3.Perform(line3,-RealLast(),+RealLast());
        //AssertThrow(Inters3.NbPnt() > 0,
        //   ExcMessage("Recovery point 3 projection on surface in given direction does not exist"));
        if (Inters3.NbPnt() > 0)
          {
            Standard_Real min_distance = 1e15;
            Standard_Real distance;
            int lowest_dist_int=0;
            for (int i=0; i<Inters3.NbPnt(); ++i)
              {
                distance = Pnt(origin).Distance(Inters3.Pnt(i+1));
                //cout<<"Point "<<i<<": "<<Pnt(Inters3.Pnt(i+1))<<"  distance: "<<distance<<endl;
                if (distance < min_distance)
                  {
                    min_distance = distance;
                    Pproj = Inters3.Pnt(i+1);
                    lowest_dist_int = i+1;
                  }
              }
            average += Pnt(Inters3.Pnt(lowest_dist_int));
            succeeded++;
          }

        gp_Pnt P4 = Pnt(origin+recovery_tolerance*Point<3>(0.0,0.0,-1.0));
        gp_Ax1 gpaxis4(P4, direction);
        gp_Lin line4(gpaxis4);
        IntCurvesFace_ShapeIntersector Inters4;
        Inters4.Load(sh,tolerance);
        Inters4.Perform(line4,-RealLast(),+RealLast());
        //AssertThrow(Inters4.NbPnt() > 0,
        //  ExcMessage("Recovery point 4 projection on surface in given direction does not exist"));
        if (Inters4.NbPnt() > 0)
          {
            Standard_Real min_distance = 1e15;
            Standard_Real distance;
            int lowest_dist_int=0;
            for (int i=0; i<Inters4.NbPnt(); ++i)
              {
                distance = Pnt(origin).Distance(Inters4.Pnt(i+1));
                //cout<<"Point "<<i<<": "<<Pnt(Inters4.Pnt(i+1))<<"  distance: "<<distance<<endl;
                if (distance < min_distance)
                  {
                    min_distance = distance;
                    Pproj = Inters4.Pnt(i+1);
                    lowest_dist_int = i+1;
                  }
              }
            average += Pnt(Inters4.Pnt(lowest_dist_int));
            succeeded++;
          }
        if (succeeded > 0)
          {
            cout<<"Recovery attempt of point projection on surface in given direction FAILED"<<endl;
            return false;
          }
        //AssertThrow(succeeded > 0,
        //   ExcMessage("Recovery attempt of point projection on surface in given direction FAILED"));

        Pproj = Pnt(average/succeeded);
      }
    else
      {
        //cout<<"Intersections found: "<<Inters.NbPnt()<<endl;
        Standard_Real min_distance = 1e15;
        Standard_Real distance;
        for (int i=0; i<Inters.NbPnt(); ++i)
          {
            distance = Pnt(origin).Distance(Inters.Pnt(i+1));
            //cout<<"Point "<<i<<": "<<Pnt(Inters.Pnt(i+1))<<"  distance: "<<distance<<endl;
            if (distance < min_distance)
              {
                min_distance = distance;
                Pproj = Inters.Pnt(i+1);
              }
          }
      }

// translating destination point
    projection = Pnt(Pproj);
    //cout<<"y dir projection: "<<projection<<endl;
    return true;
  }


  bool AxisProjection::assigned_axis_projection(Point<3> &projection,
                                                const Point<3> &origin,
                                                const Tensor<1,3> &assigned_axis) const
  {
    // translating original
    // Point<dim> to gp point
    gp_Pnt P0 = Pnt(origin);
    gp_Dir axis(assigned_axis[0],assigned_axis[1],assigned_axis[2]);
    gp_Ax1 gpaxis(P0, axis);
    gp_Lin line(gpaxis);

    // destination point
    gp_Pnt Pproj(0.0,0.0,0.0);
    // we prepare now the surface
    // for the projection we get
    // the whole shape from the
    // iges model
    IntCurvesFace_ShapeIntersector Inters;
    Inters.Load(sh,tolerance);
    Inters.Perform(line,-RealLast(),+RealLast());

    Point<3> average(0.0,0.0,0.0);
    if (Inters.NbPnt() == 0)
      {
        unsigned int succeeded = 0;
        cout<<"Axis A("<<assigned_axis <<") direction projection of point P("<<origin<<")  on shape FAILED!"<<endl;
        cout<<"Trying to fix this"<<endl;

        gp_Pnt P1 = Pnt(origin+recovery_tolerance*Point<3>(1.0,0.0,0.0));
        gp_Ax1 gpaxis1(P1, direction);
        gp_Lin line1(gpaxis1);
        IntCurvesFace_ShapeIntersector Inters1;
        Inters1.Load(sh,tolerance);
        Inters1.Perform(line1,-RealLast(),+RealLast());
        //AssertThrow(Inters1.NbPnt() > 0,
        //   ExcMessage("Recovery point 1 projection on surface in given direction does not exist"));
        if (Inters1.NbPnt() > 0)
          {
            Standard_Real min_distance = 1e15;
            Standard_Real distance;
            int lowest_dist_int=0;
            for (int i=0; i<Inters1.NbPnt(); ++i)
              {
                distance = Pnt(origin).Distance(Inters1.Pnt(i+1));
                //cout<<"Point "<<i<<": "<<Pnt(Inters1.Pnt(i+1))<<"  distance: "<<distance<<endl;
                if (distance < min_distance)
                  {
                    min_distance = distance;
                    Pproj = Inters1.Pnt(i+1);
                    lowest_dist_int = i+1;
                  }
              }
            succeeded++;
            average += Pnt(Inters1.Pnt(lowest_dist_int));
          }

        gp_Pnt P2 = Pnt(origin+recovery_tolerance*Point<3>(-1.0,0.0,0.0));
        gp_Ax1 gpaxis2(P2, direction);
        gp_Lin line2(gpaxis2);
        IntCurvesFace_ShapeIntersector Inters2;
        Inters2.Load(sh,tolerance);
        Inters2.Perform(line2,-RealLast(),+RealLast());
        //AssertThrow(Inters2.NbPnt() > 0,
        //   ExcMessage("Recovery point 2 projection on surface in given direction does not exist"));
        if (Inters2.NbPnt() > 0)
          {
            Standard_Real min_distance = 1e15;
            Standard_Real distance;
            int lowest_dist_int=0;
            for (int i=0; i<Inters2.NbPnt(); ++i)
              {
                distance = Pnt(origin).Distance(Inters2.Pnt(i+1));
                //cout<<"Point "<<i<<": "<<Pnt(Inters2.Pnt(i+1))<<"  distance: "<<distance<<endl;
                if (distance < min_distance)
                  {
                    min_distance = distance;
                    Pproj = Inters2.Pnt(i+1);
                    lowest_dist_int = i+1;
                  }
              }
            succeeded++;
            average += Pnt(Inters2.Pnt(lowest_dist_int));
          }
        gp_Pnt P3 = Pnt(origin+recovery_tolerance*Point<3>(0.0,0.0,1.0));
        gp_Ax1 gpaxis3(P3, direction);
        gp_Lin line3(gpaxis3);
        IntCurvesFace_ShapeIntersector Inters3;
        Inters3.Load(sh,tolerance);
        Inters3.Perform(line3,-RealLast(),+RealLast());
        //AssertThrow(Inters3.NbPnt() > 0,
        //  ExcMessage("Recovery point 3 projection on surface in given direction does not exist"));
        if (Inters3.NbPnt() > 0)
          {
            succeeded++;
            Standard_Real min_distance = 1e15;
            Standard_Real distance;
            int lowest_dist_int=0;
            for (int i=0; i<Inters3.NbPnt(); ++i)
              {
                distance = Pnt(origin).Distance(Inters3.Pnt(i+1));
                //cout<<"Point "<<i<<": "<<Pnt(Inters3.Pnt(i+1))<<"  distance: "<<distance<<endl;
                if (distance < min_distance)
                  {
                    min_distance = distance;
                    Pproj = Inters3.Pnt(i+1);
                    lowest_dist_int = i+1;
                  }
              }
            average += Pnt(Inters3.Pnt(lowest_dist_int));
          }
        gp_Pnt P4 = Pnt(origin+recovery_tolerance*Point<3>(0.0,0.0,-1.0));
        gp_Ax1 gpaxis4(P4, direction);
        gp_Lin line4(gpaxis4);
        IntCurvesFace_ShapeIntersector Inters4;
        Inters4.Load(sh,tolerance);
        Inters4.Perform(line4,-RealLast(),+RealLast());
        //AssertThrow(Inters4.NbPnt() > 0,
        //  ExcMessage("Recovery point 4 projection on surface in given direction does not exist"));
        if (Inters4.NbPnt() > 0)
          {
            succeeded++;
            Standard_Real min_distance = 1e15;
            Standard_Real distance;
            int lowest_dist_int=0;
            for (int i=0; i<Inters4.NbPnt(); ++i)
              {
                distance = Pnt(origin).Distance(Inters4.Pnt(i+1));
                //cout<<"Point "<<i<<": "<<Pnt(Inters4.Pnt(i+1))<<"  distance: "<<distance<<endl;
                if (distance < min_distance)
                  {
                    min_distance = distance;
                    Pproj = Inters4.Pnt(i+1);
                    lowest_dist_int = i+1;
                  }
              }
            average += Pnt(Inters4.Pnt(lowest_dist_int));
          }
        if (succeeded > 0)
          {
            cout<<"Recovery attempt of point projection on surface in given direction FAILED"<<endl;
            return false;
          }
        //AssertThrow(succeeded > 0,
        //     ExcMessage("Recovery attempt of point projection on surface in given direction FAILED"));

        Pproj = Pnt(average/succeeded);
      }
    else
      {
        //cout<<"Intersections found: "<<Inters.NbPnt()<<endl;
        Standard_Real min_distance = 1e15;
        Standard_Real distance;
        for (int i=0; i<Inters.NbPnt(); ++i)
          {
            distance = Pnt(origin).Distance(Inters.Pnt(i+1));
            //cout<<"Point "<<i<<": "<<Pnt(Inters.Pnt(i+1))<<"  distance: "<<distance<<endl;
            if (distance < min_distance)
              {
                min_distance = distance;
                Pproj = Inters.Pnt(i+1);
              }
          }
      }


    //Pproj = Inters.Pnt(1);
// translating destination point
    projection = Pnt(Pproj);
    return true;
  }



  bool AxisProjection::assigned_axis_projection_and_diff_forms(Point<3> &projection,
      Tensor<1,3> &normal,
      double &mean_curvature,
      const Point<3> &origin,
      const Tensor<1,3> &assigned_axis) const
  {
    // translating original
    // Point<dim> to gp point
    gp_Pnt P0 = Pnt(origin);
    gp_Dir axis(assigned_axis[0],assigned_axis[1],assigned_axis[2]);
    gp_Ax1 gpaxis(P0, axis);
    gp_Lin line(gpaxis);

    // destination point, normal
    // and mean curvature
    gp_Pnt Pproj(0.0,0.0,0.0);
    gp_Dir Normal;
    Standard_Real Mean_Curvature;

    // we prepare now the surface
    // for the projection we get
    // the whole shape from the
    // iges model
    IntCurvesFace_ShapeIntersector Inters;
    Inters.Load(sh,tolerance);
    //Inters.Load(sh,Precision::Confusion());
    Inters.Perform(line,-RealLast(),+RealLast());

    Point<3> average(0.0,0.0,0.0);
    Tensor<1,3> av_normal;
    double av_curvature = 0.0;
    if (Inters.NbPnt() == 0)
      {
        unsigned int succeeded = 0;
        cout<<"Axis A("<<assigned_axis <<") direction projection of point P("<<origin<<")  on shape FAILED!"<<endl;
        cout<<"Trying to fix this"<<endl;


        gp_Pnt P1 = Pnt(origin+recovery_tolerance*Point<3>(1.0,0.0,0.0));
        gp_Ax1 gpaxis1(P1, direction);
        gp_Lin line1(gpaxis1);
        IntCurvesFace_ShapeIntersector Inters1;
        Inters1.Load(sh,tolerance);
        Inters1.Perform(line1,-RealLast(),+RealLast());
        //AssertThrow(Inters1.NbPnt() > 0,
        //  ExcMessage("Recovery point 1 projection on surface in given direction does not exist"));
        if (Inters1.NbPnt() > 0)
          {
            succeeded++;
            Standard_Real min_distance = 1e15;
            Standard_Real distance;
            int lowest_dist_int=0;
            for (int i=0; i<Inters1.NbPnt(); ++i)
              {
                distance = Pnt(origin).Distance(Inters1.Pnt(i+1));
                //cout<<"Point "<<i<<": "<<Pnt(Inters1.Pnt(i+1))<<"  distance: "<<distance<<endl;
                if (distance < min_distance)
                  {
                    min_distance = distance;
                    Pproj = Inters1.Pnt(i+1);
                    lowest_dist_int = i+1;
                  }
              }
            average += Pnt(Inters1.Pnt(lowest_dist_int));
            TopoDS_Face face1 =  Inters1.Face(lowest_dist_int);
            Handle(Geom_Surface) SurfToProj1 = BRep_Tool::Surface(face1);
            GeomLProp_SLProps props1(SurfToProj1, Inters1.UParameter(lowest_dist_int), Inters1.VParameter(lowest_dist_int), 1, tolerance);
            gp_Dir Normal1 = props1.Normal();
            Standard_Real Mean_Curvature1 = props1.MeanCurvature();
            // adjusting normal orientation
            if (face1.Orientation()==TopAbs_REVERSED)
              {
                Normal1.SetCoord(-1.0*Normal1.X(),-1.0*Normal1.Y(),-1.0*Normal1.Z());
                Mean_Curvature1*=-1.0;
              }
            av_normal[0] += Normal1.X();
            av_normal[1] += Normal1.Y();
            av_normal[2] += Normal1.Z();

            av_curvature += Mean_Curvature1;
          }

        gp_Pnt P2 = Pnt(origin+recovery_tolerance*Point<3>(-1.0,0.0,0.0));
        gp_Ax1 gpaxis2(P2, direction);
        gp_Lin line2(gpaxis2);
        IntCurvesFace_ShapeIntersector Inters2;
        Inters2.Load(sh,tolerance);
        Inters2.Perform(line2,-RealLast(),+RealLast());
        //AssertThrow(Inters2.NbPnt() > 0,
        //  ExcMessage("Recovery point 2 projection on surface in given direction does not exist"));
        if (Inters2.NbPnt() > 0)
          {
            succeeded++;
            Standard_Real min_distance = 1e15;
            Standard_Real distance;
            int lowest_dist_int=0;
            for (int i=0; i<Inters2.NbPnt(); ++i)
              {
                distance = Pnt(origin).Distance(Inters2.Pnt(i+1));
                //cout<<"Point "<<i<<": "<<Pnt(Inters2.Pnt(i+1))<<"  distance: "<<distance<<endl;
                if (distance < min_distance)
                  {
                    min_distance = distance;
                    Pproj = Inters2.Pnt(i+1);
                    lowest_dist_int = i+1;
                  }
              }
            average += Pnt(Inters2.Pnt(lowest_dist_int));
            TopoDS_Face face2 =  Inters2.Face(lowest_dist_int);
            Handle(Geom_Surface) SurfToProj2 = BRep_Tool::Surface(face2);
            GeomLProp_SLProps props2(SurfToProj2, Inters2.UParameter(lowest_dist_int), Inters2.VParameter(lowest_dist_int), 1, tolerance);
            gp_Dir Normal2 = props2.Normal();
            Standard_Real Mean_Curvature2 = props2.MeanCurvature();
            // adjusting normal orientation
            if (face2.Orientation()==TopAbs_REVERSED)
              {
                Normal2.SetCoord(-1.0*Normal2.X(),-1.0*Normal2.Y(),-1.0*Normal2.Z());
                Mean_Curvature2*=-1.0;
              }
            av_normal[0] += Normal2.X();
            av_normal[1] += Normal2.Y();
            av_normal[2] += Normal2.Z();
            av_curvature += Mean_Curvature2;
          }

        gp_Pnt P3 = Pnt(origin+recovery_tolerance*Point<3>(0.0,0.0,1.0));
        gp_Ax1 gpaxis3(P3, direction);
        gp_Lin line3(gpaxis3);
        IntCurvesFace_ShapeIntersector Inters3;
        Inters3.Load(sh,tolerance);
        Inters3.Perform(line3,-RealLast(),+RealLast());
        //AssertThrow(Inters3.NbPnt() > 0,
        //  ExcMessage("Recovery point 3 projection on surface in given direction does not exist"));
        if (Inters3.NbPnt() > 0)
          {
            succeeded++;
            Standard_Real min_distance = 1e15;
            Standard_Real distance;
            int lowest_dist_int=0;
            for (int i=0; i<Inters3.NbPnt(); ++i)
              {
                distance = Pnt(origin).Distance(Inters3.Pnt(i+1));
                //cout<<"Point "<<i<<": "<<Pnt(Inters3.Pnt(i+1))<<"  distance: "<<distance<<endl;
                if (distance < min_distance)
                  {
                    min_distance = distance;
                    Pproj = Inters3.Pnt(i+1);
                    lowest_dist_int = i+1;
                  }
              }
            average += Pnt(Inters3.Pnt(lowest_dist_int));
            TopoDS_Face face3 =  Inters3.Face(lowest_dist_int);
            Handle(Geom_Surface) SurfToProj3 = BRep_Tool::Surface(face3);
            GeomLProp_SLProps props3(SurfToProj3, Inters3.UParameter(lowest_dist_int), Inters3.VParameter(lowest_dist_int), 1, tolerance);
            gp_Dir Normal3 = props3.Normal();
            Standard_Real Mean_Curvature3 = props3.MeanCurvature();
            // adjusting normal orientation
            if (face3.Orientation()==TopAbs_REVERSED)
              {
                Normal3.SetCoord(-1.0*Normal3.X(),-1.0*Normal3.Y(),-1.0*Normal3.Z());
                Mean_Curvature3*=-1.0;
              }
            av_normal[0] += Normal3.X();
            av_normal[1] += Normal3.Y();
            av_normal[2] += Normal3.Z();
            av_curvature += Mean_Curvature3;
          }

        gp_Pnt P4 = Pnt(origin+recovery_tolerance*Point<3>(0.0,0.0,-1.0));
        gp_Ax1 gpaxis4(P4, direction);
        gp_Lin line4(gpaxis4);
        IntCurvesFace_ShapeIntersector Inters4;
        Inters4.Load(sh,tolerance);
        Inters4.Perform(line4,-RealLast(),+RealLast());
        //AssertThrow(Inters4.NbPnt() > 0,
        //  ExcMessage("Recovery point 4 projection on surface in given direction does not exist"));
        if (Inters4.NbPnt() > 0)
          {
            succeeded++;
            Standard_Real min_distance = 1e15;
            Standard_Real distance;
            int lowest_dist_int=0;
            for (int i=0; i<Inters4.NbPnt(); ++i)
              {
                distance = Pnt(origin).Distance(Inters4.Pnt(i+1));
                //cout<<"Point "<<i<<": "<<Pnt(Inters4.Pnt(i+1))<<"  distance: "<<distance<<endl;
                if (distance < min_distance)
                  {
                    min_distance = distance;
                    Pproj = Inters4.Pnt(i+1);
                    lowest_dist_int = i+1;
                  }
              }
            average += Pnt(Inters4.Pnt(lowest_dist_int));
            TopoDS_Face face4 =  Inters4.Face(lowest_dist_int);
            Handle(Geom_Surface) SurfToProj4 = BRep_Tool::Surface(face4);
            GeomLProp_SLProps props4(SurfToProj4, Inters4.UParameter(lowest_dist_int), Inters4.VParameter(lowest_dist_int), 1, tolerance);
            gp_Dir Normal4 = props4.Normal();
            Standard_Real Mean_Curvature4 = props4.MeanCurvature();
            // adjusting normal orientation
            if (face4.Orientation()==TopAbs_REVERSED)
              {
                Normal4.SetCoord(-1.0*Normal4.X(),-1.0*Normal4.Y(),-1.0*Normal4.Z());
                Mean_Curvature4*=-1.0;
              }
            av_normal[0] += Normal4.X();
            av_normal[1] += Normal4.Y();
            av_normal[2] += Normal4.Z();
            av_curvature += Mean_Curvature4;
          }
        if (succeeded > 0)
          {
            cout<<"Recovery attempt of point projection on surface in given direction FAILED"<<endl;
            return false;
          }
        //AssertThrow(succeeded > 0,
        //   ExcMessage("Recovery attempt of point projection on surface in given direction FAILED"));

        Pproj = Pnt(average/succeeded);
        av_normal/=succeeded;
        Normal.SetCoord(av_normal[0],av_normal[1],av_normal[2]);
        Mean_Curvature = av_curvature/succeeded;
      }
    else
      {
        //cout<<"Intersections found: "<<Inters.NbPnt()<<endl;
        Standard_Real min_distance = 1e15;
        Standard_Real distance;
        int lowest_dist_int=0;
        for (int i=0; i<Inters.NbPnt(); ++i)
          {
            distance = Pnt(origin).Distance(Inters.Pnt(i+1));
            //cout<<"Point "<<i<<": "<<Pnt(Inters.Pnt(i+1))<<"  distance: "<<distance<<endl;
            if (distance < min_distance)
              {
                min_distance = distance;
                Pproj = Inters.Pnt(i+1);
                lowest_dist_int = i+1;
              }
          }
        TopoDS_Face face =  Inters.Face(lowest_dist_int);
        Handle(Geom_Surface) SurfToProj = BRep_Tool::Surface(face);
        GeomLProp_SLProps props(SurfToProj, Inters.UParameter(lowest_dist_int), Inters.VParameter(lowest_dist_int), 1, tolerance);
        Normal = props.Normal();
        Mean_Curvature = props.MeanCurvature();
        // adjusting normal orientation
        if (face.Orientation()==TopAbs_REVERSED)
          {
            Normal.SetCoord(-1.0*Normal.X(),-1.0*Normal.Y(),-1.0*Normal.Z());
            Mean_Curvature*=-1.0;
          }
      }


// translating destination point
    projection = Pnt(Pproj);

// translating normal vector
    normal[0] = Normal.X();
    normal[1] = Normal.Y();
    normal[2] = Normal.Z();

// translating mean curvature
    mean_curvature = double(Mean_Curvature);

    return true;
  }


  Point<3> AxisProjection::get_new_point_on_line
  (const Triangulation< 2,3 >::line_iterator &line) const
  {cout<<"Here? Do we ever pass by here?"<<endl;
    Point<3> projected_point;
    Point<3> source_point = FlatManifold<2,3>::get_new_point_on_line(line);
    axis_projection(projected_point, source_point);
    return projected_point;
  }


  Point<3> AxisProjection::get_new_point_on_quad
  (const Triangulation< 2,3 >::quad_iterator &quad) const
  {
    Point<3> projected_point;
    Point<3> source_point = FlatManifold<2,3>::get_new_point_on_quad(quad);

    axis_projection(projected_point, source_point);


    return projected_point;
  }

  Point<3> AxisProjection::project_to_surface
  (const Triangulation< 2,3 >::quad_iterator &quad, const Point<3>  &y) const
  {
    Point<3> projected_point;
    axis_projection(projected_point, y);

    return projected_point;
  }




}
