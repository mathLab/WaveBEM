#ifndef occ_axis_projection_h
#define occ_axis_projection_h

#include <deal.II/base/point.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/grid/tria_boundary.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <TopoDS.hxx>
#include <Geom_Curve.hxx>
#include <gp_Dir.hxx>

#include <deal.II/base/point.h>
#include <deal.II/grid/tria_boundary.h>



#include <stdio.h>
#include <stdlib.h>

using namespace dealii;
namespace OpenCascade 
{
  
  class AxisProjection : public  StraightBoundary<2,3> {
    public:
      AxisProjection(const TopoDS_Shape &sh, Point<3> direction, double tolerance=1e-7, double recovery_tolerance=1e-7);

      virtual Point<3> get_new_point_on_line
      (const Triangulation< 2,3 >::line_iterator &line) const;
    
      virtual Point<3> 	get_new_point_on_quad
      (const Triangulation< 2,3 >::quad_iterator &quad) const;

      virtual Point<3> 	project_to_surface
      (const Triangulation< 2,3 >::quad_iterator &quad, const Point<3> &y) const;
    
      const TopoDS_Shape &sh;
    
      bool axis_projection(Point<3> &projection,
			   const Point<3> &origin) const;

      bool assigned_axis_projection(Point<3> &projection,
			            const Point<3> &origin,
                                    const Point<3> &assigned_axis) const;


      bool assigned_axis_projection_and_diff_forms(Point<3> &projection,
                                                   Point<3> &normal,
                                                   double &mean_curvature,
                                                   const Point<3> &origin,
                                                   const Point<3> &assigned_axis) const;


      bool axis_projection_and_diff_forms(Point<3> &projection,
                                          Point<3> &normal,
                                          double &mean_curvature,
                                          const Point<3> &origin) const; 

    

    private:
      gp_Dir direction;
      Point<3> Direction;
      double tolerance;
      double recovery_tolerance;


  };


}

#endif

