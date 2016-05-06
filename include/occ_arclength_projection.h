#ifndef occ_arclength_projection_h
#define occ_arclength_projection_h

#include <deal.II/base/point.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/grid/tria_boundary.h>

#include <TopoDS.hxx>
#include <Geom_Curve.hxx>

#include <deal.II/base/point.h>
#include <deal.II/grid/tria_boundary.h>


#include <stdio.h>
#include <stdlib.h>

using namespace dealii;

namespace OpenCascade
{

  class ArclengthProjection : public  StraightBoundary<2,3>
  {
  public:

    ArclengthProjection(const TopoDS_Shape &sh,
                        double tolerance=1e-7);

    virtual Point<3> get_new_point_on_line
    (const Triangulation< 2,3 >::line_iterator &line) const;

    Point<3> arclength_projection(const Point<3> &p0,
                                  const Point<3> &p1,
                                  const double distance=.5) const;

    const TopoDS_Shape sh;


  private:
    const double tolerance;
  };

}


#endif
