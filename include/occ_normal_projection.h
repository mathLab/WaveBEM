#ifndef occ_normal_projection_h
#define occ_normal_projection_h

#include <deal.II/base/point.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/grid/manifold.h>


#include <TopoDS.hxx>
#include <Geom_Curve.hxx>



#include <stdio.h>
#include <stdlib.h>

using namespace dealii;
namespace OpenCascade
{

  template <int dim>
  class NormalProjection : public  FlatManifold<2,3>
  {
  public:
    NormalProjection(const TopoDS_Shape &sh);

    virtual Point<3> get_new_point_on_line
    (const Triangulation< 2,3 >::line_iterator &line) const;

    virtual Point<3>  get_new_point_on_quad
    (const Triangulation< 2,3 >::quad_iterator &quad) const;

    const TopoDS_Shape &sh;

    void normal_projection(Point<3> &projection,
                           const Point<3> &origin) const;



    void normal_projection_and_diff_forms(Point<3> &projection,
                                          Point<3> &normal,
                                          double &mean_curvature,
                                          const Point<3> &origin) const;

  };
}

#endif
