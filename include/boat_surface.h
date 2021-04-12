
#ifndef boat_surface_h
#define boat_surface_h

#include <deal.II/base/point.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/grid/manifold.h>


using namespace dealii;

template <int dim>
class BoatSurface : public FlatManifold<dim-1,dim>
{
public:

  BoatSurface();

  void declare_parameters(ParameterHandler &prm);

  void parse_parameters(ParameterHandler &prm);

  double HullFunction(const Point<dim> point) const;

  Point<dim> HullNormal(const Point<dim> point) const;

  double HullMeanCurvature(const Point<dim> point) const;

  virtual Point<dim> get_new_point_on_line
  (const typename Triangulation< dim-1,dim >::line_iterator &line) const;

  virtual Point<dim>  get_new_point_on_quad
  (const typename Triangulation< dim-1,dim >::quad_iterator &quad) const;

private:

};

#endif
