#ifndef occ_line_smoothing_h
#define occ_line_smoothing_h

#include <string>
#include <TopoDS_Shape.hxx>
#include <Geom_Plane.hxx>
#include <GeomLib_Tool.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <gp_Pnt.hxx>

#include <deal.II/base/point.h>
#include <deal.II/grid/tria.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/vector.h>

using namespace dealii;

#include "occ_normal_projection.h"

namespace OpenCascade
{

  struct comp_points_on_curve
  {
    bool operator() (std::pair<double, double> a,std::pair<double, double> b) const
    {
      return a.second<b.second;
    }
  };



  class LineSmoothing
  {
  public:
    //! Smooth all dofs in
    // euler_vector for which
    // smoothing_dofs is set to
    // true. The euler_vector is
    // interpreted as a collection
    // of triples associated with
    // the displacement of the
    // nodes of the mesh, like the
    // one which is used in
    // MappingQEulerian.
    //
    // A static Q1 mapping is used
    // to compute the location of
    // the reference
    // points. Exceptions are
    // thrown if the points are not
    // located on the Geometry
    // curve. The parameter
    // base_point_id, identifies
    // the 3 dofs associated with
    // the base point, which is the
    // one left untouched. All
    // other points are moved
    // proportionally to their
    // arclength distance from the
    // base point. The driving
    // point id is given by
    // drivint_point_id. These two
    // ids are such that
    // 3*base_point_id+i is the ith
    // component of the base_point,
    // contained in
    // euler_vector(3*base_point_id+i).

    LineSmoothing(Vector<double> &euler_vector,
                  Handle(Geom_Curve) ref_curve,
                  TopLoc_Location *curr_loc,
                  const DoFHandler<2,3> &dh,
                  const std::vector<bool> &smoothing_dofs,
                  const unsigned int base_point_id,
                  const unsigned int driving_point_id,
                  const double tolerance = 1e-4);

    void update_reference(unsigned int base_point_id,
                          unsigned int driving_point_id);
    /** Perform the actual
    smoothing. Notice that the
    argument decides wether or
    not the moving point is
    projected to the original
    curve. If this is not the
    case, then all points are
    maintained in their current
    location, a new curve is
    computed which passes
    through all current points,
    and their smoothing is
    computed according to their
    original location. */
    void smooth(bool maintain_on_original_curve=true);

    inline Handle(Geom_Curve) get_curve()
    {
      return this->curve;
    }

    inline Vector<double> &get_lengths_before_smoothing()
    {
      return this->lengths_before_smoothing;
    }

    inline Vector<double> &get_lengths_after_smoothing()
    {
      return this->lengths_after_smoothing;
    }

    inline std::vector<unsigned int> &get_node_indices()
    {
      return this->node_indices;
    }


    void get_curve_tangent_vectors_at_smoothing_dofs(Vector<double> &tangents);

    void get_curve_length_ratios_at_smoothing_dofs(Vector<double> &length_ratios);


    typedef std::map<std::pair<double, double>, unsigned int, comp_points_on_curve >::iterator iterator;

  private:
    Vector<double> &euler_vector;
    Handle(Geom_Curve) ref_curve;
    TopLoc_Location ref_location;
    TopLoc_Location used_location;
    TopLoc_Location *curr_location;
    Handle(Geom_Curve) curve;
    NormalProjection<1> projection;

    const DoFHandler<2,3> &dh;
    const std::vector<bool> &smoothing_dofs;
    unsigned int base_point_id;
    unsigned int driving_point_id;
    double occ_base_t;
    double occ_driving_t;

    std::vector<Point<3> > support_points;

    double ref_L;
    double tolerance;

    Vector<double> fixed_length_ratios;
    Vector<double> lengths_before_smoothing;
    Vector<double> lengths_after_smoothing;
    std::vector<unsigned int> node_indices;
    std::map<std::pair<double, double>, unsigned int, comp_points_on_curve > smoothing_list;
    GeomLib_Tool tool;
  };


}

#endif
