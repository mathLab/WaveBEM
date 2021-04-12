
#ifndef boat_model_h
#define boat_model_h

#include <fstream>
#include <iostream>
#include <cmath>
#include <limits>
#include <vector>

#include <deal.II/base/logstream.h>
#include <deal.II/opencascade/manifold_lib.h>
#include <deal.II/opencascade/utilities.h>

#include "occ_utilities.h"
#include "occ_normal_projection.h"
#include "occ_arclength_projection.h"
#include "occ_axis_projection.h"
#include "occ_line_smoothing.h"


class BoatModel
{
public:

  BoatModel();

  ~BoatModel();

  void start_iges_model(std::string igesFileName,
                        double scale,
                        double displacement,
                        double assigned_sink,
                        double assigned_trim,
                        double back_keel_length=0.1,
                        double front_keel_length=0.05,
                        double middle_keel_length=.47,
                        unsigned int number_of_transom_edges=1);

  TopoDS_Shape ReverseFaceOrientation(const TopoDS_Shape &shape,
                                      const TopoDS_Face &face);

  void compute_hydrostatic_sink(double &sink, const double &weight);

  Point<3> compute_hydrostatic_force(const double &sink);

  Point<3> compute_hydrostatic_moment(const double &sink, const Point<3> moment_center);

  gp_Trsf set_current_position(const Point<3> &translation_vect,
                               const double &quaternion_scalar,
                               const Point<3> &quaternion_vect);

//private:
  // keel intersection with undisturbed free surface at bow
  Point<3>  PointFrontTop;
  // this point will be "cornerpoint" of boat mesh at bow
  Point<3> PointFrontBot;
  // point on the undisturbed free surface and on boat surface
  // located roughly midway from bow to stern
  Point<3> PointMidTop;
  // point on the boat keel located roughly midway from bow to stern
  Point<3> PointMidBot;
  // keel intersection with undisturbed free surface at stern
  Point<3>  PointBackTop;
  // this point will be "cornerpoint" of boat mesh at stern
  Point<3> PointBackBot;
  // this point is intersection of left transom edge with undisturbed free surface
  Point<3> PointLeftTransom;
  // this point is intersection of right transom edge with undisturbed free surface
  Point<3> PointRightTransom;
  // this point is the base point of the transom edge (intersection with symmetry plane)
  Point<3> CurrentPointCenterTransom;
  // this point is intersection of left transom edge with undisturbed free surface
  // in current hull configuration
  Point<3> CurrentPointLeftTransom;
  // this point is intersection of right transom edge with undisturbed free surface
  // in current hull configuration
  Point<3> CurrentPointRightTransom;
  // this point is the base point of the transom edge (intersection with symmetry plane)
  // in current hull configuration
  Point<3> PointCenterTransom;
  // this is the wet length of the boat
  double boatWetLength;
  // this is the wet surface of the boat
  double boatWetSurface;
  // this is the whole geometric model of the boat (right side)
  TopoDS_Shape sh;
  // this is the whole geometric model of the boat (left side)
  TopoDS_Shape refl_sh;
  // this is the keel of the boat
  TopoDS_Shape keel_edge;
  // this is the right transom edge of the boat
  TopoDS_Shape right_transom_edge;
  // this is the left transom edge of the boat
  TopoDS_Shape left_transom_edge;
  // this is the undisturbed right water line
  TopoDS_Shape right_undist_water_line;
  // this is the undisturbed left water line
  TopoDS_Shape left_undist_water_line;
  // this is the undisturbed left water surface
  TopoDS_Shape undisturbed_water_surface_face;
  // location of the translated shapes (for hull roto-translated
  // in HYDROSTATIC position)
  TopLoc_Location reference_loc;
  // location of the translated curves (for hull roto-translated
  //in CURRENT position)
  TopLoc_Location current_loc;
  // surface normal projector on boat surface right side
  std::shared_ptr<dealii::OpenCASCADE::NormalProjectionManifold<2,3> > Boat_surface_right;
  // surface normal projector on boat surface left side
  std::shared_ptr<dealii::OpenCASCADE::NormalProjectionManifold<2,3> > Boat_surface_left;
  // y-axis direction projector on boat surface right side
  std::shared_ptr<dealii::OpenCASCADE::DirectionalProjectionManifold<2,3> > Boat_water_line_right;
  // y-axis direction projector on boat surface left side
  std::shared_ptr<dealii::OpenCASCADE::DirectionalProjectionManifold<2,3> > Boat_water_line_left;
  // z-axis direction projector on undisturbed free surface
  std::shared_ptr<dealii::OpenCASCADE::DirectionalProjectionManifold<2,3> > Undist_water_surf;
  // arclength projection on keel
  std::shared_ptr<dealii::OpenCASCADE::ArclengthProjectionLineManifold<2,3> > Boat_keel;
  // normal projection on keel
  std::shared_ptr<dealii::OpenCASCADE::NormalProjectionManifold<1,3> > Boat_keel_norm;
  // arclength projection on transom left edge
  std::shared_ptr<dealii::OpenCASCADE::ArclengthProjectionLineManifold<2,3> > Boat_transom_left;
  // arclength projection on transom right edge
  std::shared_ptr<dealii::OpenCASCADE::ArclengthProjectionLineManifold<2,3> > Boat_transom_right;
  // arclength projection on left wake line
  std::shared_ptr<dealii::OpenCASCADE::ArclengthProjectionLineManifold<2,3> > Wake_line_left;
  // arclength projection on right wake line
  std::shared_ptr<dealii::OpenCASCADE::ArclengthProjectionLineManifold<2,3> > Wake_line_right;
  // arclength projection on boat surface right side
  std::shared_ptr<dealii::OpenCASCADE::ArclengthProjectionLineManifold<2,3> > Water_line_right;
  // arclength projectionon boat surface left side
  std::shared_ptr<dealii::OpenCASCADE::ArclengthProjectionLineManifold<2,3> > Water_line_left;
  // keel curve
  Handle(Geom_Curve) equiv_keel_bspline;
  // undisturbed right water line curve
  Handle(Geom_Curve) right_undisturbed_waterline_curve;
  // undisturbed left water line curve
  Handle(Geom_Curve) left_undisturbed_waterline_curve;
  // left transom edge bspline
  Handle(Geom_Curve) left_transom_bspline;
  // right transom edge bspline
  Handle(Geom_Curve) right_transom_bspline;
  // right left bspline
  Handle(Geom_Curve) left_wake_bspline;
  // right wake bspline
  Handle(Geom_Curve) right_wake_bspline;
  // flag to determine if we have a transom stern or not
  bool is_transom;

  double boat_mass;

  Point<3> current_left_transom_tangent;

  Point<3> current_right_transom_tangent;
  // the position of the hull baricenter
  // in hydrostatic configuration
  Point<3> hydrostatic_hull_baricenter;
  // the position of the hull baricenter
  // in reference configuration
  Point<3> reference_hull_baricenter;
  // the position of the hull baricenter
  // in current configuration
  Point<3> current_hull_baricenter;
  // initial (assigned) sink
  double initial_sink;
  // initial (assigned) trim angle
  double initial_trim;
  // Euler angles (only for output purposes)
  double yaw_angle;
  double  pitch_angle;
  double  roll_angle;

};

#endif
