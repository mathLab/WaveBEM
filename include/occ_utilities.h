#ifndef occ_utilities_h
#define occ_utilities_h				 

#include <string>
#include <TopoDS_Shape.hxx>

#include <Geom_Plane.hxx>
#include <Geom_Curve.hxx>
#include <gp_Pnt.hxx>

#include <deal.II/base/point.h>


/** We collect in this namespace all utilities which operate on
    OpenCascade entities which don't need classes of their own.
*/

namespace OpenCascade 
{
    struct edge_ordering_rule
    {
    bool operator() (std::pair<unsigned int, double> a,std::pair<unsigned int, double> b) const
    {return a.second<b.second;}
    };

				   //! Read IGES files and translate
				   // their content into openCascade
				   // topological entities. The option
				   // scale_factor is used to
				   // compensate for different units
				   // being used in the IGES files and
				   // in the target application. 
  TopoDS_Shape read_IGES(std::string filename, double scale_factor=1e-3);
  
				   //! Perform the intersection of the
				   //  given topological shape with
				   //  the given plane. The returned
				   //  topological shape will contain
				   //  as few bsplines as possible.
  void intersect_plane(const TopoDS_Shape &in_shape,
                       TopoDS_Shape &out_shape,
		       const double c_x,
		       const double c_y,
		       const double c_z,
		       const double c,
		       const double tolerance=1e-7);

  TopoDS_Shape extract_xz_edges(const TopoDS_Shape &in_shape,
				const double tolerance=1e-7,
				const unsigned int max_num_edges=20);

				   // identifies transom edge, which is the
                                   // edge with the highest mean x coordinate
                                   // fails when transom is composed by more than
                                   // one face/edge, so this method needs a trim
  TopoDS_Shape extract_transom_edges(const TopoDS_Shape &in_shape,
                                     const unsigned int num_transom_edges,
				     const double tolerance=1e-7);

				   // creates a 3D smooth curve (bspline) passing
                                   // through the points in the
                                   // assigned vector. The points are
                                   // reordered internally according
                                   // to their scalar product with the
                                   // direction, if direction is
                                   // different from zero, otherwise
                                   // they are used as passed.
  Handle_Geom_Curve interpolation_curve_points_sort(std::vector<dealii::Point<3> >  &curve_points,
					            dealii::Point<3> direction=dealii::Point<3>());

  Handle_Geom_Curve  interpolation_curve(std::vector<dealii::Point<3> > &curve_points);
				   //! Convert OpenCascade point into
				   // Deal.II point.
  inline dealii::Point<3> Pnt(const gp_Pnt &p)
  {
    dealii::Point<3> P(p.X(), p.Y(), p.Z());
    return P;
  }

				   //! Convert Deal.II point into
				   // OpenCascade point.
  inline gp_Pnt Pnt(const dealii::Point<3> &p)
  {
    gp_Pnt P(p(0), p(1), p(2));
    return P;
  } 
  

  inline bool point_compare(const dealii::Point<3> &p1, const dealii::Point<3> &p2,
			    const dealii::Point<3> &direction) 
  {
    return (p1*direction < p2*direction);
  }

  void feature_edges_detection(); 
}

#endif
