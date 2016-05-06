#ifndef bem_utilities_h
#define bem_utilities_h

#include <string>
#include <deal.II/base/point.h>
#include <deal.II/grid/tria.h>

using namespace dealii;

/** We collect in this namespace all general BEM utilities.
*/

namespace BEMUtilities
{
  void remove_mesh_anisotropy(Triangulation<2,3> &tria);
}

#endif
