#include "bem_utilities.h"

#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>

using namespace std;

namespace BEMUtilities 
{
  void remove_mesh_anisotropy(Triangulation<2,3> &tria)
  {
    Triangulation<2,3>::active_cell_iterator
      cell = tria.begin_active(), endc = tria.end();
    unsigned int refinedCellCounter = 1;
    while(refinedCellCounter)
      {
	refinedCellCounter = 0;
	for (cell=tria.begin_active(); cell!= endc;++cell)
	  {
	    if (cell->extent_in_direction(0) > 1.5*cell->extent_in_direction(1))
	      {
		cell->set_refine_flag(RefinementCase<2>::cut_axis(0));
		refinedCellCounter++;
	      }
	    else
	      {
		if (cell->extent_in_direction(1) > 1.5*cell->extent_in_direction(0))
		  {
		    cell->set_refine_flag(RefinementCase<2>::cut_axis(1));
		    refinedCellCounter++;
		  }
	      }
	  }
	tria.execute_coarsening_and_refinement();
      }
  }
}
