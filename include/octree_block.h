#ifndef octree_block_h
#define octree_block_h

#include <base/smartpointer.h>
#include <base/convergence_table.h>
#include <base/quadrature_lib.h>
#include <base/quadrature_selector.h>
#include <base/parsed_function.h>
#include <base/utilities.h>
#include <base/point.h>

#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/matrix_lib.h>
#include <lac/vector.h>
#include <lac/solver_control.h>
#include <lac/solver_gmres.h>
#include <lac/precondition.h>

#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/grid_generator.h>
#include <grid/grid_in.h>
#include <grid/grid_out.h>
#include <grid/tria_boundary_lib.h>

#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <dofs/dof_renumbering.h>
#include <fe/fe_q.h>
#include <fe/fe_values.h>
#include <fe/fe_system.h>
#include <fe/mapping_q1_eulerian.h>
#include <fe/mapping_q1.h>

#include <numerics/data_out.h>
#include <numerics/vector_tools.h>
#include <numerics/solution_transfer.h>

				 // And here are a few C++ standard header
				 // files that we will need:
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <map>


using namespace dealii;

template <int dim>
class OctreeBlock
{

	public :
	
	        typedef typename DoFHandler<dim-1,dim>::active_cell_iterator cell_it; 
		
		
	private :	
		
		unsigned int level;

		unsigned int parentId;

		unsigned int numChildren;

		unsigned int childrenId[8];
		
		std::vector <std::set <unsigned int> > nearNeigh;
		
		std::vector <std::set <unsigned int> > intList;
		
		std::vector <std::set <unsigned int> > nonIntList;
		
		Point<dim> pMin;
		
		double delta;
		
		std::vector <unsigned int> nodesId;
		
		std::map <cell_it, std::vector<unsigned int> > quadPointsId;
		
	
public: 

        OctreeBlock();

	OctreeBlock(unsigned int level, unsigned int parent, Point<dim> pMin, double delta);
	
        OctreeBlock(const OctreeBlock<dim> &other);
	
	~OctreeBlock();
	
	void CopyContent(const  OctreeBlock *other);
	
        void AddNode(unsigned int nodeId);

        void AddQuadPoint(cell_it elemPointer, unsigned int quadPointId);
	
	std::vector <unsigned int> GetBlockNodeList() const;
	
	void DelNodeList();

        std::map <cell_it, std::vector<unsigned int> > GetBlockQuadPointsList() const;
	
	void DelQuadPointsList();

	unsigned int GetBlockNodesNum() const;
	
	unsigned int GetBlockChildrenNum() const;
	
	unsigned int GetParentId() const;
		
	void AddChild(unsigned int childId);
	
	unsigned int GetChildId(unsigned int idInList) const;
		
	Point<dim> GetPMin() const;
		
	double GetDelta() const;
		
	void AddNearNeigh(unsigned int sublevel, const unsigned int nnBlockId);
		
	unsigned int NumNearNeigh(unsigned int sublevel) const;
		
	unsigned int NumNearNeighLevels() const;
		
	std::set <unsigned int> GetNearNeighs(unsigned int sublevel) const;
		
	void AddBlockToIntList(unsigned int sublevel, const unsigned int intListBlockId);
		
	unsigned int NumIntList(unsigned int sublevel) const;
		
	unsigned int NumIntListLevels() const;
		
	std::set <unsigned int> GetIntList(unsigned int sublevel) const;

        std::vector<std::set <unsigned int> > GetIntList() const;
		
	void AddBlockToNonIntList(unsigned int sublevel, const unsigned int intListBlockId);
		
	unsigned int NumNonIntList(unsigned int sublevel) const;
		
	unsigned int NumNonIntListLevels() const;
		
	std::set <unsigned int> GetNonIntList(unsigned int sublevel) const;

	void SetNearNeighSize(unsigned int sublevels);

	void SetIntListSize(unsigned int sublevels);

	void SetNonIntListSize(unsigned int sublevels);
	
	unsigned int GetNearNeighSize() const;

	unsigned int GetIntListSize() const;

	unsigned int GetNonIntListSize() const;
	
};


#endif
