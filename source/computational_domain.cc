//----------------------------  step-34.cc  ---------------------------
//    $Id: step-34.cc 18734 2009-04-25 13:36:48Z heltai $
//    Version: $Name$ 
//
//    Copyright (C) 2009, 2010, 2011 by the deal.II authors 
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//    Authors: Luca Heltai, Cataldo Manigrasso, Andrea Mola
//
//----------------------------  step-34.cc  ---------------------------

#define TOLL 0.001
#define MAXELEMENTSPERBLOCK 1

#include "../include/computational_domain.h"
#include <dofs/dof_renumbering.h>
#include <grid/grid_refinement.h>

				 // @sect4{ComputationalDomain::ComputationalDomain and
				 // ComputationalDomain::read_parameters}
				 // The constructor initializes the
				 // variuous object in much the same
				 // way as done in the finite element
				 // programs such as step-4 or
				 // step-6. The only new ingredient
				 // here is the ParsedFunction object,
				 // which needs, at construction time,
				 // the specification of the number of
				 // components.
				 //
				 // For the exact solution the number
				 // of vector components is one, and
				 // no action is required since one is
				 // the default value for a
				 // ParsedFunction object. The wind,
				 // however, requires dim components
				 // to be specified. Notice that when
				 // declaring entries in a parameter
				 // file for the expression of the
				 // Functions::ParsedFunction, we need
				 // to specify the number of
				 // components explicitly, since the
				 // function
				 // Functions::ParsedFunction::declare_parameters
				 // is static, and has no knowledge of
				 // the number of components.
template <int dim>
ComputationalDomain<dim>::ComputationalDomain(const unsigned int fe_degree,
	                                      const unsigned int mapping_degree)
		:
                mapping_degree(mapping_degree),
		tria(coarse_tria),
		fe(fe_degree),
		dh(tria),
                vector_fe(FE_Q<dim-1,dim>(fe_degree), dim),
                vector_dh(tria),
                mapping(NULL)
{}

template <int dim>
ComputationalDomain<dim>::~ComputationalDomain()
{
if (blocks.size() > 0)
   {
   for (unsigned int ii = 0; ii < num_blocks;  ii++)
        delete blocks[ii];
   }
 if(mapping != NULL)
   {
     delete mapping;
     mapping = NULL;
   }
}



template <int dim> 
void ComputationalDomain<dim>::declare_parameters (ParameterHandler &prm)
{    
  prm.enter_subsection("Quadrature rules");
  {
    prm.declare_entry("Quadrature type", "gauss", 
		      Patterns::Selection(QuadratureSelector<(dim-1)>::get_quadrature_names()));
    prm.declare_entry("Quadrature order", "4", Patterns::Integer());
    prm.declare_entry("Singular quadrature order", "5", Patterns::Integer());
  }
  prm.leave_subsection();
    

  
  prm.enter_subsection("Boundary Conditions ID Numbers");
  {    
   prm.declare_entry("Free Surface 1 ID", "1", Patterns::Integer());
   prm.declare_entry("Free Surface 2 ID", "0", Patterns::Integer());
   prm.declare_entry("Free Surface 3 ID", "0", Patterns::Integer());
   prm.declare_entry("Wall Surface 1 ID", "0", Patterns::Integer());
   prm.declare_entry("Wall Surface 2 ID", "0", Patterns::Integer());
   prm.declare_entry("Wall Surface 3 ID", "0", Patterns::Integer());
   prm.declare_entry("Inflow Surface 1 ID", "0", Patterns::Integer());
   prm.declare_entry("Inflow Surface 2 ID", "0", Patterns::Integer());
   prm.declare_entry("Inflow Surface 3 ID", "0", Patterns::Integer());
   prm.declare_entry("Free Surface Edge On Boat ID", "0", Patterns::Integer());
  }
  prm.leave_subsection();
  
  prm.enter_subsection("Octree Params");
  {    
   prm.declare_entry("Number of Octree Levels", "1", Patterns::Integer());
  }
  prm.leave_subsection();// to be moved

}

template <int dim> 
void ComputationalDomain<dim>::parse_parameters (ParameterHandler &prm)
{
   
  prm.enter_subsection("Quadrature rules");
  {
    quadrature =
      std_cxx1x::shared_ptr<Quadrature<dim-1> >
      (new QuadratureSelector<dim-1> (prm.get("Quadrature type"),
				      prm.get_integer("Quadrature order")));
    singular_quadrature_order = prm.get_integer("Singular quadrature order");
  }
  prm.leave_subsection();
    
  prm.enter_subsection("Boundary Conditions ID Numbers");
  {    
   free_sur_ID1 = prm.get_integer("Free Surface 1 ID");
   free_sur_ID2 = prm.get_integer("Free Surface 2 ID");
   free_sur_ID3 = prm.get_integer("Free Surface 3 ID");
   wall_sur_ID1 = prm.get_integer("Wall Surface 1 ID");
   wall_sur_ID2 = prm.get_integer("Wall Surface 2 ID");
   wall_sur_ID3 = prm.get_integer("Wall Surface 3 ID");
   inflow_sur_ID1 = prm.get_integer("Inflow Surface 1 ID");
   inflow_sur_ID2 = prm.get_integer("Inflow Surface 2 ID");
   inflow_sur_ID3 = prm.get_integer("Inflow Surface 3 ID");
   free_sur_edge_on_boat_ID = prm.get_integer("Free Surface Edge On Boat ID"); 
  }
  prm.leave_subsection();
  
  prm.enter_subsection("Octree Params");
  {
   num_octree_levels = prm.get_integer("Number of Octree Levels");
  }
  prm.leave_subsection();

}


				 // @sect4{ComputationalDomain::read_domain}
    
				 // A boundary element method
				 // triangulation is basically the
				 // same as a (dim-1) dimensional
				 // triangulation, with the difference
				 // that the vertices belong to a
				 // (dim) dimensional space.
				 //
				 // Some of the mesh formats supported
				 // in deal.II use by default three
				 // dimensional points to describe
				 // meshes. These are the formats
				 // which are compatible with the
				 // boundary element method
				 // capabilities of deal.II. In
				 // particular we can use either UCD
				 // or GMSH formats. In both cases, we
				 // have to be particularly careful
				 // with the orientation of the mesh,
				 // because, unlike in the standard
				 // finite element case, no reordering
				 // or compatibility check is
				 // performed here.  All meshes are
				 // considered as oriented, because
				 // they are embedded in a higher
				 // dimensional space. (See the
				 // documentation of the GridIn and of
				 // the Triangulation for further
				 // details on orientation of cells in
				 // a triangulation.) In our case, the
				 // normals to the mesh are external
				 // to both the circle in 2d or the
				 // sphere in 3d.
				 //
				 // The other detail that is required
				 // for appropriate refinement of the
				 // boundary element mesh, is an
				 // accurate description of the
				 // manifold that the mesh is
				 // approximating. We already saw this
				 // several times for the boundary of
				 // standard finite element meshes
				 // (for example in step-5 and
				 // step-6), and here the principle
				 // and usage is the same, except that
				 // the HyperBallBoundary class takes
				 // an additional template parameter
				 // that specifies the embedding space
				 // dimension. The function object
				 // still has to be static to live at
				 // least as long as the triangulation
				 // object to which it is attached.
        
template <int dim>
void ComputationalDomain<dim>::read_domain()
{
  tria.set_mesh_smoothing (Triangulation<dim-1,dim>::do_not_produce_unrefined_islands );
  coarse_tria.set_mesh_smoothing (Triangulation<dim-1,dim>::do_not_produce_unrefined_islands );
  tria.set_mesh_smoothing (Triangulation<dim-1,dim>::eliminate_refined_boundary_islands );
  coarse_tria.set_mesh_smoothing (Triangulation<dim-1,dim>::eliminate_refined_boundary_islands );

  tria.restore();

}


				 // @sect4{ComputationalDomain::refine_and_resize}

				 // This function globally refines the
				 // mesh, distributes degrees of
				 // freedom, and resizes matrices and
				 // vectors.

template <int dim>
void ComputationalDomain<dim>::refine_and_resize()
{}    

template <int dim>
void ComputationalDomain<dim>::compute_min_diameter()
{
  min_diameter = 10000;
  typename Triangulation<dim-1,dim>::active_cell_iterator
    cell = tria.begin_active(), endc = tria.end();

  for( ; cell != endc; ++cell)
    {
      min_diameter = std::min(min_diameter,cell->diameter());
    }  
  std::cout << "Min diameter: << " << min_diameter << std::endl;
}


template <int dim>
bool ComputationalDomain<dim>::mesh_check_and_update()
{
  return false;
}

template <int dim>
void ComputationalDomain<dim>::generate_double_nodes_set()
{ 
  std::cout<<"Generating double nodes set..."<<std::endl;

				 // The following is the function
				 // which creates a set containing
				 // the double nodes.
				 

  std::vector<bool> boundary_dofs(vector_dh.n_dofs(),false);
  std::vector< bool > comp_sel(dim,true);   
  DoFTools::extract_boundary_dofs(vector_dh,comp_sel,boundary_dofs);

  double tol = 1e-8;
  double_nodes_set.clear();
  double_nodes_set.resize(dh.n_dofs());
  vector_double_nodes_set.clear();
  vector_double_nodes_set.resize(vector_dh.n_dofs());

  update_support_points();
  
  for (unsigned int i=0; i<dh.n_dofs(); ++i)  
      {
      double_nodes_set[i].insert(i);
      for (unsigned int k=0; k<dim; ++k)
          vector_double_nodes_set[dim*i+k].insert(dim*i+k);
      if (boundary_dofs[dim*i] == true)
         {
         for (unsigned int j=0; j<dh.n_dofs(); ++j) 
	     {
	  //std::cout<<"i "<<i<<" ("<<support_points[i]<<")  j "<<j<<" ("<<support_points[j]<<")  distance "<<support_points[i].distance(support_points[j])<<std::endl;
	     if (support_points[i].distance(support_points[j]) < tol)
	        {
             	//std::cout<<"i "<<i<<" ("<<support_points[i]<<")  j "<<j<<" ("<<support_points[j]<<")  distance "<<support_points[i].distance(support_points[j])<<std::endl;
	        double_nodes_set[i].insert(j);
	        //std::cout<<"i "<<i<<" double "<<j<<" "<<tol<<std::endl;
                for (unsigned int k=0; k<dim; ++k)
                    {
                    vector_double_nodes_set[dim*i+k].insert(dim*j+k);
                    //std::cout<<"dim*i+k "<<dim*i+k<<" ("<<vector_support_points[dim*i+k]<<")   dim*j+k "<<dim*j+k<<" ("<<vector_support_points[dim*j+k]<<")  distance "<<vector_support_points[dim*i+k].distance(vector_support_points[dim*j+k])<<std::endl;
                    }
	        }
	     }
         }
	
      }
      

std::cout<<"...done"<<std::endl;
}



template <int dim>
void ComputationalDomain<dim>::generate_octree_blocking()
{

std::cout<<"Generating octree blocking... "<<std::endl;

				 // @sect5{BEMProblem::generate_double_nodes_set}

				 // The following is the function
				 // which creates the octree blocking
				 // for the fast multipole algorithm


std::vector<Point<dim> > support_points(dh.n_dofs());
DoFTools::map_dofs_to_support_points<dim-1, dim>( *mapping,
                                                  dh, support_points);

FEValues<dim-1,dim> fe_v(*mapping,fe, *quadrature,
			 update_values |
			 update_cell_normal_vectors |
			 update_quadrature_points |
			 update_JxW_values);

double max_coor_value = 0;

for (unsigned int i=0; i < dh.n_dofs(); i++)
    {
    //for printout
    //std::cout<<"Node "<<i<< "["<<support_points[i]<<"] "<<std::endl;
    for (unsigned int j=0; j < dim; j++)
        {
	max_coor_value = std::max(max_coor_value,std::abs(support_points[i](j)));
	}
    }
    
if (blocks.size() > 0)
   {
   for (unsigned int ii = 0; ii < num_blocks;  ii++)
        delete blocks[ii];
   }
    
unsigned int maxNumBlocks = num_octree_levels*tria.n_active_cells()*fe_v.n_quadrature_points;
//unsigned int maxNumBlocks = 0;
//for (unsigned int ii = 0; ii < num_octree_levels + 1;  ii++)
//	{
//	 maxNumBlocks += int(pow(8.,double(ii)));
//	}

blocks.clear();
blocks.reserve(maxNumBlocks);
blocks.resize(maxNumBlocks);

unsigned int blocksCount = 0;
startLevel.resize(num_octree_levels+1);
endLevel.resize(num_octree_levels+1);

//for (unsigned int j=0; j < num_octree_levels + 1; j++)
//     parentList[j].clear();
parentList.clear();
parentList.resize(num_octree_levels+1);
parentList[0].push_back(0);


childlessList.clear();
unsigned int numChildless = 0;
numParent.resize(num_octree_levels+1);

//qui di seguito vengono reinizializzate strutture utili al multipolo

// mappa che associa ad ogni dof un vettore con i blocchi cui essa appartiene per ogni livello
dof_to_block.clear();

// mappa che associa ad ogni quad point un vettore con i blocchi cui essa appartiene per ogni livello
quad_point_to_block.clear();

// vettore di vettori contenente per ogni livello, gli ids dei blocchi
// contenenti almeno un dof
dofs_filled_blocks.clear();

// vettore di vettori contenente per ogni livello, gli ids dei blocchi
// contenenti almeno un quad point
quad_points_filled_blocks.clear();

quadPoints.clear();
quadNormals.clear();
quadShapeFunValues.clear();
quadJxW.clear();

dofs_filled_blocks.resize(num_octree_levels+1);

quad_points_filled_blocks.resize(num_octree_levels+1);



for (unsigned int ii = 0; ii < num_octree_levels + 1 ;  ii++)
	{
	 numParent[ii] = 0;
	}



Point<dim> pMin;
for (int i=0; i<dim; i++)
   pMin(i) = -1.1*max_coor_value;

				 // delta e' il lato del kazzo di kubo...
double delta = 2.2*max_coor_value;

OctreeBlock<dim>* block = new OctreeBlock<dim>(0, 0, pMin, delta);

std::vector<unsigned int> local_dof_indices(fe.dofs_per_cell);
cell_it
cell = dh.begin_active(),
endc = dh.end();
for (cell = dh.begin_active(); cell != endc; ++cell)    
    {
    fe_v.reinit(cell);
    const unsigned int n_q_points = fe_v.n_quadrature_points;
    quadPoints[cell] = fe_v.get_quadrature_points();
    quadNormals[cell] = fe_v.get_normal_vectors();
    quadJxW[cell].resize(n_q_points);
    quadShapeFunValues[cell].resize(n_q_points);
    for(unsigned int q=0; q<n_q_points; ++q)
       {
       quadJxW[cell][q] = fe_v.JxW(q);
       for(unsigned int j=0; j<fe.dofs_per_cell; ++j)
          quadShapeFunValues[cell][q].push_back(fe_v.shape_value(j,q));
       }
    
    quad_point_to_block[cell].resize(n_q_points);		     
    for (unsigned int j=0; j<n_q_points; ++j)
        {
	block->AddQuadPoint(cell,j);
	quad_point_to_block[cell][j].push_back(0);
	}
    
    cell->get_dof_indices(local_dof_indices);
    for(unsigned int j=0; j<fe.dofs_per_cell; ++j)
        {
	dof_to_elems[local_dof_indices[j]].push_back(cell);
	}
    }
    
for (unsigned int ii = 0; ii < dh.n_dofs(); ii++)
    {
    block->AddNode(ii);
    dof_to_block[ii].push_back(0);
    }    


// just for output    
/*for (cell = dh.begin_active(); cell != endc; ++cell)    
    {
    std::set<cell_it> surr_elems = elem_to_surr_elems[cell];
    std::cout<<std::endl<<"cell "<<cell<<"  surrounded by: "; 
    for (typename std::set<cell_it>::iterator pos = surr_elems.begin(); pos !=surr_elems.end(); pos++)
         std::cout<<" "<<*pos;
    }*/
    
blocks[0] = block;
numParent[0] = 1;

//std::cout<<"blocks[0].GetBlockChildrenNum() "<<blocks[0].GetBlockChildrenNum()<<std::endl;

/*std::cout<<std::endl;
std::cout<<blocks[0].GetPMin()<<std::endl;
std::cout<<pMin<<std::endl;
std::cout<<block.GetDelta()<<std::endl;
std::cout<<block.GetBlockNodeList()[0]<<std::endl;
std::cout<<block.GetBlockElementsList()[1]<<std::endl;
std::cout<<delta<<std::endl;
std::cout<<std::endl;//*/

unsigned int quadPointsInChildless = 0;
unsigned int nodesInChildless = 0;

for (unsigned int level = 1; level < num_octree_levels + 1;  level++)

    {
    unsigned int quadPointsCheck = quadPointsInChildless;
    unsigned int nodesCheck = nodesInChildless;
    delta /= 2.;

    for (unsigned int kk = 0; kk < numParent[level-1];  kk++)

        {
	unsigned int jj = parentList[level-1][kk];
        //std::cout<<" level "<<level<<"     block "<<jj<<std::endl;
	OctreeBlock<dim> *parent = blocks[jj];
        //std::cout<<"parent.GetBlockChildrenNum() "<<parent.GetBlockChildrenNum()<<std::endl;
        //std::cout<<" Pmin "<<parent.GetPMin()(0)<<", "<<parent.GetPMin()(1)<<", "<<parent.GetPMin()(2)<<" "<<std::endl;
        //std::cout<<" delta "<<parent.GetDelta()<<" "<<std::endl;

        pMin = parent->GetPMin();
	unsigned int num_children_per_block = int(pow((double)2,(double)dim));
        std::vector<OctreeBlock<dim> *> children(num_children_per_block);
	
	if (dim == 3)
	   {
           children[0] = new OctreeBlock<dim>(level, jj, pMin                              , delta);
           children[1] = new OctreeBlock<dim>(level, jj, pMin+Point<dim>(delta   ,0.,   0.), delta);
           children[2] = new OctreeBlock<dim>(level, jj, pMin+Point<dim>(delta,delta,   0.), delta);
           children[3] = new OctreeBlock<dim>(level, jj, pMin+Point<dim>(0.   ,delta,   0.), delta);
           children[4] = new OctreeBlock<dim>(level, jj, pMin+Point<dim>(0.   ,   0.,delta), delta);
           children[5] = new OctreeBlock<dim>(level, jj, pMin+Point<dim>(delta,   0.,delta), delta);
           children[6] = new OctreeBlock<dim>(level, jj, pMin+Point<dim>(delta,delta,delta), delta);
           children[7] = new OctreeBlock<dim>(level, jj, pMin+Point<dim>(   0.,delta,delta), delta);
           }

	if (dim == 2)
	   {
           children[0] = new OctreeBlock<dim>(level, jj, pMin                        , delta);
           children[1] = new OctreeBlock<dim>(level, jj, pMin+Point<dim>(delta   ,0.), delta);
           children[2] = new OctreeBlock<dim>(level, jj, pMin+Point<dim>(delta,delta), delta);
           children[3] = new OctreeBlock<dim>(level, jj, pMin+Point<dim>(0.   ,delta), delta);
           
           }

           std::map <cell_it, std::vector <unsigned int> > blockQuadPointsList =
	                                   parent->GetBlockQuadPointsList();
					   
           std::vector <unsigned int> blockNodeList = parent->GetBlockNodeList();
        
	if (dim == 3) 
	{
        for (unsigned int i = 0; i < blockNodeList.size(); i++)    
            {
	    Point <dim> node = support_points[blockNodeList[i]];
	    // assegnamento nodi del blocco padre ai blocchi figli
	    
	    if (node(2) <= parent->GetPMin()(2)+delta)
	       {
	       if (node(1) <= parent->GetPMin()(1)+delta)
		  {
		  if (node(0) <= parent->GetPMin()(0)+delta)
		     {
		     //std::cout<<" Sono in 1 "<<std::endl;
		     children[0]->AddNode(blockNodeList[i]);
		     }
		  else
		     {
		     //std::cout<<" Sono in 2 "<<std::endl;
		     children[1]->AddNode(blockNodeList[i]);
		     }
		  }
	       else
		  {
                  if (node(0) <= parent->GetPMin()(0)+delta)
		     {
		     //std::cout<<" Sono in 4 "<<std::endl;
		     children[3]->AddNode(blockNodeList[i]);
		     }
		  else
		     {
		     //std::cout<<" Sono in 3 "<<std::endl;
		     children[2]->AddNode(blockNodeList[i]);
		     }
		  }
	       }
            else
	       {
	       if (node(1) <= parent->GetPMin()(1)+delta)
		  {
	          if (node(0) <= parent->GetPMin()(0)+delta)
		     {
		     //std::cout<<" Sono in 5 "<<std::endl;
		     children[4]->AddNode(blockNodeList[i]);
		     }
		  else
		     {
		     //std::cout<<" Sono in 6 "<<std::endl;
		     children[5]->AddNode(blockNodeList[i]);
		     }
		  }
               else
		  {
		  if (node(0) <= parent->GetPMin()(0)+delta)
		     {
		     //std::cout<<" Sono in 8 "<<std::endl;
		     children[7]->AddNode(blockNodeList[i]);
		     }
	          else
	             {
		     //std::cout<<" Sono in 7 "<<std::endl;
		     children[6]->AddNode(blockNodeList[i]);
		     }
		  }
	       } //fine assegnazione nodi del padre ai blocchi figli
	       
	    } //fine loop nodi del blocco 
	
	typename std::map <cell_it, std::vector<unsigned int> >::iterator it;    
        for (it = blockQuadPointsList.begin(); it != blockQuadPointsList.end(); it++)    
            {
	    for (unsigned int pp = 0; pp < (*it).second.size(); pp++)
	        {
	        Point<dim> quadPoint = quadPoints[(*it).first][(*it).second[pp]];
	        // assegnamento punti quadratura del blocco padre ai blocchi figli
	        if (quadPoint(2) <= parent->GetPMin()(2)+delta)
	           {
	           if (quadPoint(1) <= parent->GetPMin()(1)+delta)
		      {
		      if (quadPoint(0) <= parent->GetPMin()(0)+delta)
		         {
		         //std::cout<<" Sono in 1 "<<std::endl;
		         children[0]->AddQuadPoint((*it).first,(*it).second[pp]);
		         }
		      else
		         {
		         //std::cout<<" Sono in 2 "<<std::endl;
		         children[1]->AddQuadPoint((*it).first,(*it).second[pp]);
		         }
		      }
	           else
		      {
                      if (quadPoint(0) <= parent->GetPMin()(0)+delta)
		         {
		         //std::cout<<" Sono in 4 "<<std::endl;
		         children[3]->AddQuadPoint((*it).first,(*it).second[pp]);
		         }
		      else
		         {
		         //std::cout<<" Sono in 3 "<<std::endl;
		         children[2]->AddQuadPoint((*it).first,(*it).second[pp]);
		         }
		      }
	           }
                else
	           {
	           if (quadPoint(1) <= parent->GetPMin()(1)+delta)
		      {
	              if (quadPoint(0) <= parent->GetPMin()(0)+delta)
		         {
		         //std::cout<<" Sono in 5 "<<std::endl;
		         children[4]->AddQuadPoint((*it).first,(*it).second[pp]);
		         }
		      else
		         {
		         //std::cout<<" Sono in 6 "<<std::endl;
		         children[5]->AddQuadPoint((*it).first,(*it).second[pp]);
		         }
		      }
                   else
		      {
		      if (quadPoint(0) <= parent->GetPMin()(0)+delta)
		         {
		         //std::cout<<" Sono in 8 "<<std::endl;
		         children[7]->AddQuadPoint((*it).first,(*it).second[pp]);
		         }
	              else
	                 {
		         //std::cout<<" Sono in 7 "<<std::endl;
		         children[6]->AddQuadPoint((*it).first,(*it).second[pp]);
		         }
		      }
	           } //fine assegnazione punti quadratura del padre ai blocchi figli
	        }
	    } //fine loop punti quadratura del blocco

        for (unsigned int j=0; j < num_children_per_block; j++ )
	    {   
	    if (children[j]->GetBlockNodeList().size() + 
	        children[j]->GetBlockQuadPointsList().size()> 0)
	       {
	       blocksCount += 1;
	       blocks[blocksCount] = new OctreeBlock<dim>();
	       blocks[blocksCount]->CopyContent(children[j]);
	       delete children[j];
	       
	       parent->AddChild(blocksCount);
               std::map <cell_it, std::vector<unsigned int> >
	       blockQuadPointsList = blocks[blocksCount]->GetBlockQuadPointsList();
	       typename std::map <cell_it, std::vector<unsigned int> >::iterator it;    
               for (it = blockQuadPointsList.begin(); it != blockQuadPointsList.end(); it++)
	           {
		   cell_it cell = (*it).first;
		   for(unsigned int kk = 0; kk < (*it).second.size(); kk++)
		      {
		      quad_point_to_block[(*it).first][(*it).second[kk]].push_back(blocksCount);
		      }
		   }
               std::vector<unsigned int> blockNodesList = blocks[jj]->GetBlockNodeList();
	       for(unsigned int k = 0; k < blockNodesList.size(); k++)
	          dof_to_block[blockNodesList[k]].push_back(jj);
		   
	       }
	    else
	       {
	       delete children[j];
	       }

	    } // fine loop sui blocchi figlio appena creati
	
	} //fine ramo dim = 3 dell'if
	else
	{
	for (unsigned int i = 0; i < blockNodeList.size(); i++)    
            {
	    
	    // assegnamento nodi del blocco padre ai blocchi figli
            Point <dim> node = support_points[blockNodeList[i]];
	    
	       if (node(1) <= parent->GetPMin()(1)+delta)
		  {
		  if (node(0) <= parent->GetPMin()(0)+delta)
		     {
		     //std::cout<<" Sono in 1 "<<std::endl;
		     children[0]->AddNode(blockNodeList[i]);
		     }
		  else
		     {
		     //std::cout<<" Sono in 2 "<<std::endl;
		     children[1]->AddNode(blockNodeList[i]);
		     }
		  }
	       else
		  {
                  if (node(0) <= parent->GetPMin()(0)+delta)
		     {
		     //std::cout<<" Sono in 4 "<<std::endl;
		     children[3]->AddNode(blockNodeList[i]);
		     }
		  else
		     {
		     //std::cout<<" Sono in 3 "<<std::endl;
		     children[2]->AddNode(blockNodeList[i]);
		     }
		  }//fine assegnazione blocchi del padre ai blocchi figli
	       
	    } //fine loop nodi del blocco	
	
	typename std::map <cell_it, std::vector<unsigned int> >::iterator it;    
        for (it = blockQuadPointsList.begin(); it != blockQuadPointsList.end(); it++)    
            {
	    for (unsigned int pp = 0; pp < (*it).second.size(); pp++)
	        {
                // assegnamento quad points del blocco padre ai blocchi figli
               Point<dim> quadPoint = quadPoints[(*it).first][(*it).second[pp]];
	       if (quadPoint(1) <= parent->GetPMin()(1)+delta)
		  {
		  if (quadPoint(0) <= parent->GetPMin()(0)+delta)
		     {
		     //std::cout<<" Sono in 1 "<<std::endl;
		     children[0]->AddQuadPoint((*it).first,(*it).second[pp]);
		     }
		  else
		     {
		     //std::cout<<" Sono in 2 "<<std::endl;
		     children[1]->AddQuadPoint((*it).first,(*it).second[pp]);
		     }
		  }
	       else
		  {
                  if (quadPoint(0) <= parent->GetPMin()(0)+delta)
		     {
		     //std::cout<<" Sono in 4 "<<std::endl;
		     children[3]->AddQuadPoint((*it).first,(*it).second[pp]);
		     }
		  else
		     {
		     //std::cout<<" Sono in 3 "<<std::endl;
		     children[2]->AddQuadPoint((*it).first,(*it).second[pp]);
		     }
		  }//fine assegnazione blocchi del padre ai blocchi figli

		}    
  	    }
	    
        for (unsigned int j=0; j < num_children_per_block; j++ )
	    {   
	    if (children[j]->GetBlockNodeList().size() + 
	        children[j]->GetBlockQuadPointsList().size()> 0)
	       {
	       blocksCount += 1;
	       blocks[blocksCount] = new OctreeBlock<dim>();
	       blocks[blocksCount]->CopyContent(children[j]);
	       delete children[j];
	       
	       parent->AddChild(blocksCount);
               std::map <cell_it, std::vector<unsigned int> >
	       blockQuadPointsList = blocks[blocksCount]->GetBlockQuadPointsList();
	       typename std::map <cell_it, std::vector<unsigned int> >::iterator it;    
               for (it = blockQuadPointsList.begin(); it != blockQuadPointsList.end(); it++)
	           {
		   cell_it cell = (*it).first;
		   for(unsigned int kk = 0; kk < (*it).second.size(); kk++)
		      {
		      quad_point_to_block[(*it).first][(*it).second[kk]].push_back(blocksCount);
		      }
		   }
               std::vector<unsigned int> blockNodesList = blocks[jj]->GetBlockNodeList();
	       for(unsigned int k = 0; k < blockNodesList.size(); k++)
	          dof_to_block[blockNodesList[k]].push_back(jj);
		   
	       }
	    else
	       {
	       delete children[j];
	       }

	    } // fine loop sui blocchi figlio appena creati
	    
	
	} // fine ramo dim == 2 dell'if
	
	} //fine loop blocchi livello precedente	
    
	
    //double elemCheck = numChildless;

    startLevel[level] = endLevel[level-1] + 1;
    endLevel[level] = blocksCount;
	
	// here we loop over the blocks of the newly created level and
	// we decide if each block is to be split again in the next level:
	// if it contains more
	// than a node or quad point, it will be placed in the parent list.
	// Instead, if it only contains a node or quad point, it goes in the
	// childless list, and not be refined any more. It is important to
	// account for the presence of double nodes: if not, blocks will be
	// always refined    
    for (unsigned int jj = startLevel[level]; jj < endLevel[level]+1;  jj++)
	{

	// here we get the number of nodes in the block 
	std::vector<unsigned int> nodesId = blocks[jj]->GetBlockNodeList();
	int blockNumNodes = (int) nodesId.size();
	
	// now we compute the number of the nodes that are double of others
	int numDoubleNodes = 0;  
	for (unsigned int kk = 0; kk < nodesId.size();  kk++)
	    {
	    int a = (int) double_nodes_set[nodesId[kk]].size();
	    numDoubleNodes += a - 1;
	    }   
	
	// here we compute the number of quad points in the block
	int blockNumQuadPoints = 0;
	std::map <cell_it, std::vector<unsigned int> >
	blockQuadPointsList = blocks[jj]->GetBlockQuadPointsList();
	typename std::map <cell_it, std::vector<unsigned int> >::iterator it;    
        for (it = blockQuadPointsList.begin(); it != blockQuadPointsList.end(); it++)
	     blockNumQuadPoints += (int) (*it).second.size();
	//std::cout<<"Level "<<level<<" Block "<<jj<<"  nodes "<<blockNumNodes<<" + quadPoints "<<blockNumQuadPoints<<std::endl;     
	
	quadPointsCheck += blockNumQuadPoints;
	nodesCheck += blockNumNodes;
	// here we decide if a block is to be placed in the parent
	// or childless list
	//if (blockNumNodes + blockNumQuadPoints - numDoubleNodes < 2)
	if (blockNumNodes - numDoubleNodes < 2)
	   {
	   numChildless += 1;
	   childlessList.push_back(jj);
	   quadPointsInChildless += blockNumQuadPoints; 
	   nodesInChildless += blockNumNodes;
	   
	   
	   // if a block is childless, we must assign now the nodes and quad points
	   // that belong to it for all the next levels
	    
	for (unsigned int kk = 0; kk < nodesId.size();  kk++)   
	    for (unsigned int j = level+1; j < num_octree_levels+1; j++)
	        dof_to_block[nodesId[kk]].push_back(jj);
	    
	for (it = blockQuadPointsList.begin(); it != blockQuadPointsList.end(); it++)
	    for (unsigned int i = 0; i < (*it).second.size(); i++)    
	        for (unsigned int j = level+1; j < num_octree_levels+1; j++)
	             quad_point_to_block[(*it).first][(*it).second[i]].push_back(jj);      
	   
	   }
	else
	   {
	   numParent[level] += 1;
	   parentList[level].push_back(jj);
	   }
	
	   // let's update the list of node filled block
	if (blockNumNodes > 0)
	   dofs_filled_blocks[level].push_back(jj);

	   // let's update the list of quad point filled block	   
	if (blockNumQuadPoints > 0)
	   quad_points_filled_blocks[level].push_back(jj);   
	//elemCheck += blockNumNodes + blockNumQuadPoints;  
        }


     std::cout<<" Total nodes at level "<<level<<" of "<<num_octree_levels<<" are "<<nodesCheck<<std::endl;
     std::cout<<" Total quad points at level "<<level<<" of "<<num_octree_levels<<" are "<<quadPointsCheck<<std::endl;
     std::cout<<" Blocks at level "<<level<<" of "<<num_octree_levels<<" are "<<endLevel[level]-endLevel[level-1]<<std::endl;
     std::cout<<" Total blocks at level "<<level<<" of "<<num_octree_levels<<" are "<<endLevel[level] + 1<<std::endl;
     std::cout<<std::endl;	   
	
     } //fine loop livelli

childlessList.resize(childlessList.size()+parentList[num_octree_levels].size());

for (unsigned int jj = 0; jj < parentList[num_octree_levels].size();  jj++)
	{
	 childlessList[numChildless + jj] = parentList[num_octree_levels][jj];
	}




num_blocks = blocksCount+1;

	
std::cout<<"...done generating octree blocking"<<std::endl;

std::cout<<"Computing proximity lists for blocks"<<std::endl;

//just for output
/*for (cell = dh.begin_active(); cell != endc; ++cell)    
    {
    unsigned int levelCheck = elem_to_blocks[cell].size();
    std::cout<<std::endl<<"Elem "<<cell<<" is in the "<<levelCheck<<" blocks: ";
    for (unsigned int zz = 0; zz < levelCheck; zz++)
         std::cout<<elem_to_blocks[cell][zz]<<" ";
    }*/

/*for (cell_it cell = dh.begin_active(); cell != endc; ++cell)
    for (unsigned int j=0; j < quadPoints[cell].size(); j++)
        std::cout<<"Cell "<<cell<<"  QP "<<j<<"  of "<<quadPoints[cell].size()<<": "<<quadPoints[cell][j]<<std::endl;//*/

/*for (cell_it cell = dh.begin_active(); cell != endc; ++cell)
    for (unsigned int j=0; j < quad_point_to_block[cell].size(); j++)
        {
	std::cout<<"Cell "<<cell<<"  QP "<<j<<"  of "<<quad_point_to_block[cell].size()<<": ";
        for (unsigned int i=0; i < quad_point_to_block[cell][j].size(); i++)
            std::cout<<quad_point_to_block[cell][j][i]<<" ";
	std::cout<<std::endl;    
	}  //*/  


/*for (unsigned int i=0; i < dh.n_dofs(); i++)
    {
    std::cout<<"Node "<<i<<"  doubles: ";
    std::set <unsigned int> doubleNodes = double_nodes_set[i]; 
    for (std::set<unsigned int>::iterator pos = doubleNodes.begin(); pos != doubleNodes.end(); pos++)
        {
	std::cout<<*pos<<"( ";
	for (unsigned int j=0; j < dof_to_elems[*pos].size(); j++)
	    std::cout<<" "<<dof_to_elems[*pos][j];
	std::cout<<") ";
        }
    std::cout<<std::endl;
    } //*/



// ricerca blocchi nearest neighbors

for (unsigned int ii = startLevel[1]; ii < endLevel[1] + 1;  ii++)
	{
	for (unsigned int jj = startLevel[1]; jj < endLevel[1] + 1;  jj++)
		{
		blocks[ii]->AddNearNeigh(0,jj);
		}
	}


for (unsigned int level = 2; level < num_octree_levels + 1;  level++)

	{
	for (unsigned int kk = startLevel[level]; kk < endLevel[level]+1;  kk++)

		{
		OctreeBlock<dim> *block1 = blocks[kk];
		block1->AddNearNeigh(0,kk); // a block is NearNeigh of itself

		double delta1 = block1->GetDelta();
		Point<dim> PMin1 = block1->GetPMin();
		Point<dim> Center1;
		for (unsigned int iii = 0; iii < dim; iii++)
		     Center1(iii) = delta1/2.;
		Point<dim> PMax1 = 2*Center1;     
		PMax1 += PMin1;
		Center1 += PMin1;
		unsigned int parentId = block1->GetParentId();
		std::set <unsigned int> parentNNeighs = blocks[parentId]->GetNearNeighs(0);

                // the nearest neighbors are searched among the father's nearest neighbors children
		for (std::set <unsigned int>::iterator pos = parentNNeighs.begin(); pos != parentNNeighs.end();  pos++)

			{
		        if (blocks[*pos]->GetBlockChildrenNum() == 0) // if a parent's near neigh is childless, he can be a near neigh: let's check
				{
				unsigned int block2Id = *pos;
				OctreeBlock<dim> *block2 = blocks[block2Id];
				double delta2 = block2->GetDelta();
		                Point<dim> PMin2 = block2->GetPMin();
		                Point<dim> Center2;
		                for (unsigned int iii = 0; iii < dim; iii++)
		                    Center2(iii) = delta2/2.;
		                Point<dim> PMax2 = 2*Center2;     
		                PMax2 += PMin2;
		                Center2 += PMin2;
				
				if (dim == 3)
				{
				if 	((fabs(PMin1(0) - PMax2(0)) <= TOLL) || (fabs(PMax1(0) - PMin2(0)) <= TOLL))
					{
					if ((PMin1(1)-TOLL <= PMax2(1)) && (PMax1(1)+TOLL >= PMin2(1)))
						{
						if ((PMin1(2)-TOLL <= PMax2(2)) && (PMax1(2)+TOLL >= PMin2(2)))
							{
							block1->AddNearNeigh(0,block2Id);
							//std::cout<<" *"<<block2Id;
							}
						}
					}

				if ((fabs(PMin1(1) - PMax2(1)) <= TOLL) || (fabs(PMax1(1) - PMin2(1)) <= TOLL))
					{
					if ((PMin1(0)-TOLL <= PMax2(0)) && (PMax1(0)+TOLL >= PMin2(0)))
						{
						if ((PMin1(2)-TOLL <= PMax2(2)) && (PMax1(2)+TOLL >= PMin2(2)))
							{
							block1->AddNearNeigh(0,block2Id);
							//std::cout<<" *"<<block2Id;
							}
						}
					}

				if ((fabs(PMin1(2) - PMax2(2)) <= TOLL) || (fabs(PMax1(2) - PMin2(2)) <= TOLL))
					{
					if ((PMin1(1)-TOLL <= PMax2(1)) && (PMax1(1)+TOLL >= PMin2(1)))
						{
						if ((PMin1(0)-TOLL <= PMax2(0)) && (PMax1(0)+TOLL >= PMin2(0)))
							{
							block1->AddNearNeigh(0,block2Id);
							//std::cout<<" *"<<block2Id;
							}
						}
					}
				} //fine caso dim ==3
						
				else if (dim == 2)
				{
				if 	((fabs(PMin1(0) - PMax2(0)) <= TOLL) || (fabs(PMax1(0) - PMin2(0)) <= TOLL))
					{
					if ((PMin1(1)-TOLL <= PMax2(1)) && (PMax1(1)+TOLL >= PMin2(1)))
						{
						block1->AddNearNeigh(0,block2Id);
						//std::cout<<block2Id<<" ";
						}
					}

				if ((fabs(PMin1(1) - PMax2(1)) <= TOLL) || (fabs(PMax1(1) - PMin2(1)) <= TOLL))
					{
					if ((PMin1(0)-TOLL <= PMax2(0)) && (PMax1(0)+TOLL >= PMin2(0)))
						{
						block1->AddNearNeigh(0,block2Id);
						//std::cout<<block2Id<<" ";
						}
					}

                                } // fine caso dim == 2 
				
			        }

			for (unsigned int ii = 0; ii < blocks[*pos]->GetBlockChildrenNum();  ii++)
				{
				unsigned int block2Id = blocks[*pos]->GetChildId(ii);
				OctreeBlock<dim> *block2 = blocks[block2Id];
				double delta2 = block2->GetDelta();
		                Point<dim> PMin2 = block2->GetPMin();
		                Point<dim> Center2;
		                for (unsigned int iii = 0; iii < dim; iii++)
		                    Center2(iii) = delta2/2.;
		                Point<dim> PMax2 = 2*Center2;     
		                PMax2 += PMin2;
		                Center2 += PMin2;
				
				if (dim == 3)
				{
				if 	((fabs(PMin1(0) - PMax2(0)) <= TOLL) || (fabs(PMax1(0) - PMin2(0)) <= TOLL))
					{
					if ((PMin1(1)-TOLL <= PMax2(1)) && (PMax1(1)+TOLL >= PMin2(1)))
						{
						if ((PMin1(2)-TOLL <= PMax2(2)) && (PMax1(2)+TOLL >= PMin2(2)))
							{
							block1->AddNearNeigh(0,block2Id);
							//std::cout<<" "<<block2Id;
							}
						}
					}

				if ((fabs(PMin1(1) - PMax2(1)) <= TOLL) || (fabs(PMax1(1) - PMin2(1)) <= TOLL))
					{
					if ((PMin1(0)-TOLL <= PMax2(0)) && (PMax1(0)+TOLL >= PMin2(0)))
						{
						if ((PMin1(2)-TOLL <= PMax2(2)) && (PMax1(2)+TOLL >= PMin2(2)))
							{
							block1->AddNearNeigh(0,block2Id);
							//std::cout<<" "<<block2Id;
							}
						}
					}

				if ((fabs(PMin1(2) - PMax2(2)) <= TOLL) || (fabs(PMax1(2) - PMin2(2)) <= TOLL))
					{
					if ((PMin1(1)-TOLL <= PMax2(1)) && (PMax1(1)+TOLL >= PMin2(1)))
						{
						if ((PMin1(0)-TOLL <= PMax2(0)) && (PMax1(0)+TOLL >= PMin2(0)))
							{
							block1->AddNearNeigh(0,block2Id);
							//std::cout<<" "<<block2Id;
							}
						}
					} 
				} //fine caso dim ==3
						
				else if (dim == 2)
				{
				if 	((fabs(PMin1(0) - PMax2(0)) <= TOLL) || (fabs(PMax1(0) - PMin2(0)) <= TOLL))
					{
					if ((PMin1(1)-TOLL <= PMax2(1)) && (PMax1(1)+TOLL >= PMin2(1)))
						{
						block1->AddNearNeigh(0,block2Id);
						//std::cout<<block2Id<<" ";
						}
					}

				if ((fabs(PMin1(1) - PMax2(1)) <= TOLL) || (fabs(PMax1(1) - PMin2(1)) <= TOLL))
					{
					if ((PMin1(0)-TOLL <= PMax2(0)) && (PMax1(0)+TOLL >= PMin2(0)))
						{
						block1->AddNearNeigh(0,block2Id);
						//std::cout<<block2Id<<" ";
						}
					}

                                } // fine caso dim == 2 
				
				} // fine loop sui figli di un nearest neighbor del padre


			} // fine loop sui nearest neighbors del padre

                        

		if ((block1->GetBlockChildrenNum() == 0))  // if the block is childless we must compute now its nearneigh at all residual levels
			{
			block1->SetNearNeighSize(num_octree_levels-level+1);
			block1->SetIntListSize(num_octree_levels-level+1); // intList is a vector of sets with the same number of members of nearNeigh
			block1->SetNonIntListSize(num_octree_levels-level+1); // nonIntList is a vector of sets with the same number of members of nearNeigh
			
			for (unsigned int subLevel = 1; subLevel < num_octree_levels - level + 1; subLevel++)

				{
				
				std::set <unsigned int> upperLevelNNeighs = block1->GetNearNeighs(subLevel-1);
				for (std::set <unsigned int>::iterator pos = upperLevelNNeighs.begin(); pos != upperLevelNNeighs.end();  pos++)

					{
					if (blocks[*pos]->GetBlockChildrenNum() == 0) // if nearneigh is childless, it will stay a near neigh
						block1->AddNearNeigh(subLevel,*pos);

					for (unsigned int ii = 0; ii < blocks[*pos]->GetBlockChildrenNum();  ii++)
					    {
						unsigned int block2Id = blocks[*pos]->GetChildId(ii);
						OctreeBlock<dim>* block2 = blocks[block2Id];
						double delta2 = block2->GetDelta();
		                		Point<dim> PMin2 = block2->GetPMin();
		                		Point<dim> Center2;
		                		for (unsigned int iii = 0; iii < dim; iii++)
		                    		    Center2(iii) = delta2/2.;
		                		Point<dim> PMax2 = 2*Center2;     
		                		PMax2 += PMin2;
		                		Center2 += PMin2;
						
						if (dim == 3)
						{
						if 	((fabs(PMin1(0) - PMax2(0)) <= TOLL) || (fabs(PMax1(0) - PMin2(0)) <= TOLL))
							{
							if ((PMin1(1)-TOLL <= PMax2(1)) && (PMax1(1)+TOLL >= PMin2(1)))
								{
								if ((PMin1(2)-TOLL <= PMax2(2)) && (PMax1(2)+TOLL >= PMin2(2)))
									{
									block1->AddNearNeigh(subLevel,block2Id);
									//std::cout<<block2Id<<" ";
									}
								}
							}

						if ((fabs(PMin1(1) - PMax2(1)) <= TOLL) || (fabs(PMax1(1) - PMin2(1)) <= TOLL))
							{
							if ((PMin1(0)-TOLL <= PMax2(0)) && (PMax1(0)+TOLL >= PMin2(0)))
								{
								if ((PMin1(2)-TOLL <= PMax2(2)) && (PMax1(2)+TOLL >= PMin2(2)))
									{
									block1->AddNearNeigh(subLevel,block2Id);
									//std::cout<<block2Id<<" ";
									}
								}
							}

						if ((fabs(PMin1(2) - PMax2(2)) <= TOLL) || (fabs(PMax1(2) - PMin2(2)) <= TOLL))
							{
							if ((PMin1(1)-TOLL <= PMax2(1)) && (PMax1(1)+TOLL >= PMin2(1)))
								{
								if ((PMin1(0)-TOLL <= PMax2(0)) && (PMax1(0)+TOLL >= PMin2(0)))
									{
									block1->AddNearNeigh(subLevel,block2Id);
									//std::cout<<block2Id<<" ";
									}
								}
							} 
                                                } //fine caso dim ==3
						
						else if (dim == 2)
						{
						if 	((fabs(PMin1(0) - PMax2(0)) <= TOLL) || (fabs(PMax1(0) - PMin2(0)) <= TOLL))
							{
							if ((PMin1(1)-TOLL <= PMax2(1)) && (PMax1(1)+TOLL >= PMin2(1)))
								{
								block1->AddNearNeigh(subLevel,block2Id);
								//std::cout<<block2Id<<" ";
								}
							}

						if ((fabs(PMin1(1) - PMax2(1)) <= TOLL) || (fabs(PMax1(1) - PMin2(1)) <= TOLL))
							{
							if ((PMin1(0)-TOLL <= PMax2(0)) && (PMax1(0)+TOLL >= PMin2(0)))
								{
								block1->AddNearNeigh(subLevel,block2Id);
								//std::cout<<block2Id<<" ";
								}
							}

                                                } // fine caso dim == 2 

					    } // fine loop sui figli di ciascun nearest neighbor del blocco childless


					} // fine loop sui nearest neighbors del blocco childless


				} // fine loop sui subLevels (da quello del blocco childless all'ultimo)


			} // fine if (il blocco e' childless?)



	} // fine loop sui blocchi di un livello

} // fine loop sui livelli

//for printout                      

/*std::cout<<std::endl;
std::cout<<"-------------------------------- "<<std::endl;
std::cout<<"-------------------------------- "<<std::endl;
std::cout<<std::endl;
		
//for printout
for (cell=dh.begin_active();cell!=endc; cell++)
    {
    std::cout<<std::endl;
    std::cout<<"-------------------------------- "<<std::endl;
    std::cout<<"Cell "<<cell<<"  elementPlot(";
    cell->get_dof_indices(local_dof_indices);
    for (unsigned int j = 0; j<local_dof_indices.size(); j++)
        std::cout<<"["<<support_points[local_dof_indices[j]]<<"],";
        std::cout<<"'r')";
    }
			
std::cout<<std::endl;
std::cout<<"-------------------------------- "<<std::endl;
std::cout<<"-------------------------------- "<<std::endl;
std::cout<<std::endl;
std::cout<<std::endl;//*/


// search for interaction list blocks (NearNeigh + NearNeighOfNearNeigh)
// and for non interaction list blocks (nonIntList for a block B is composed by blocks that are children
// of blocks being in intList of B's parent, but are not in intList of B)

for (unsigned int ii = startLevel[1]; ii < endLevel[1] + 1;  ii++) // at level 1, add all blocks to intList
	{
	for (unsigned int jj = startLevel[1]; jj < endLevel[1] + 1;  jj++)
		{
		blocks[ii]->AddBlockToIntList(0,jj);
		}
	}


for (unsigned int level = 2; level < num_octree_levels + 1;  level++) // loop over levels

	{
	for (unsigned int jj = startLevel[level]; jj < endLevel[level] + 1;  jj++) // loop over blocks of each level
		{
		OctreeBlock<dim>* block1 = blocks[jj];
		//std::cout<<"??Out "<<jj<<" "<<block1.GetNearNeighSize()<<std::endl;
		for (unsigned int subLevel = 0; subLevel < block1->NumNearNeighLevels(); subLevel++)
                    {
		    std::set <unsigned int> NNList = block1->GetNearNeighs(subLevel);
		    
		    for (std::set <unsigned int>::iterator pos1 = NNList.begin(); pos1 != NNList.end();  pos1++) //loop over blocks in NN list and get their NNs
		        {
			block1->AddBlockToIntList(subLevel,*pos1);
			}
		    //std::cout<<std::endl<<"Sublevel "<<subLevel<<" elem("<<block1.GetBlockElementsList()[0]<<") NearNeighs: ";
		    std::vector <unsigned int> nodeIds = block1->GetBlockNodeList();
		    //std::cout<<std::endl<<"Level "<<level<<"  Block1: "<<kk<<"  NumElements: "<<block1.GetBlockElementsNum()<<" "<<elemIds.size()<<std::endl;
		    //std::cout<<"Nearest Neighbors Found:"<<std::endl;
		    for (unsigned int pp = 0; pp < nodeIds.size();  pp++)
		        {
			//std::cout<<"Node "<<nodeIds[pp]<<std::endl;
			std::set <unsigned int> doubleNodes = double_nodes_set[nodeIds[pp]];
			for (std::set<unsigned int>::iterator pos = doubleNodes.begin();
			     pos != doubleNodes.end(); pos++)
			     {
			     //std::cout<<"Node Double"<<*pos<<std::endl;
			     //std::vector<cell_it > surrCellIds = dof_to_elems[*pos];
			     for (unsigned int k=0; k < dof_to_elems[*pos].size(); k++)
			         {
				 cell_it cell = dof_to_elems[*pos][k];
				 //std::cout<<cell<<std::endl;
				 for (unsigned int j=0; j < quadPoints[cell].size(); j++)
				     {
				     block1->AddBlockToIntList(subLevel,quad_point_to_block[cell][j][level+subLevel]);
				     }
				 }  
			     }
		        }
                    }
		for (unsigned int subLevel = 0; subLevel < block1->GetNearNeighSize();  subLevel++) // for each block, loop over all sublevels in his NN list (to account for childless blocks)
			{
			// now use intList to compute nonIntList
			std::set <unsigned int> intList = block1->GetIntList(subLevel);
			std::set <unsigned int> parentIntList; // intList at the  previous level
			if (subLevel == 0) // if a block is childless we get its intList at the previous level, otherwise we get its parent's intList
				parentIntList = blocks[block1->GetParentId()]->GetIntList(0);
			else
				parentIntList = block1->GetIntList(subLevel-1);

			for (std::set <unsigned int>::iterator pos1 = parentIntList.begin(); pos1 != parentIntList.end();  pos1++) // loop over blocks in parentIntList
				{
			 	OctreeBlock<dim>* block2 = blocks[*pos1];
			 	if (block2->GetBlockChildrenNum() == 0) // if blocks in parentIntList are childless, don't look for their children, but see if they are in nonIntList
			 		{
			 		if (intList.count(*pos1) == 0) // if these blocks are not in intList
			 			block1->AddBlockToNonIntList(subLevel,*pos1); // then they go in nonIntList
			 		}
			 	else // if blocks in parentIntList are not childless, do the same test on all their children
			 		{
			 		for (unsigned int kk = 0; kk < block2->GetBlockChildrenNum();  kk++) // loop over children of blocks in parentIntList
			 			{
						//std::cout<<"Sublevel "<<subLevel<<" Block1 "<<jj<<"  Block2 "<<*pos1<<" child(kk) "<<block2.GetChildId(kk)<<std::endl;
			 			if (intList.count(block2->GetChildId(kk)) == 0)   // if these blocks are not in intList
			 				block1->AddBlockToNonIntList(subLevel,block2->GetChildId(kk)); // then they go in nonIntList
			 			} // end loop over children of blocks in parentIntList
			 		}
				}	// loop over blocks in parentIntList

			}	// end loop over subLevels of each block's intList

//			for printout
                        
/*			
			std::cout<<std::endl;
                        std::cout<<"-------------------------------- "<<std::endl;
			std::cout<<"Block "<<jj<<": Parent "<<block1->GetParentId();
			std::cout<<"  cubePlot(["<<block1->GetPMin()<<"],"<<block1->GetDelta()<<",'m')"<<std::endl;
			std::cout<<"Quad points of block "<<jj<<": "<<std::endl;
			std::map <cell_it,std::vector<unsigned int> > quadPointsList = block1->GetBlockQuadPointsList();
			typename std::map <cell_it, std::vector<unsigned int> >::iterator it;    
                        for (it = quadPointsList.begin(); it != quadPointsList.end(); it++)
	                    {
			    std::cout<<(*it).first<<"( ";
			    for (unsigned int zz = 0; zz < (*it).second.size();  zz++)
			        std::cout<<(*it).second[zz]<<" ";
			    std::cout<<")  ";
			    }
			std::cout<<std::endl;

			std::cout<<"Nodes of block "<<jj<<" :"<<std::endl;
			std::vector <unsigned int> nodeList = block1->GetBlockNodeList();
			for (unsigned int zz = 0; zz < block1->GetBlockNodesNum();  zz++)
				{
				std::cout<<nodeList.at(zz)<<" ";
				}
			std::cout<<std::endl;

			for (unsigned int subLevel = 0; subLevel < block1->GetNearNeighSize();  subLevel++)
				{
				std::set <unsigned int> NNList = block1->GetNearNeighs(subLevel);
				std::cout<<"NearNeigh for block "<<jj<<" at level "<<level+subLevel<<":"<<std::endl;
				for (std::set <unsigned int>::iterator pos1 = NNList.begin(); pos1 != NNList.end();  pos1++)
					{
					std::cout<<*pos1<<" ";
					}
				std::cout<<std::endl;

				std::set <unsigned int> intList = block1->GetIntList(subLevel);
				std::cout<<"IntList for block "<<jj<<" at level "<<level+subLevel<<":"<<std::endl;
				for (std::set <unsigned int>::iterator pos1 = intList.begin(); pos1 != intList.end();  pos1++)
					{
					std::cout<<*pos1<<" ";
					}
				std::cout<<std::endl;

				std::set <unsigned int> nonIntList = block1->GetNonIntList(subLevel);
				std::cout<<"NonIntList for block "<<jj<<" at level "<<level+subLevel<<":"<<std::endl;
				for (std::set <unsigned int>::iterator pos1 = nonIntList.begin(); pos1 != nonIntList.end();  pos1++)
					{
					std::cout<<*pos1<<" ";
					}
				std::cout<<std::endl;
				}
//*/



		} // end loop over blocks of a level

	} // end loop over levels

//if (blocks.size() > 0)
//   {
//   for (unsigned int ii = 0; ii < num_blocks;  ii++)
//        delete blocks[ii];
//   }

integralCheck.clear(); 
for (unsigned int i = 0; i < dh.n_dofs(); i++)
    {
    for (cell_it cell = dh.begin_active(); cell != endc; ++cell)
        {
	//std::cout<<i<<" "<<cell<<" "<<integralCheck[i][cell]<<std::endl;
	integralCheck[i][cell] = 0;
	}
    }



std::cout<<"Done computing proximity lists for blocks"<<std::endl;
} //end method for octree blocking generation



template <int dim>
void ComputationalDomain<dim>::inflow_mesh_smoothing()

{
std::cout<<"Smoothing inflow surface"<<std::endl;

std::cout<<"done smoothing inflow surface"<<std::endl;
}


template <int dim>
void ComputationalDomain<dim>::restore_tria(std::string fname)
{
  std::ifstream ifile(fname.c_str());
  tria.clear();

  tria.read_flags(ifile);



  tria.restore();
  ifile.close();
  
  dh.distribute_dofs(fe);
  vector_dh.distribute_dofs(vector_fe);
  //DoFRenumbering::Cuthill_McKee(dh);
  //DoFRenumbering::Cuthill_McKee(vector_dh);

  map_points.reinit(vector_dh.n_dofs());

  if(mapping == NULL)
    mapping = new MappingQEulerian<dim-1, Vector<double>, dim>
	      (mapping_degree, map_points, vector_dh);
  
  generate_double_nodes_set();
  compute_min_diameter();
  std::cout<<"Total number of dofs after restore: "<<dh.n_dofs()<<std::endl;
}

template <int dim>
void ComputationalDomain<dim>::dump_tria(std::string fname) const
{
  std::ofstream ofile(fname.c_str());
  tria.write_flags(ofile);
  ofile.close();
}  



template <int dim>
void ComputationalDomain<dim>::update_support_points()
{
  ref_points.resize(vector_dh.n_dofs());
  DoFTools::map_dofs_to_support_points<dim-1, dim>(StaticMappingQ1<2,3>::mapping,
					    vector_dh, ref_points);

  vector_support_points.resize(vector_dh.n_dofs());
  DoFTools::map_dofs_to_support_points<dim-1, dim>( *mapping, vector_dh, vector_support_points);
  
  support_points.resize(dh.n_dofs());
  DoFTools::map_dofs_to_support_points<dim-1, dim>( *mapping, dh, support_points);
}





  template <int dim>
  void ComputationalDomain<dim>::compute_normals()
  {
   typedef typename DoFHandler<dim-1,dim>::active_cell_iterator cell_it;

   SparsityPattern      normals_sparsity_pattern;
   normals_sparsity_pattern.reinit(vector_dh.n_dofs(),
                                   vector_dh.n_dofs(),
                                   vector_dh.max_couplings_between_dofs());
   ConstraintMatrix  vector_constraints;
   vector_constraints.clear();
   DoFTools::make_hanging_node_constraints (vector_dh,vector_constraints);
   vector_constraints.close();
   DoFTools::make_sparsity_pattern (vector_dh, normals_sparsity_pattern, vector_constraints);
   normals_sparsity_pattern.compress();
   Vector<double> vector_normals_solution(vector_dh.n_dofs());
                                   
   SparseMatrix<double> vector_normals_matrix;
   Vector<double> vector_normals_rhs;

   vector_normals_matrix.reinit (normals_sparsity_pattern);
   vector_normals_rhs.reinit(vector_dh.n_dofs());
   vector_normals_solution.reinit(vector_dh.n_dofs());


   FEValues<dim-1,dim> vector_fe_v(*mapping, vector_fe, *quadrature,
			     update_values | update_cell_normal_vectors |  
			     update_JxW_values);

   const unsigned int vector_n_q_points = vector_fe_v.n_quadrature_points;
   const unsigned int   vector_dofs_per_cell   = vector_fe.dofs_per_cell;
   std::vector<unsigned int> vector_local_dof_indices (vector_dofs_per_cell);

   FullMatrix<double>   local_normals_matrix (vector_dofs_per_cell, vector_dofs_per_cell);
   Vector<double>       local_normals_rhs (vector_dofs_per_cell);

   cell_it
   vector_cell = vector_dh.begin_active(),
   vector_endc = vector_dh.end();
            
   for (; vector_cell!=vector_endc; ++vector_cell)
     {
       vector_fe_v.reinit (vector_cell);
       local_normals_matrix = 0;
       local_normals_rhs = 0;
       const std::vector<Point<dim> > &vector_node_normals = vector_fe_v.get_normal_vectors();
       unsigned int comp_i, comp_j;
       
       for (unsigned int q=0; q<vector_n_q_points; ++q)
	 for (unsigned int i=0; i<vector_dofs_per_cell; ++i)
	   {
	     comp_i = vector_fe.system_to_component_index(i).first;
	     for (unsigned int j=0; j<vector_dofs_per_cell; ++j)
	       {
		 comp_j = vector_fe.system_to_component_index(j).first;
		 if (comp_i == comp_j) 
		   {
		     local_normals_matrix(i,j) += vector_fe_v.shape_value(i,q)*
						  vector_fe_v.shape_value(j,q)*
						  vector_fe_v.JxW(q);
		   }
	       }
	   local_normals_rhs(i) += (vector_fe_v.shape_value(i, q)) *
                                    vector_node_normals[q](comp_i) * vector_fe_v.JxW(q);
	   }
        
       vector_cell->get_dof_indices (vector_local_dof_indices);
       
       vector_constraints.distribute_local_to_global
       (local_normals_matrix,
	local_normals_rhs,
	vector_local_dof_indices,
	vector_normals_matrix,
	vector_normals_rhs);
     }



   SparseDirectUMFPACK normals_inverse;
   normals_inverse.initialize(vector_normals_matrix);
   normals_inverse.vmult(vector_normals_solution, vector_normals_rhs);
   vector_constraints.distribute(vector_normals_solution);

   nodes_normals.resize(dh.n_dofs());
 
   for (unsigned int i=0; i<vector_dh.n_dofs()/dim; ++i)
       { 
       for (unsigned int d=0; d<dim; d++)
           nodes_normals[i](d) = vector_normals_solution(3*i+d);
       nodes_normals[i]/= nodes_normals[i].distance(Point<dim>(0.0,0.0,0.0));
       //cout<<i<<" Gradient: "<<node_normals[i]<<endl;
       for (unsigned int d=0; d<dim; d++)
           vector_normals_solution(3*i+d) = nodes_normals[i](d);
       }


  }



template class ComputationalDomain<3>;
