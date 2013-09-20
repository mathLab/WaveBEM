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
//    Authors: Luca Heltai, Cataldo Manigrasso
//
//----------------------------  step-34.cc  ---------------------------


#include "../include/bem_problem.h"
#include "../include/laplace_kernel.h"

				 // @sect4{BEMProblem::BEMProblem and
				 // BEMProblem::read_parameters}
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
BEMProblem<dim>::BEMProblem(ComputationalDomain<dim> &comp_dom,
                            BEMFMA<dim> &fma)
		:
		comp_dom(comp_dom), fma(fma)
{}

template <int dim>
void BEMProblem<dim>::reinit()
{
  const unsigned int n_dofs =  comp_dom.dh.n_dofs();
  
  neumann_matrix.reinit(n_dofs, n_dofs);
  dirichlet_matrix.reinit(n_dofs, n_dofs);  
  system_rhs.reinit(n_dofs);
  sol.reinit(n_dofs);
  alpha.reinit(n_dofs);
  serv_phi.reinit(n_dofs);
  serv_dphi_dn.reinit(n_dofs);
  serv_tmp_rhs.reinit(n_dofs);
  preconditioner_band = 100;
  preconditioner_sparsity_pattern.reinit (n_dofs,n_dofs,preconditioner_band);
  is_preconditioner_initialized = false;
}

template <int dim> 
void BEMProblem<dim>::declare_parameters (ParameterHandler &prm)
{
		    		    
				   // In the solver section, we set
				   // all SolverControl
				   // parameters. The object will then
				   // be fed to the GMRES solver in
				   // the solve_system() function.

  prm.enter_subsection("Solver");
  SolverControl::declare_parameters(prm);
  prm.leave_subsection();

  prm.declare_entry("Solution method", "Direct", 
                    Patterns::Selection("Direct|FMA")); 
}

template <int dim> 
void BEMProblem<dim>::parse_parameters (ParameterHandler &prm)
{

  prm.enter_subsection("Solver");
  solver_control.parse_parameters(prm);
  prm.leave_subsection();
  
  solution_method = prm.get("Solution method");
  
  
}



template <int dim>
void BEMProblem<dim>::assemble_system()
{    
std::cout<<"(Directly) Assembling system matrices"<<std::endl;

  neumann_matrix = 0;
  dirichlet_matrix = 0;  



  std::vector<QGaussOneOverR<2> > sing_quadratures_3d; 
  for(unsigned int i=0; i<comp_dom.fe.dofs_per_cell; ++i)
    sing_quadratures_3d.push_back
      (QGaussOneOverR<2>(comp_dom.singular_quadrature_order,
			 comp_dom.fe.get_unit_support_points()[i], true));
  
    
				   // Next, we initialize an FEValues
				   // object with the quadrature
				   // formula for the integration of
				   // the kernel in non singular
				   // cells. This quadrature is
				   // selected with the parameter
				   // file, and needs to be quite
				   // precise, since the functions we
				   // are integrating are not
				   // polynomial functions.
  FEValues<dim-1,dim> fe_v(*comp_dom.mapping,comp_dom.fe, *comp_dom.quadrature,
			   update_values |
			   update_cell_normal_vectors |
			   update_quadrature_points |
			   update_JxW_values);
    
  const unsigned int n_q_points = fe_v.n_quadrature_points;
    
  std::vector<unsigned int> local_dof_indices(comp_dom.fe.dofs_per_cell);
    
				   // Unlike in finite element
				   // methods, if we use a collocation
				   // boundary element method, then in
				   // each assembly loop we only
				   // assemble the information that
				   // refers to the coupling between
				   // one degree of freedom (the
				   // degree associated with support
				   // point $i$) and the current
				   // cell. This is done using a
				   // vector of fe.dofs_per_cell
				   // elements, which will then be
				   // distributed to the matrix in the
				   // global row $i$. The following
				   // object will hold this
				   // information:
  Vector<double>      local_neumann_matrix_row_i(comp_dom.fe.dofs_per_cell);
  Vector<double>      local_dirichlet_matrix_row_i(comp_dom.fe.dofs_per_cell);
    
				   // Now that we have checked that
				   // the number of vertices is equal
				   // to the number of degrees of
				   // freedom, we construct a vector
				   // of support points which will be
				   // used in the local integrations:
  std::vector<Point<dim> > support_points(comp_dom.dh.n_dofs());
  DoFTools::map_dofs_to_support_points<dim-1, dim>( *comp_dom.mapping, comp_dom.dh, support_points);

  
				   // After doing so, we can start the
				   // integration loop over all cells,
				   // where we first initialize the
				   // FEValues object and get the
				   // values of $\mathbf{\tilde v}$ at
				   // the quadrature points (this
				   // vector field should be constant,
				   // but it doesn't hurt to be more
				   // general):


  cell_it
    cell = comp_dom.dh.begin_active(),
    endc = comp_dom.dh.end();
    
  Point<dim> D;
  double s;
    
  for(cell = comp_dom.dh.begin_active(); cell != endc; ++cell)
    {
      fe_v.reinit(cell);
      cell->get_dof_indices(local_dof_indices);
        
      const std::vector<Point<dim> > &q_points = fe_v.get_quadrature_points();
      const std::vector<Point<dim> > &normals = fe_v.get_normal_vectors();
        
      				       // We then form the integral over
				       // the current cell for all
				       // degrees of freedom (note that
				       // this includes degrees of
				       // freedom not located on the
				       // current cell, a deviation from
				       // the usual finite element
				       // integrals). The integral that
				       // we need to perform is singular
				       // if one of the local degrees of
				       // freedom is the same as the
				       // support point $i$. A the
				       // beginning of the loop we
				       // therefore check wether this is
				       // the case, and we store which
				       // one is the singular index:
      for(unsigned int i=0; i<comp_dom.dh.n_dofs() ; ++i)
	{
          
	  local_neumann_matrix_row_i = 0;
	  local_dirichlet_matrix_row_i = 0;
            
	  bool is_singular = false; 
	  unsigned int singular_index = numbers::invalid_unsigned_int;
            
	  for(unsigned int j=0; j<comp_dom.fe.dofs_per_cell; ++j) 
	    //if(local_dof_indices[j] == i)
	    if(comp_dom.double_nodes_set[i].count(local_dof_indices[j]) > 0)
	      {
		singular_index = j;
		is_singular = true;
		break;
	      }

					   // We then perform the
					   // integral. If the index $i$
					   // is not one of the local
					   // degrees of freedom, we
					   // simply have to add the
					   // single layer terms to the
					   // right hand side, and the
					   // double layer terms to the
					   // matrix:
	  if(is_singular == false)
	    {
	      for(unsigned int q=0; q<n_q_points; ++q)
		{
		  const Point<dim> R = q_points[q] - support_points[i];
                  LaplaceKernel::kernels(R, D, s);
    
		  for (unsigned int j=0; j<comp_dom.fe.dofs_per_cell; ++j)
		      {
		    local_neumann_matrix_row_i(j) += ( ( D * 
							 normals[q] ) *
						       fe_v.shape_value(j,q) *
						       fe_v.JxW(q)       );
		    local_dirichlet_matrix_row_i(j) += ( s * 
							 fe_v.shape_value(j,q) *
							 fe_v.JxW(q) );
						 
		      }
		}
	    } else {
					     // Now we treat the more
					     // delicate case. If we
					     // are here, this means
					     // that the cell that
					     // runs on the $j$ index
					     // contains
					     // support_point[i]. In
					     // this case both the
					     // single and the double
					     // layer potential are
					     // singular, and they
					     // require special
					     // treatment.
					     //
					     // Whenever the
					     // integration is
					     // performed with the
					     // singularity inside the
					     // given cell, then a
					     // special quadrature
					     // formula is used that
					     // allows one to
					     // integrate arbitrary
					     // functions against a
					     // singular weight on the
					     // reference cell.
				     	     // Notice that singular
				     	     // integration requires a
				     	     // careful selection of
				     	     // the quadrature
				     	     // rules. In particular
				     	     // the deal.II library
				     	     // provides quadrature
				     	     // rules which are
				     	     // taylored for
				     	     // logarithmic
				     	     // singularities
				     	     // (QGaussLog,
				     	     // QGaussLogR), as well
				     	     // as for 1/R
				     	     // singularities
				     	     // (QGaussOneOverR).
				     	     //
				     	     // Singular integration
				     	     // is typically obtained
				     	     // by constructing
				     	     // weighted quadrature
				     	     // formulas with singular
				     	     // weights, so that it is
				     	     // possible to write
				     	     //
				     	     // \f[
				     	     //   \int_K f(x) s(x) dx = \sum_{i=1}^N w_i f(q_i)
				     	     // \f]
				     	     //
				     	     // where $s(x)$ is a given
				     	     // singularity, and the weights
				     	     // and quadrature points
				     	     // $w_i,q_i$ are carefully
				     	     // selected to make the formula
				     	     // above an equality for a
				     	     // certain class of functions
				     	     // $f(x)$.
				     	     //
				     	     // In all the finite
				     	     // element examples we
				     	     // have seen so far, the
				     	     // weight of the
				     	     // quadrature itself
				     	     // (namely, the function
				     	     // $s(x)$), was always
				     	     // constantly equal to 1.
				     	     // For singular
				     	     // integration, we have
				     	     // two choices: we can
				     	     // use the definition
				     	     // above, factoring out
				     	     // the singularity from
				     	     // the integrand (i.e.,
				     	     // integrating $f(x)$
				     	     // with the special
				     	     // quadrature rule), or
				     	     // we can ask the
				     	     // quadrature rule to
				     	     // "normalize" the
				     	     // weights $w_i$ with
				     	     // $s(q_i)$:
				     	     //
				     	     // \f[
				     	     //   \int_K f(x) s(x) dx =
				     	     //   \int_K g(x) dx = \sum_{i=1}^N \frac{w_i}{s(q_i)} g(q_i)
				     	     // \f]
				     	     //
				     	     // We use this second
				     	     // option, through the @p
				     	     // factor_out_singularity
				     	     // parameter of both
				     	     // QGaussLogR and
				     	     // QGaussOneOverR.
				     	     //
				     	     // These integrals are
				     	     // somewhat delicate,
				     	     // especially in two
				     	     // dimensions, due to the
				     	     // transformation from
				     	     // the real to the
				     	     // reference cell, where
				     	     // the variable of
				     	     // integration is scaled
				     	     // with the determinant
				     	     // of the transformation.
				     	     //
				     	     // In two dimensions this
				     	     // process does not
				     	     // result only in a
				     	     // factor appearing as a
				     	     // constant factor on the
				     	     // entire integral, but
				     	     // also on an additional
				     	     // integral alltogether
				     	     // that needs to be
				     	     // evaluated:
				     	     // 
				     	     // \f[
				     	     //  \int_0^1 f(x)\ln(x/\alpha) dx =
				     	     //  \int_0^1 f(x)\ln(x) dx - \int_0^1 f(x) \ln(\alpha) dx.
				     	     // \f]
				     	     //
				     	     // This process is taken care of by
				     	     // the constructor of the QGaussLogR
				     	     // class, which adds additional
				     	     // quadrature points and weights to
				     	     // take into consideration also the
				     	     // second part of the integral.
				     	     // 
				     	     // A similar reasoning
				     	     // should be done in the
				     	     // three dimensional
				     	     // case, since the
				     	     // singular quadrature is
				     	     // taylored on the
				     	     // inverse of the radius
				     	     // $r$ in the reference
				     	     // cell, while our
				     	     // singular function
				     	     // lives in real space,
				     	     // however in the three
				     	     // dimensional case
				     	     // everything is simpler
				     	     // because the
				     	     // singularity scales
				     	     // linearly with the
				     	     // determinant of the
				     	     // transformation. This
				     	     // allows us to build the
				     	     // singular two
				     	     // dimensional quadrature
				     	     // rules once and for all
				     	     // outside the loop over
				     	     // all cells, using only
				     	     // a pointer where needed.
					     //
					     // Notice that in one
					     // dimensional
					     // integration this is
					     // not possible, since we
					     // need to know the
					     // scaling parameter for
					     // the quadrature, which
					     // is not known a
					     // priori. Here, the
					     // quadrature rule itself
					     // depends also on the
					     // size of the current
					     // cell. For this reason,
					     // it is necessary to
					     // create a new
					     // quadrature for each
					     // singular
					     // integration. Since we
					     // create it using the
					     // new operator of C++,
					     // we also need to
					     // destroy it using the
					     // dual of new:
					     // delete. This is done
					     // at the end, and only
					     // if dim == 2.
					     //
					     // Putting all this into a
					     // dimension independent
					     // framework requires a little
					     // trick. The problem is that,
					     // depending on dimension, we'd
					     // like to either assign a
					     // QGaussLogR<1> or a
					     // QGaussOneOverR<2> to a
					     // Quadrature<dim-1>. C++
					     // doesn't allow this right
					     // away, and neither is a
					     // static_cast
					     // possible. However, we can
					     // attempt a dynamic_cast: the
					     // implementation will then
					     // look up at run time whether
					     // the conversion is possible
					     // (which we <em>know</em> it
					     // is) and if that isn't the
					     // case simply return a null
					     // pointer. To be sure we can
					     // then add a safety check at
					     // the end:
	    Assert(singular_index != numbers::invalid_unsigned_int,
		   ExcInternalError());

	    const Quadrature<dim-1> *
	      singular_quadrature
	      = (dim == 2
		 ?
		 dynamic_cast<Quadrature<dim-1>*>(
		   new QGaussLogR<1>(comp_dom.singular_quadrature_order,
				     Point<1>((double)singular_index),
				     1./cell->measure(), true))
		 :
		 (dim == 3
		  ?
		  dynamic_cast<Quadrature<dim-1>*>(
		    &sing_quadratures_3d[singular_index])
		  :
		  0));
	    Assert(singular_quadrature, ExcInternalError());
                        
	    FEValues<dim-1,dim> fe_v_singular (*comp_dom.mapping, comp_dom.fe, *singular_quadrature, 
					       update_jacobians |
					       update_values |
					       update_cell_normal_vectors |
					       update_quadrature_points );

	    fe_v_singular.reinit(cell);
                    
	    const std::vector<Point<dim> > &singular_normals = fe_v_singular.get_normal_vectors();
	    const std::vector<Point<dim> > &singular_q_points = fe_v_singular.get_quadrature_points();
                           
	    for(unsigned int q=0; q<singular_quadrature->size(); ++q)
	      { 
		const Point<dim> R = singular_q_points[q] - support_points[i];
                LaplaceKernel::kernels(R, D, s);   
                		       
		for(unsigned int j=0; j<comp_dom.fe.dofs_per_cell; ++j) {
		    local_neumann_matrix_row_i(j) += (( D *
							singular_normals[q])                *
						      fe_v_singular.shape_value(j,q)        *
						      fe_v_singular.JxW(q)       );

		    local_dirichlet_matrix_row_i(j) += ( s   * 
							 fe_v_singular.shape_value(j,q)		  *
							 fe_v_singular.JxW(q) );	
		}
	      }
	    if(dim==2) 
	      delete singular_quadrature;
	  }

					   // Finally, we need to add
					   // the contributions of the
					   // current cell to the
					   // global matrix.
	  for(unsigned int j=0; j<comp_dom.fe.dofs_per_cell; ++j) {
	      neumann_matrix(i,local_dof_indices[j]) 
		+= local_neumann_matrix_row_i(j);
	      dirichlet_matrix(i,local_dof_indices[j]) 
		+= local_dirichlet_matrix_row_i(j);
	  }
	}
    }

				   // The second part of the integral
				   // operator is the term
				   // $\alpha(\mathbf{x}_i)
				   // \phi_j(\mathbf{x}_i)$. Since we
				   // use a collocation scheme,
				   // $\phi_j(\mathbf{x}_i)=\delta_{ij}$
				   // and the corresponding matrix is
				   // a diagonal one with entries
				   // equal to $\alpha(\mathbf{x}_i)$.
    
				   // One quick way to compute this
				   // diagonal matrix of the solid
				   // angles, is to use the Neumann
				   // matrix itself. It is enough to
				   // multiply the matrix with a
				   // vector of elements all equal to
				   // -1, to get the diagonal matrix
				   // of the alpha angles, or solid
				   // angles (see the formula in the
				   // introduction for this). The
				   // result is then added back onto
				   // the system matrix object to
				   // yield the final form of the
				   // matrix:

/*
  std::cout<<"Neumann"<<std::endl;
  for (unsigned int i = 0; i < comp_dom.dh.n_dofs(); i++)
      {
      for (unsigned int j = 0; j < comp_dom.dh.n_dofs(); j++)
          {
          std::cout<<neumann_matrix(i,j)<<" ";
          }
      std::cout<<std::endl;
      }
        


  std::cout<<"Dirichlet"<<std::endl;
  for (unsigned int i = 0; i < comp_dom.dh.n_dofs(); i++)
      {
      for (unsigned int j = 0; j < comp_dom.dh.n_dofs(); j++)
          {
          std::cout<<dirichlet_matrix(i,j)<<" ";
          }
      std::cout<<std::endl;
      }
      //*/
      std::cout<<"done assembling system matrices"<<std::endl;
}



template <int dim>
void BEMProblem<dim>::compute_alpha()
{
  static Vector<double> ones, zeros, dummy;
  if(ones.size() != comp_dom.dh.n_dofs()) 
    {
      ones.reinit(comp_dom.dh.n_dofs());
      ones.add(-1.);
      zeros.reinit(comp_dom.dh.n_dofs());
      dummy.reinit(comp_dom.dh.n_dofs());
    }


  if (solution_method == "Direct")
     {
     neumann_matrix.vmult(alpha, ones);
     }
  else
     {
     fma.generate_multipole_expansions(ones,zeros);
     fma.multipole_matr_vect_products(ones,zeros,alpha,dummy);
     }   

  //alpha.print(std::cout);
}

template <int dim>
void BEMProblem<dim>::vmult(Vector<double> &dst, const Vector<double> &src) const {
  static Vector<double> phi(src.size());
  static Vector<double> dphi_dn(src.size());

  Vector <double> matrVectProdN;
  Vector <double> matrVectProdD;
  
  matrVectProdN.reinit(src.size());
  matrVectProdD.reinit(src.size());
  
  dst = 0;
  
  if(phi.size() != src.size()) {
    phi.reinit(src.size());
    dphi_dn.reinit(src.size());
  }
  phi = src;
  dphi_dn = src;
  //phi.scale(comp_dom.surface_nodes);
  //dphi_dn.scale(comp_dom.other_nodes);
  phi.scale(comp_dom.other_nodes);
  dphi_dn.scale(comp_dom.surface_nodes);
         
  if (solution_method == "Direct")
     {
     dirichlet_matrix.vmult(dst, dphi_dn);
     dst *= -1;
     neumann_matrix.vmult_add(dst, phi);
     phi.scale(alpha);
     dst += phi;
  
     }
  else
     {
     fma.generate_multipole_expansions(phi,dphi_dn);
     fma.multipole_matr_vect_products(phi,dphi_dn,matrVectProdN,matrVectProdD);
     phi.scale(alpha);
     dst += matrVectProdD;
     dst *= -1;
     dst += matrVectProdN;
     dst += phi;
     }
  // in fully neumann bc case, we have to rescale the vector to have a zero mean
  // one
  if (comp_dom.surface_nodes.linfty_norm() < 1e-10)
     dst.add(-dst.l2_norm());

}


template <int dim>
void BEMProblem<dim>::compute_rhs(Vector<double> &dst, const Vector<double> &src) const {
  static Vector<double> phi(src.size());
  static Vector<double> dphi_dn(src.size());
  static Vector<double> rhs(src.size());
  
  static Vector <double> matrVectProdN;
  static Vector <double> matrVectProdD;
  
  if(phi.size() != src.size()) {
    phi.reinit(src.size());
    dphi_dn.reinit(src.size());
    matrVectProdN.reinit(src.size());
    matrVectProdD.reinit(src.size());  
  }  
  
  phi = src;
  dphi_dn = src;

  phi.scale(comp_dom.surface_nodes);
  dphi_dn.scale(comp_dom.other_nodes);
  
  //for (unsigned int ppp = 0; ppp < src.size(); ++ppp)
      //std::cout<<ppp<<"  phi(ppp)="<<phi(ppp)<<"  |||   dphi_dn(ppp)="<<dphi_dn(ppp)<<std::endl;
    
  if (solution_method == "Direct")
     {
     neumann_matrix.vmult(dst, phi);
     phi.scale(alpha);
     dst += phi;
     dst *= -1;
     dirichlet_matrix.vmult_add(dst, dphi_dn);
     }
  else
     {
     fma.generate_multipole_expansions(phi,dphi_dn);
     fma.multipole_matr_vect_products(phi,dphi_dn,matrVectProdN,matrVectProdD);
     phi.scale(alpha);
     dst += matrVectProdN;
     dst += phi;
     dst *= -1;
     dst += matrVectProdD;
     }   
     
/*  
  // here we correct the right hand side so as to account for the presence of double and
  // triple dofs
  
  // we start looping on the dofs   
  for (unsigned int i=0; i <src.size(); i++)
      {
      // in the next line we compute the "first" among the set of double nodes: this node
      // is the first dirichlet node in the set, and if no dirichlet node is there, we get the
      // first neumann node
       
      std::set<unsigned int> doubles = comp_dom.double_nodes_set[i];
      unsigned int firstOfDoubles = *doubles.begin();
      for (std::set<unsigned int>::iterator it = doubles.begin() ; it != doubles.end(); it++ )
          {
          if (comp_dom.surface_nodes(*it) == 1)
             {
             firstOfDoubles = *it;
             break;
             }
          }
      
      // for each set of double nodes, we will perform the correction only once, and precisely
      // when the current node is the first of the set
      if (i == firstOfDoubles)
	 {
	 // the vector entry corresponding to the first node of the set does not need modification,
	 // thus we erase ti form the set
	 doubles.erase(i);
	 
	 // if the current (first) node is a dirichlet node, for all its neumann doubles we will
	 // impose that the potential is equal to that of the first node: this means that in the
	 // rhs we will put the potential value of the first node (which is prescribed by bc) 
	 if (comp_dom.surface_nodes(i) == 1)
	    {
	    for (std::set<unsigned int>::iterator it = doubles.begin() ; it != doubles.end(); it++ )
	        {
		if (comp_dom.surface_nodes(*it) == 1)
		   {
		   
		   }
		else
		   {
	           dst(*it) = phi(i)/alpha(i);
		   }
	        }
	    }
	    
	 // if the current (first) node is a neumann node, for all its doubles we will impose that
	 // the potential is equal to that of the first node: this means that in the rhs we will
	 // put 0  	    
	 if (comp_dom.surface_nodes(i) == 0)
	    {
	    for (std::set<unsigned int>::iterator it = doubles.begin() ; it != doubles.end(); it++ )
	        {
	        dst(*it) = 0;
	        }
	    }
	   
	   
	 }    
      }
    
    
  //*/  
/*    for (unsigned int i=0; i <src.size(); i++)
     {
     if (comp_dom.double_nodes_set[i].size() > 1 && comp_dom.surface_nodes(i) == 0)
        {
	std::set<unsigned int> doubles = comp_dom.double_nodes_set[i];
	doubles.erase(i);
	unsigned int doubleSurIndex = 0;
	double a = 0;
        for (std::set<unsigned int>::iterator it = doubles.begin() ; it != doubles.end(); it++ )
	    {
	    if (comp_dom.surface_nodes(*it) == 1)
	        {
		doubleSurIndex = *it;
		a += comp_dom.surface_nodes(*it);
		}
	    }
	if (a > 0)
	    dst(i) = phi(doubleSurIndex)/alpha(doubleSurIndex);
	else
	    {
	    std::set<unsigned int>::iterator it = doubles.begin();
	    if ( i > *it)
	       dst(i) = 0;
	    }
        }
     }//*/

}





				 // @sect4{BEMProblem::solve_system}

				 // The next function simply solves
				 // the linear system.
template <int dim>
void BEMProblem<dim>::solve_system(Vector<double> &phi, Vector<double> &dphi_dn,
                                   const Vector<double> &tmp_rhs)
{	

   SolverGMRES<Vector<double> > solver (solver_control,
        SolverGMRES<Vector<double> >::AdditionalData(100));
   
   system_rhs = 0;
   sol = 0; 
   alpha = 0;

   compute_alpha();
   
   //std::cout<<"alpha "<<std::endl;
   //for (unsigned int i = 0; i < alpha.size(); i++)
   //std::cout<<i<<" "<<alpha(i)<<std::endl;
            
   compute_rhs(system_rhs, tmp_rhs);

   compute_constraints(constraints, tmp_rhs);
   ConstrainedOperator<Vector<double>, BEMProblem<dim> >
       cc(*this, constraints);

   cc.distribute_rhs(system_rhs);   

   if (solution_method == "Direct")
      {
      //SparseDirectUMFPACK &inv = fma.FMA_preconditioner(alpha);
      //solver.solve (*this, sol, system_rhs, inv);
      assemble_preconditioner();
					   //solver.solve (cc, sol, system_rhs, PreconditionIdentity());
      solver.solve (cc, sol, system_rhs, preconditioner);
      }
    else
      {
      SparseDirectUMFPACK &inv = fma.FMA_preconditioner(alpha);
      solver.solve (*this, sol, system_rhs, inv);
				// solver.solve (*this, sol, system_rhs, PreconditionIdentity());
      }   
    
    //std::cout<<"sol = ["; 
    //for (unsigned int i = 0; i < comp_dom.dh.n_dofs(); i++)
    //    std::cout<<sol(i)<<"; ";
    //std::cout<<"];"<<std::endl; 

   
   
  for (unsigned int i=0; i <comp_dom.surface_nodes.size(); i++)
      {
      if (comp_dom.surface_nodes(i) == 0)
         {
         phi(i) = sol(i);
	 }
	 else
	 {
         dphi_dn(i) = sol(i); 
         }
      }

      //for (unsigned int i=0;i<comp_dom.dh.n_dofs();++i)
      // cout<<i<<" "<<tmp_rhs(i)<<" "<<dphi_dn(i)<<" "<<phi(i)<<" "<<comp_dom.surface_nodes(i)<<endl;

   //std::cout<<"sol "<<std::endl;
   //for (unsigned int i = 0; i < sol.size(); i++)
   //    {
       //std::cout<<i<<" "<<sol(i)<<" ";
       //std::set<unsigned int> doubles = comp_dom.double_nodes_set[i];
       // for (std::set<unsigned int>::iterator it = doubles.begin() ; it != doubles.end(); it++ )
       //    std::cout<<*it<<"("<<comp_dom.surface_nodes(*it)<<") ";
       //std::cout<<"phi "<<phi(i)<<"  dphi_dn "<<dphi_dn(i);    
       //std::cout<<std::endl;
   //    }
       
}



				 // @sect4{BEMProblem::solve_system}

				 // The next function simply solves
				 // the linear system.
template <int dim>
void BEMProblem<dim>::residual(Vector<double> &res, const Vector<double> &phi,
                                                    const Vector<double> &dphi_dn)
{	

   alpha = 0;
   res = 0;
   sol = 0;
   serv_phi = phi;
   serv_dphi_dn = dphi_dn;
   serv_tmp_rhs = 0;

   compute_alpha();

   serv_dphi_dn.scale(comp_dom.other_nodes);
   serv_phi.scale(comp_dom.surface_nodes);
   serv_tmp_rhs.add(serv_dphi_dn);
   serv_tmp_rhs.add(serv_phi);

   //for (unsigned int i=0;i<comp_dom.dh.n_dofs();++i)
   //    cout<<i<<" "<<serv_tmp_rhs(i)<<" "<<serv_dphi_dn(i)<<" "<<serv_phi(i)<<" "<<comp_dom.surface_nodes(i)<<endl;
            
   Vector<double> rrhs(comp_dom.dh.n_dofs());
   compute_rhs(rrhs, serv_tmp_rhs);
   //for (unsigned int i=0;i<comp_dom.dh.n_dofs();++i)
   //    cout<<i<<" "<<rrhs(i)<<endl;
   compute_constraints(constraints, serv_tmp_rhs);
   ConstrainedOperator<Vector<double>, BEMProblem<dim> >
       cc(*this, constraints);
   cc.distribute_rhs(rrhs);   


//for (unsigned int i=0;i<comp_dom.dh.n_dofs();++i)
//       cout<<i<<"--->  "<<rrhs(i)<<" vs "<<system_rhs(i)<<"  diff: "<<rrhs(i)-system_rhs(i)<<endl;

   rrhs *= -1;
   
   sol = 0;

   serv_phi = phi;
   serv_dphi_dn = dphi_dn;
   serv_dphi_dn.scale(comp_dom.surface_nodes);  
   serv_phi.scale(comp_dom.other_nodes);  
   sol.add(serv_dphi_dn);
   sol.add(serv_phi);  

   if (solution_method == "Direct")
      {
      cc.vmult(res,sol);
      res += rrhs;
      }
    else
      {
      //SparseDirectUMFPACK &inv = fma.FMA_preconditioner(alpha);
      //solver.solve (*this, sol, system_rhs, inv);
				// solver.solve (*this, sol, system_rhs, PreconditionIdentity());
      }   
   //cout<<"Residual test: "<<res.l2_norm()<<" £££££ "<<rrhs.l2_norm()<<endl;          
}



				 // @sect4{BEMProblem::compute_errors}

				 // This method performs a Bem resolution,
				 // either in a direct or multipole method
template <int dim>
void BEMProblem<dim>::solve(Vector<double> &phi, Vector<double> &dphi_dn,
                            const Vector<double> &tmp_rhs)
{
  
  if (solution_method == "Direct")
     {
     assemble_system();
     }
  else
     {
     comp_dom.generate_octree_blocking();
     fma.direct_integrals();
     fma.multipole_integrals();
     }
  
  solve_system(phi,dphi_dn,tmp_rhs);
  
}


template <int dim>
void BEMProblem<dim>::compute_constraints(ConstraintMatrix &c, const Vector<double> &tmp_rhs)

{
  // we start clearing the constraint matrix
  c.clear();

  // here we prepare the constraint matrix so as to account for the presence hanging
  // nodes

  DoFTools::make_hanging_node_constraints (comp_dom.dh,c);


  // here we prepare the constraint matrix so as to account for the presence of double and
  // triple dofs
  
  // we start looping on the dofs   
  for (unsigned int i=0; i <tmp_rhs.size(); i++)
      {
      // in the next line we compute the "first" among the set of double nodes: this node
      // is the first dirichlet node in the set, and if no dirichlet node is there, we get the
      // first neumann node
       
      std::set<unsigned int> doubles = comp_dom.double_nodes_set[i];
      unsigned int firstOfDoubles = *doubles.begin();
      for (std::set<unsigned int>::iterator it = doubles.begin() ; it != doubles.end(); it++ )
          {
          if (comp_dom.surface_nodes(*it) == 1)
             {
             firstOfDoubles = *it;
             break;
             }
          }
      
      // for each set of double nodes, we will perform the correction only once, and precisely
      // when the current node is the first of the set
      if (i == firstOfDoubles)
	 {
	 // the vector entry corresponding to the first node of the set does not need modification,
	 // thus we erase ti form the set
	 doubles.erase(i);
	 
	 // if the current (first) node is a dirichlet node, for all its neumann doubles we will
	 // impose that the potential is equal to that of the first node: this means that in the
	 // matrix vector product we will put the potential value of the double node 
	 if (comp_dom.surface_nodes(i) == 1)
	    {
	    for (std::set<unsigned int>::iterator it = doubles.begin() ; it != doubles.end(); it++ )
	        {
		if (comp_dom.surface_nodes(*it) == 1)
		   {
                   // to be checked seven(hundred) more times
                   // this is the dirichlet-dirichlet case: here we impose that
                   // dphi_dn on the two (or more) sides is equal. thus we can only
                   // deal with flat dirichlet-dirichlet boundaries (good enough
                   // for simmetry/wake line on free surface)  
                   c.add_line(*it);
                   c.add_entry(*it,i,1);		   
		   }
		else
		   {
                   c.add_line(*it);
                   c.set_inhomogeneity(*it,tmp_rhs(i));
	           //dst(*it) = phi(*it)/alpha(*it);
		   }
	        }
	    }
	    
	 // if the current (first) node is a neumann node, for all its doubles we will impose that
	 // the potential is equal to that of the first node: this means that in the matrix vector
	 // product we will put the difference between the potential at the fist node in the doubles
	 // set, and the current double node  	    
	 if (comp_dom.surface_nodes(i) == 0)
	    {
	    for (std::set<unsigned int>::iterator it = doubles.begin() ; it != doubles.end(); it++ )
	        {
                c.add_line(*it);
                c.add_entry(*it,i,1);
	        //dst(*it) = phi(*it)/alpha(*it)-phi(i)/alpha(i);
	        }
	    }
	   
	   
	 }    
      }

  c.close();
}

template<int dim>
void BEMProblem<dim>::assemble_preconditioner()
{

  
  static SparseMatrix<double> band_system;

  if(is_preconditioner_initialized == false)
    {
    for (int i=0; (unsigned int) i<comp_dom.dh.n_dofs(); ++i)
        for (int j=std::max(i-preconditioner_band/2,0); j<std::min(i+preconditioner_band/2,(int)comp_dom.dh.n_dofs()); ++j)
	    preconditioner_sparsity_pattern.add((unsigned int) j,(unsigned int) i);
    preconditioner_sparsity_pattern.compress();
    band_system.reinit(preconditioner_sparsity_pattern);
    is_preconditioner_initialized = true;
    }
  else
    band_system = 0;
  
  for(int i=0; (unsigned int) i<comp_dom.dh.n_dofs(); ++i) 
    {
      if(constraints.is_constrained(i))
	band_system.set((unsigned int)i,(unsigned int) i, 1);
      for(int j=std::max(i-preconditioner_band/2,0); j<std::min(i+preconditioner_band/2,(int)comp_dom.dh.n_dofs()); ++j)
	{
	  if(constraints.is_constrained(j) == false) 
	    {
	      if(comp_dom.surface_nodes(i) == 0) 
		{
						   // Nodo di Dirichlet
		  band_system.set(j,i,neumann_matrix(j,i));
		  if(i == j)
		    band_system.add((unsigned int) j,(unsigned int) i, alpha(i));
		}
	      else
		band_system.set(j,i,-dirichlet_matrix(j,i));
	    }
	}
    }
  
  preconditioner.initialize(band_system);
      
}

template class BEMProblem<3>;
