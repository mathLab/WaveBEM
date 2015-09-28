				 // The program starts with including a bunch
				 // of include files that we will use in the
				 // various parts of the program. Most of them
				 // have been discussed in previous tutorials
				 // already:
#ifndef restart_nonlinear_problem_diff_h
#define restart_nonlinear_problem_diff_h


#include <free_surface.h>
#include <numerical_towing_tank.h>
#include <newton_argument.h>


class RestartNonlinearProblemDiff:
public NewtonArgument 
{
public:
  RestartNonlinearProblemDiff(FreeSurface<3> &free_surface,
                              const NumericalTowingTank &comp_dom,
                              const double &time, 
                              Vector<double> &free_surf_y,
                              Vector<double> &free_surf_y_dot,
                              const SparseMatrix<double> &free_surf_jacobian_dot) :
		              free_surface(free_surface),
                              comp_dom(comp_dom),
                              time(time),
                              free_surf_y(free_surf_y),
                              free_surf_y_dot(free_surf_y_dot),
                              free_surf_jacobian_dot(free_surf_jacobian_dot)
  {
  // we first need to count the nodes of the water line. for each of them, we will have a vertical (z)
  // component, and a horizontal (y) component

  water_line_indices.clear();
  bow_stern_indices.clear();
  water_indices.clear();
  phi_indices.clear();
  rigid_modes_indices.clear();

  if (!comp_dom.no_boat)
     for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
         {
         if ( (comp_dom.flags[i] & water) &&
              (comp_dom.flags[i] & near_boat) &&
              !(comp_dom.flags[i] & transom_on_water) &&
              (comp_dom.moving_point_ids[3] != i) &&
              (comp_dom.moving_point_ids[4] != i) &&
              (comp_dom.moving_point_ids[5] != i) &&
              (comp_dom.moving_point_ids[6] != i) )
            {
            //cout<<"WL "<<i<<" ("<<3*i+1<<" "<<3*i+2<<")"<<endl;
            water_line_indices.insert(i);
            }
         }

     //this takes care of the bow and stern nodes 
  if (!comp_dom.no_boat)
     for (unsigned int k=3; k<7; ++k)
         {
         unsigned int i = comp_dom.moving_point_ids[k];
         bow_stern_indices.insert(i);
         //cout<<"BS "<<i<<" ("<<3*i<<" "<<3*i+1<<" "<<3*i+2<<")"<<endl;
         }              

  for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
      {
      if ( (comp_dom.flags[i] & water) &&
           (comp_dom.flags[i] & near_boat)==0 &&
           (comp_dom.flags[i] & transom_on_water)==0 &&
           (comp_dom.vector_constraints.is_constrained(3*i+2))==0 )
         {
         //cout<<"W "<<i<<" ("<<3*i+2<<")"<<endl;
         water_indices.insert(i);
         }
      }

        // let's take care of the potential dofs
      for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
	  {
          if (comp_dom.flags[i] & water)
             // if node is on free surface, component is differential
             {// unless node is a transom_on_water node . then it's algebraic
             if ( (comp_dom.vector_flags[3*i] & transom_on_water) ||
                  (comp_dom.vector_constraints.is_constrained(3*i)) )
                {
                }
             else
                {
                phi_indices.insert(i);
                }
             }
          }
  
        // we now add the rigid modes dofs
  for (unsigned int d=0; d<6; ++d)
      {
      rigid_modes_indices.insert(d);
      }

  free_surf_jac_x_delta.reinit(free_surf_y.size());
  free_surf_res.reinit(free_surf_y.size());



  unsigned int count = 0;
  for (std::set <unsigned int>::iterator pos = water_line_indices.begin(); pos != water_line_indices.end(); ++pos)
      {
      unsigned int i=*pos;
      indices_map[3*i+1] = count;
      count++;
      indices_map[3*i+2] = count;
      count++;
      }
  for (std::set <unsigned int>::iterator pos = bow_stern_indices.begin(); pos != bow_stern_indices.end(); ++pos)
      {
      unsigned int i=*pos;
      Point<3> t(comp_dom.edges_tangents[3*i],comp_dom.edges_tangents[3*i+1],comp_dom.edges_tangents[3*i+2]);
      indices_map[3*i] = count;
      count++;
      indices_map[3*i+1] = count;
      count++;
      indices_map[3*i+2] = count;
      count++;
      }  
  for (std::set <unsigned int>::iterator pos = water_indices.begin(); pos != water_indices.end(); ++pos)
      {
      unsigned int i=*pos;
      indices_map[3*i+2] = count;
      count++;
      }

  for (std::set <unsigned int>::iterator pos = phi_indices.begin(); pos != phi_indices.end(); ++pos)
      {
      unsigned int i=*pos;
      indices_map[i+comp_dom.vector_dh.n_dofs()] = count;
      count++;
      }  

  for (std::set <unsigned int>::iterator pos = rigid_modes_indices.begin(); pos != rigid_modes_indices.end(); ++pos)
      {
      unsigned int i=*pos;
      indices_map[i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()] = count;
      count++;
      }


  //cout<<" MAP "<<endl;
  //for (std::map <unsigned int,unsigned int>::iterator pos = indices_map.begin(); pos != indices_map.end(); ++pos)
  //    {
  //    cout<<"Restart: "<<pos->second<<" --->  DAE: "<<pos->first<<endl;
  //    }
  //cout<<" REPORT "<<endl;
  //cout<<count<<endl;
  //cout<<2*water_line_indices.size()<<endl;
  //cout<<3*bow_stern_indices.size()<<endl;
  //cout<<water_indices.size()<<endl;
  //cout<<phi_indices.size()<<endl;
  //cout<<2*water_line_indices.size()+3*bow_stern_indices.size()+water_indices.size()+phi_indices.size()+dphi_dn_indices.size()<<endl;
  std::vector<unsigned int> line_lengths(2*water_line_indices.size()+
                                         3*bow_stern_indices.size()+
                                         water_indices.size()+
                                         phi_indices.size()+
                                         rigid_modes_indices.size());
  for (unsigned int i=0; i<2*water_line_indices.size()+3*bow_stern_indices.size()+water_indices.size()+phi_indices.size(); ++i)
      line_lengths[i] = 30;
  for (unsigned int i=0; i<rigid_modes_indices.size(); ++i)
      line_lengths[i+2*water_line_indices.size()+3*bow_stern_indices.size()+water_indices.size()+phi_indices.size()] = 
                   2*water_line_indices.size()+3*bow_stern_indices.size()+water_indices.size()+phi_indices.size()+rigid_modes_indices.size();

  jacobian_sparsity_pattern.reinit(2*water_line_indices.size()+3*bow_stern_indices.size()+water_indices.size()+phi_indices.size()+rigid_modes_indices.size(),
                                   2*water_line_indices.size()+3*bow_stern_indices.size()+water_indices.size()+phi_indices.size()+rigid_modes_indices.size(),
                 line_lengths );

  count = 0;
  for (std::set <unsigned int>::iterator pos = water_line_indices.begin(); pos != water_line_indices.end(); ++pos)
      {
      unsigned int i=*pos;
      jacobian_sparsity_pattern.add(count,count);
      jacobian_sparsity_pattern.add(count,count+1);
      count++;
      for (SparsityPattern::iterator col=free_surf_jacobian_dot.get_sparsity_pattern().begin(3*i+2);
           col!=free_surf_jacobian_dot.get_sparsity_pattern().end(3*i+2); ++col)
          if ( indices_map.count(col->column()) )
             jacobian_sparsity_pattern.add(count,indices_map.find(col->column())->second);
      count++;
      }
  for (std::set <unsigned int>::iterator pos = bow_stern_indices.begin(); pos != bow_stern_indices.end(); ++pos)
      {
      unsigned int i=*pos;
      Point<3> t(comp_dom.edges_tangents[3*i],comp_dom.edges_tangents[3*i+1],comp_dom.edges_tangents[3*i+2]);
      jacobian_sparsity_pattern.add(count,count);
      jacobian_sparsity_pattern.add(count,count+2);      
      count++;
      jacobian_sparsity_pattern.add(count,count);
      jacobian_sparsity_pattern.add(count,count+1);      
      count++;
      for (SparsityPattern::iterator col=free_surf_jacobian_dot.get_sparsity_pattern().begin(3*i+2);
           col!=free_surf_jacobian_dot.get_sparsity_pattern().end(3*i+2); ++col)
          if ( indices_map.count(col->column()) )
             jacobian_sparsity_pattern.add(count,indices_map.find(col->column())->second);
      count++;
      }  
  for (std::set <unsigned int>::iterator pos = water_indices.begin(); pos != water_indices.end(); ++pos)
      {
      unsigned int i=*pos;
      for (SparsityPattern::iterator col=free_surf_jacobian_dot.get_sparsity_pattern().begin(3*i+2);
           col!=free_surf_jacobian_dot.get_sparsity_pattern().end(3*i+2); ++col)
          if ( indices_map.count(col->column()) )
             jacobian_sparsity_pattern.add(count,indices_map.find(col->column())->second);
      count++;
      } 

  for (std::set <unsigned int>::iterator pos = phi_indices.begin(); pos != phi_indices.end(); ++pos)
      {
      unsigned int i=*pos;
      for (SparsityPattern::iterator col=free_surf_jacobian_dot.get_sparsity_pattern().begin(i+comp_dom.vector_dh.n_dofs());
           col!=free_surf_jacobian_dot.get_sparsity_pattern().end(i+comp_dom.vector_dh.n_dofs()); ++col)
          if ( indices_map.count(col->column()) )
             jacobian_sparsity_pattern.add(count,indices_map.find(col->column())->second);
      count++;
      }

  for (std::set <unsigned int>::iterator pos = rigid_modes_indices.begin(); pos != rigid_modes_indices.end(); ++pos)
      {
      unsigned int i=*pos;
      for (SparsityPattern::iterator col=free_surf_jacobian_dot.get_sparsity_pattern().begin(i+comp_dom.vector_dh.n_dofs());
           col!=free_surf_jacobian_dot.get_sparsity_pattern().end(i+comp_dom.vector_dh.n_dofs()); ++col)
          if ( indices_map.count(col->column()) )
             jacobian_sparsity_pattern.add(count,indices_map.find(col->column())->second);
      count++;
      }
  jacobian_sparsity_pattern.compress();
  jacobian_matrix.reinit(jacobian_sparsity_pattern);



  }

  virtual unsigned int n_dofs() const;

  virtual void output_step(Vector<double> & solution,
                           const unsigned int step_number);

				     /** For newton solver, we need a
					 residual function. This one is computed with
                                         sacado to allow computations of derivatives for jacobian*/
  virtual int residual(Vector<double> &dst,  
	               const Vector<double> &src_yy);

				     /** Jacobian vector product for newton solver. */
  virtual int jacobian(Vector<double> &dst,  
	               const Vector<double> &src_yy,
		       const Vector<double> &src);

				     /** Setup Jacobian preconditioner for Newton. */
  virtual int setup_jacobian_prec(const Vector<double> &src_yy);    
    
				     /** Jacobian inverse preconditioner
					 vector product for newton solver. */
  virtual int jacobian_prec(Vector<double> &dst,  
			    const Vector<double> &src_yy,
			    const Vector<double> &src);

				     /** Jacobian preconditioner
					 vector product for newton solver. */
  virtual int jacobian_prec_prod(Vector<double> &dst,  
			    const Vector<double> &src_yy,
			    const Vector<double> &src);

  private:

  FreeSurface<3> &free_surface;
  const NumericalTowingTank &comp_dom;
  const double &time;
  Vector<double> &free_surf_y;
  Vector<double> &free_surf_y_dot;
  const SparseMatrix<double> &free_surf_jacobian_dot;

  Vector<double> free_surf_jac_x_delta;
  Vector<double> free_surf_res;
  std::set<unsigned int> water_line_indices;
  std::set<unsigned int> bow_stern_indices;
  std::set<unsigned int> water_indices;
  std::set<unsigned int> phi_indices;
  std::set<unsigned int> rigid_modes_indices;
  SparsityPattern      jacobian_sparsity_pattern;
  SparseMatrix<double> jacobian_matrix;

  public: 
  std::map<unsigned int,unsigned int> indices_map;

};

#endif
