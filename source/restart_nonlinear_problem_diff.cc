#include "../include/restart_nonlinear_problem_diff.h"

unsigned int RestartNonlinearProblemDiff::n_dofs() const
{

return 2*water_line_indices.size()+3*bow_stern_indices.size()+water_indices.size()+phi_indices.size()+rigid_modes_indices.size();
}



void RestartNonlinearProblemDiff::output_step(Vector<double> & solution,
                 const unsigned int step_number)
{
}


				     /** For newton solver, we need a
					 residual function. This one is computed with
                                         sacado to allow computations of derivatives for jacobian*/
int RestartNonlinearProblemDiff::residual(Vector<double> &dst,  
	                                      const Vector<double> &src_yy)
{

  // we first put all the "proposed" velocity values in the free surface y_dot vector (with proper positions)
  for (std::map <unsigned int, unsigned int>::iterator pos = indices_map.begin(); pos != indices_map.end(); ++pos)
      {
      free_surf_y_dot(pos->first) = src_yy(pos->second);
      }



  // now we can call the residual funtion of free surface class
  free_surf_jac_x_delta=0;
  free_surface.residual_and_jacobian(time,free_surf_res,free_surf_y,free_surf_y_dot,free_surf_jac_x_delta,1.0,false);


  // the resulting free surface residual must be "distributed" on the present problem degrees of freedom (with proper ordering)      
  cout<<"Restart problem residual "<<endl;
  unsigned int count = 0;
  for (std::set <unsigned int>::iterator pos = water_line_indices.begin(); pos != water_line_indices.end(); ++pos)
      {
      unsigned int i=*pos;
      dst(count) = src_yy(count) + src_yy(count+1)*comp_dom.iges_normals[i](2)/comp_dom.iges_normals[i](1);
      //cout<<count<<"("<<3*i+1<<")WL  "<<dst(count)<<"    "<<src_yy(count)<<endl;
      count++;
      dst(count) = free_surf_res(3*i+2);
      //cout<<count<<"("<<3*i+2<<")WL  "<<dst(count)<<"    "<<src_yy(count)<<endl;
      count++;
      }
  for (std::set <unsigned int>::iterator pos = bow_stern_indices.begin(); pos != bow_stern_indices.end(); ++pos)
      {
      unsigned int i=*pos;
      Point<3> t(comp_dom.edges_tangents[3*i],comp_dom.edges_tangents[3*i+1],comp_dom.edges_tangents[3*i+2]);
      dst(count) = src_yy(count) - src_yy(count+2)*t(0)/t(2);
      //cout<<count<<"("<<3*i<<")BS  "<<dst(count)<<"    "<<src_yy(count)<<endl;
      count++;
      dst(count) = src_yy(count) - src_yy(count+1)*t(1)/t(2);
      //cout<<count<<"("<<3*i+1<<")BS  "<<dst(count)<<"    "<<src_yy(count)<<endl;
      count++;
      dst(count) = free_surf_res(3*i+2);
      //cout<<count<<"("<<3*i+2<<")BS  "<<dst(count)<<"    "<<src_yy(count)<<endl;
      count++;
      }  
  for (std::set <unsigned int>::iterator pos = water_indices.begin(); pos != water_indices.end(); ++pos)
      {
      unsigned int i=*pos;
      dst(count) = free_surf_res(3*i+2);
      //cout<<count<<"("<<3*i+2<<")W  "<<dst(count)<<"    "<<src_yy(count)<<endl;
      count++;
      } 
  for (std::set <unsigned int>::iterator pos = phi_indices.begin(); pos != phi_indices.end(); ++pos)
      {
      unsigned int i=*pos;
      dst(count) = free_surf_res(i+comp_dom.vector_dh.n_dofs());
      //cout<<count<<"("<<3*i+2<<")W  "<<dst(count)<<"    "<<src_yy(count)<<endl;
      count++;
      }
  for (std::set <unsigned int>::iterator pos = rigid_modes_indices.begin(); pos != rigid_modes_indices.end(); ++pos)
      {
      unsigned int i=*pos;
      dst(count) = free_surf_res(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs());
      //cout<<count<<"("<<3*i+2<<")W  "<<dst(count)<<"    "<<src_yy(count)<<endl;
      count++;
      }

  cout<<"Restart problem nonlin residual: "<<dst.l2_norm()<<endl;
  return 0;
}

				     /** Jacobian vector product for newton solver. */
int RestartNonlinearProblemDiff::jacobian(Vector<double> &dst,  
	                              const Vector<double> &src_yy,
		                      const Vector<double> &src)
{


  //cout<<"Original jacobian"<<endl;
  //for (unsigned int i=0; i<free_surf_res.size(); ++i)
  //    for (SparsityPattern::iterator c=free_surf_jacobian.get_sparsity_pattern().begin(i); c!=free_surf_jacobian.get_sparsity_pattern().end(i); ++c)
  //        {
  //        unsigned int j = c->column();
  //        cout<<" "<<i<<" "<<j<<" "<<free_surf_jacobian(i,j)<<endl;
  //        }
  cout<<"Restart problem jacobian"<<endl;
  jacobian_matrix = 0;

  unsigned int count = 0;
  for (std::set <unsigned int>::iterator pos = water_line_indices.begin(); pos != water_line_indices.end(); ++pos)
      {
      unsigned int i=*pos;
      jacobian_matrix.set(count,count,1.0);
      jacobian_matrix.set(count,count+1,comp_dom.iges_normals[i](2)/comp_dom.iges_normals[i](1));
      //cout<<count<<" "<<count<<" "<<1.0<<endl;
      //cout<<count<<" "<<count+1<<" "<<comp_dom.iges_normals[i](2)/comp_dom.iges_normals[i](1)<<endl;
      count++;
      for (SparsityPattern::iterator col=free_surf_jacobian_dot.get_sparsity_pattern().begin(3*i+2);
           col!=free_surf_jacobian_dot.get_sparsity_pattern().end(3*i+2); ++col)
          if ( indices_map.count(col->column()) )
             {
             jacobian_matrix.set(count,indices_map.find(col->column())->second,free_surf_jacobian_dot(3*i+2,col->column()));
             //cout<<count<<" "<<indices_map.find(col->column())->second<<" "<<free_surf_jacobian_dot(3*i+2,col->column())<<endl;
             }
      count++;
      }
  for (std::set <unsigned int>::iterator pos = bow_stern_indices.begin(); pos != bow_stern_indices.end(); ++pos)
      {
      unsigned int i=*pos;
      Point<3> t(comp_dom.edges_tangents[3*i],comp_dom.edges_tangents[3*i+1],comp_dom.edges_tangents[3*i+2]);
      jacobian_matrix.set(count,count,1.0);
      jacobian_matrix.set(count,count+2,-t(0)/t(2));
      //cout<<count<<" "<<count<<" "<<1.0<<endl;
      //cout<<count<<" "<<count+2<<" "<<-t(0)/t(2)<<endl;
      count++;
      jacobian_matrix.set(count,count,1.0);
      jacobian_matrix.set(count,count+1,-t(1)/t(2));
      //cout<<count<<" "<<count<<" "<<1.0<<endl;
      //cout<<count<<" "<<count+1<<" "<<-t(1)/t(2)<<endl;
      count++;
      for (SparsityPattern::iterator col=free_surf_jacobian_dot.get_sparsity_pattern().begin(3*i+2);
           col!=free_surf_jacobian_dot.get_sparsity_pattern().end(3*i+2); ++col)
          if ( indices_map.count(col->column()) )
             {
             jacobian_matrix.set(count,indices_map.find(col->column())->second,free_surf_jacobian_dot(3*i+2,col->column()));
             //cout<<count<<" "<<indices_map.find(col->column())->second<<" "<<free_surf_jacobian_dot(3*i+2,col->column())<<endl;
             }
      count++;
      }  
  for (std::set <unsigned int>::iterator pos = water_indices.begin(); pos != water_indices.end(); ++pos)
      {
      unsigned int i=*pos;
      for (SparsityPattern::iterator col=free_surf_jacobian_dot.get_sparsity_pattern().begin(3*i+2);
           col!=free_surf_jacobian_dot.get_sparsity_pattern().end(3*i+2); ++col)
          if ( indices_map.count(col->column()) )
             {
             jacobian_matrix.set(count,indices_map.find(col->column())->second,free_surf_jacobian_dot(3*i+2,col->column()));
             //cout<<count<<" "<<indices_map.find(col->column())->second<<" "<<free_surf_jacobian_dot(3*i+2,col->column())<<endl;
             }
      count++;
      } 
  for (std::set <unsigned int>::iterator pos = phi_indices.begin(); pos != phi_indices.end(); ++pos)
      {
      unsigned int i=*pos;
      for (SparsityPattern::iterator col=free_surf_jacobian_dot.get_sparsity_pattern().begin(i+comp_dom.vector_dh.n_dofs());
           col!=free_surf_jacobian_dot.get_sparsity_pattern().end(i+comp_dom.vector_dh.n_dofs()); ++col)
          if ( indices_map.count(col->column()) )
             {
             jacobian_matrix.set(count,indices_map.find(col->column())->second,free_surf_jacobian_dot(i+comp_dom.vector_dh.n_dofs(),col->column()));
             //cout<<count<<" "<<indices_map.find(col->column())->second<<" "<<free_surf_jacobian_dot(i+comp_dom.vector_dh.n_dofs(),col->column())<<endl;
             }
      count++;
      }

  for (std::set <unsigned int>::iterator pos = rigid_modes_indices.begin(); pos != rigid_modes_indices.end(); ++pos)
      {
      unsigned int i=*pos;
      for (SparsityPattern::iterator col=free_surf_jacobian_dot.get_sparsity_pattern().begin(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs());
           col!=free_surf_jacobian_dot.get_sparsity_pattern().end(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs()); ++col)
          if ( indices_map.count(col->column()) )
             {
             jacobian_matrix.set(count,indices_map.find(col->column())->second,free_surf_jacobian_dot(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),col->column()));
             //cout<<count<<" "<<indices_map.find(col->column())->second<<" "<<free_surf_jacobian_dot(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()+comp_dom.dh.n_dofs(),col->column())<<endl;
             }
      count++;
      }

  jacobian_matrix.vmult(dst,src);

  return 0;
}


				     /** Inverse preconditioner vector product for newton solver. */
int RestartNonlinearProblemDiff::jacobian_prec(Vector<double> &dst,  
	                                   const Vector<double> &src_yy,
		                           const Vector<double> &src)
{

  SparseDirectUMFPACK prec_direct;
  prec_direct.initialize(jacobian_matrix);
  prec_direct.vmult(dst,src);

  return 0;
}

				     /** Preconditioner vector product for newton solver. */
int RestartNonlinearProblemDiff::jacobian_prec_prod(Vector<double> &dst,  
	                                        const Vector<double> &src_yy,
		                                const Vector<double> &src)
{
  jacobian_matrix.vmult(dst,src);

  return 0;
}

				     /** Preconditioner is the jacobian itself, so this matrix has nothing to do */
int RestartNonlinearProblemDiff::setup_jacobian_prec(const Vector<double> &src_yy)
{
  return 0;
}
