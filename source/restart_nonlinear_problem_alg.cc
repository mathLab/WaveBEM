#include "../include/restart_nonlinear_problem_alg.h"

unsigned int RestartNonlinearProblemAlg::n_dofs() const
{

  return dphi_dn_indices.size();
}



void RestartNonlinearProblemAlg::output_step(Vector<double> &solution,
                                             const unsigned int step_number)
{
}


/** For newton solver, we need a
residual function. This one is computed with
                            sacado to allow computations of derivatives for jacobian*/
int RestartNonlinearProblemAlg::residual(Vector<double> &dst,
                                         const Vector<double> &src_yy)
{

  // we first put all the "proposed" dphi_dn values in the free surface y_dot vector (with proper positions)
  for (std::map <unsigned int, unsigned int>::iterator pos = indices_map.begin(); pos != indices_map.end(); ++pos)
    {
      free_surf_y(pos->first) = src_yy(pos->second);
    }

  // now we can call the residual funtion of free surface class
  free_surf_jac_x_delta=0;
  free_surface.residual_and_jacobian(time,free_surf_res,free_surf_y,free_surf_y_dot,free_surf_jac_x_delta,1.0,false);


  // the resulting free surface residual must be "distributed" on the present problem degrees of freedom (with proper ordering)
  cout<<"Restart problem residual "<<endl;
  unsigned int count = 0;
  for (std::set <unsigned int>::iterator pos = dphi_dn_indices.begin(); pos != dphi_dn_indices.end(); ++pos)
    {
      unsigned int i=*pos;
      dst(count) = free_surf_res(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs());
      //cout<<count<<"("<<3*i+2<<")W  "<<dst(count)<<"    "<<src_yy(count)<<endl;
      count++;
    }

  cout<<"Restart problem nonlin residual: "<<dst.l2_norm()<<endl;
  return 0;
}

/** Jacobian vector product for newton solver. */
int RestartNonlinearProblemAlg::jacobian(Vector<double> &dst,
                                         const Vector<double> &src_yy,
                                         const Vector<double> &src)
{


  //cout<<"Original jacobian"<<endl;
  //for (unsigned int i=0; i<free_surf_res.size(); ++i)
  //    for (SparsityPattern::iterator c=free_surf_jacobian.get_sparsity_pattern().row_begin(i); c!=free_surf_jacobian.get_sparsity_pattern().row_end(i); ++c)
  //        {
  //        unsigned int j = *c;
  //        cout<<" "<<i<<" "<<j<<" "<<free_surf_jacobian(i,j)<<endl;
  //        }
  cout<<"Restart problem jacobian"<<endl;
  jacobian_matrix = 0;

  unsigned int count = 0;
  for (std::set <unsigned int>::iterator pos = dphi_dn_indices.begin(); pos != dphi_dn_indices.end(); ++pos)
    {
      unsigned int i=*pos;
      for (SparsityPattern::iterator col=free_surf_jacobian.get_sparsity_pattern().begin(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs());
           col!=free_surf_jacobian.get_sparsity_pattern().end(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()); ++col)
        if ( (indices_map.count(col->column())) && ((col->column()) > (comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()-1)) )
          {
            jacobian_matrix.set(count,indices_map.find(col->column())->second,free_surf_jacobian(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs(),col->column()));
            //cout<<"dopo "<<count<<","<<indices_map.find(col->column())->second<<endl;
            //cout<<count<<" "<<indices_map.find(col->column())->second<<" "<<free_surf_jacobian(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs(),col->column())<<endl;
          }
      count++;
    }
  jacobian_matrix.vmult(dst,src);

  return 0;
}


/** Inverse preconditioner vector product for newton solver. */
int RestartNonlinearProblemAlg::jacobian_prec(Vector<double> &dst,
                                              const Vector<double> &src_yy,
                                              const Vector<double> &src)
{
  SparseDirectUMFPACK prec_direct;
  prec_direct.initialize(jacobian_matrix);
  prec_direct.vmult(dst,src);

  return 0;
}

/** Preconditioner vector product for newton solver. */
int RestartNonlinearProblemAlg::jacobian_prec_prod(Vector<double> &dst,
                                                   const Vector<double> &src_yy,
                                                   const Vector<double> &src)
{
  jacobian_matrix.vmult(dst,src);

  return 0;
}

/** Preconditioner is the jacobian itself, so this matrix has nothing to do */
int RestartNonlinearProblemAlg::setup_jacobian_prec(const Vector<double> &src_yy)
{
  return 0;
}
