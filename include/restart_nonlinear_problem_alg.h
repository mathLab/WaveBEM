// The program starts with including a bunch
// of include files that we will use in the
// various parts of the program. Most of them
// have been discussed in previous tutorials
// already:
#ifndef restart_nonlinear_problem_alg_h
#define restart_nonlinear_problem_alg_h


#include <free_surface.h>
#include <numerical_towing_tank.h>
#include <newton_argument.h>


class RestartNonlinearProblemAlg:
  public NewtonArgument
{
public:
  RestartNonlinearProblemAlg(FreeSurface<3> &free_surface,
                             const NumericalTowingTank &comp_dom,
                             const double &time,
                             Vector<double> &free_surf_y,
                             Vector<double> &free_surf_y_dot,
                             const SparseMatrix<double> &free_surf_jacobian) :
    free_surface(free_surface),
    comp_dom(comp_dom),
    time(time),
    free_surf_y(free_surf_y),
    free_surf_y_dot(free_surf_y_dot),
    free_surf_jacobian(free_surf_jacobian)
  {

    dphi_dn_indices.clear();


    // let's take care of the potential normal derivative dofs
    for (unsigned int i=0; i<comp_dom.dh.n_dofs(); ++i)
      {
        // if node is not on water all dofs are algebraic components
        if (!(comp_dom.flags[i] & water))
          {
            // neumann boundary condition on boat surface is what we're interested in
            //if (comp_dom.vector_constraints.is_constrained(3*i) )
            //   {// if node is constrained, it is disregarded
            //   }
            //else
            //   {// otherwise, it is the kind of dof we're looking for
            dphi_dn_indices.insert(i);
            //   }
          }
      }

    free_surf_jac_x_delta.reinit(free_surf_y.size());
    free_surf_res.reinit(free_surf_y.size());



    unsigned int count = 0;

    for (std::set <unsigned int>::iterator pos = dphi_dn_indices.begin(); pos != dphi_dn_indices.end(); ++pos)
      {
        unsigned int i=*pos;
        indices_map[i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()] = count;
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
    //cout<<dphi_dn_indices.size()<<endl;
    //cout<<2*water_line_indices.size()+3*bow_stern_indices.size()+water_indices.size()+phi_indices.size()+dphi_dn_indices.size()<<endl;
    jacobian_sparsity_pattern.reinit(dphi_dn_indices.size(),
                                     dphi_dn_indices.size(),
                                     30);

    count = 0;
    for (std::set <unsigned int>::iterator pos = dphi_dn_indices.begin(); pos != dphi_dn_indices.end(); ++pos)
      {
        unsigned int i=*pos;
        for (SparsityPattern::iterator col=free_surf_jacobian.get_sparsity_pattern().begin(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs());
             col!=free_surf_jacobian.get_sparsity_pattern().end(i+comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()); ++col)
          if ( (indices_map.count(col->column())) && ((col->column())>= comp_dom.vector_dh.n_dofs()+comp_dom.dh.n_dofs()) )
            {
              jacobian_sparsity_pattern.add(count,indices_map.find(col->column())->second);
              //cout<<"prima "<<count<<","<<indices_map.find(col->column())->second<<endl;
            }
        count++;
      }

    jacobian_sparsity_pattern.compress();
    jacobian_matrix.reinit(jacobian_sparsity_pattern);



  }

  virtual unsigned int n_dofs() const;

  virtual void output_step(Vector<double> &solution,
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
  const SparseMatrix<double> &free_surf_jacobian;

  Vector<double> free_surf_jac_x_delta;
  Vector<double> free_surf_res;
  std::set<unsigned int> dphi_dn_indices;
  SparsityPattern      jacobian_sparsity_pattern;
  SparseMatrix<double> jacobian_matrix;

public:
  std::map<unsigned int,unsigned int> indices_map;

};

#endif
