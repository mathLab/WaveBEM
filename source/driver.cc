#include "../include/driver.h"




template <int dim>
Driver<dim>::Driver(int argc, char **argv) :
                computational_domain(1,1),
		fma(computational_domain),
		bem_problem(computational_domain, fma),
	        free_surface(computational_domain, bem_problem),
	        dae_time_int(free_surface),
	        prm()
{

                                   // Declare the parameter entries..
  DeclareParameters();
    
  deallog.depth_console(3);

  std::list<string> args;
  for (int i=0; i<argc; ++i)
    args.push_back (argv[i]);

  string default_prm;
                                   // The default parameter file is the name of the application plus prm
  default_prm = args.back() + ".prm";
    
  prm.read_input(default_prm, false, true);

                                   // Now that we have the final version of the parameters, parse them.
  ParseParameters();
    
                                   // And write the used ones.
  default_prm = args.front() + "_used.prm";
  ofstream outprm(default_prm.c_str());
  prm.print_parameters(outprm, ParameterHandler::ShortText);
}

template <int dim>
Driver<dim>::~Driver()
{}


template <int dim>
void Driver<dim>::run()
{

 computational_domain.read_domain();
 Vector<double> y, yp, res, diff_compos,sys_comp;
 if(free_surface.initial_condition_from_dump) 
   {
     free_surface.restore_solution(y, yp,
				   free_surface.restore_filename);
     bem_problem.reinit();
     free_surface.reinit();
   }
 else
   {
     computational_domain.refine_and_resize();
     bem_problem.reinit(); 
     free_surface.reinit();
     y.reinit(free_surface.n_dofs());
     diff_compos.reinit(free_surface.n_dofs());
     diff_compos = free_surface.differential_components();
     sys_comp = free_surface.sys_comp;
     free_surface.initial_conditions(y);
     yp.reinit(free_surface.n_dofs());
     res.reinit(free_surface.n_dofs());
   }

 //free_surface.prepare_restart(0, y, yp);
 //computational_domain.generate_octree_blocking();

 dae_time_int.start_ode(y, yp, maxNumSteps);

 
 cout<<"Number of rhs evaluations is "<<free_surface.Rhs_evaluations_counter()<<endl; //*/

}

template <int dim>
void Driver<dim>::DeclareParameters()
{
  prm.declare_entry("Maximum number of steps", "1", Patterns::Integer());
  computational_domain.declare_parameters(prm);
  free_surface.declare_parameters(prm);
  dae_time_int.declare_parameters(prm);
  bem_problem.declare_parameters(prm);
  fma.declare_parameters(prm);    
}

template <int dim>
void Driver<dim>::ParseParameters()
{    
  maxNumSteps = prm.get_integer("Maximum number of steps");
  computational_domain.parse_parameters(prm);
  free_surface.parse_parameters(prm);
  dae_time_int.parse_parameters(prm);
  bem_problem.parse_parameters(prm);
  fma.parse_parameters(prm);
}


template class Driver<3>;
