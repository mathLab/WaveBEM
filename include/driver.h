#ifndef driver_h
#define driver_h


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

#include "bem_problem.h"
#include "bem_fma.h"
#include "free_surface.h"
#include "numerical_towing_tank.h"

#include <dae_time_integrator.h>

using namespace dealii;


template <int dim>
class Driver 
{
  public:
   
   Driver(int argc, char **argv);
   
   ~Driver();

   void run();
    
   void DeclareParameters(); 

   void ParseParameters();

   static string get_library_names(); 
   
   private:
   
   NumericalTowingTank computational_domain;
   

   BEMFMA<dim> fma;

   BEMProblem<dim> bem_problem;
      
   FreeSurface<dim> free_surface;

   DAETimeIntegrator dae_time_int;
   
   ParameterHandler prm;

   string library_str;

   unsigned int maxNumSteps;

};

#endif

