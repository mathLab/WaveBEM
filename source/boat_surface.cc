//----------------------------  step-34.cc  ---------------------------
//    $Id: step-34.cc 18734 2009-04-25 13:36:48Z heltai $
//    Version: $Name$ 
//
//    Copyright (C) 2009, 2011 by the deal.II authors 
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//    Authors: Luca Heltai, Cataldo Manigrasso
//
//----------------------------  step-34.cc  ---------------------------


				 // @sect3{Include files}

				 // The program starts with including a bunch
				 // of include files that we will use in the
				 // various parts of the program. Most of them
				 // have been discussed in previous tutorials
				 // already:
				 

#include "../include/boat_surface.h"


template <int dim>
BoatSurface<dim>::BoatSurface()
{}

template <int dim>
void BoatSurface<dim>::declare_parameters(ParameterHandler &prm) {

  }

template <int dim>
void BoatSurface<dim>::parse_parameters(ParameterHandler &prm) {

  

}

template <int dim>
double BoatSurface<dim>::HullFunction(const Point<dim> point) const
{


    double L = 2.5; //wigley
    double B = L/10.0; //wigley
    double T = B/1.6; //wigley

    //double L = 2.5; // ellipse
    //double B = 0.5; // ellipse
    //double T = 1; //ellipse

    double x = point(0);
    double y = point(1);
    double z = point(2);
    
    double sign_y;
    
    if (fabs(y)/y > 0)
       sign_y = 1;
    else if (std::abs(y)/y < 0)
       sign_y = -1;
    else
       sign_y = 0;   
       
    if (std::abs(x) > L/2.0 || z < -T)
    	{
    	return 0.0;
    	}
    else
    	{
    	if (z > 0.0)
    		{
                //return sign_y*B/2.0*sqrt(1.0-pow(x/(L/2.0),2.0)); //ellipse
   		return sign_y*B/2.0*(1.0-pow(2.0*x/L,2.0)); //wigley
                //return 0; //tank
    		}
        else if (z == -T) //ellipse
                {         // ellipse
                return 0; // ellipse --->  return y;
                }         // ellipse
    	else
    		{
                //return sign_y*B/2.0*sqrt(1.0-pow(x/(L/2.0),2.0)); //ellipse
    		return sign_y*B/2.0*(1.0-pow(2.0*x/L,2.0))*(1.0-pow(z/T,2.0)); //wigley
                //return 0; //tank
    		}
    	}

}

template <int dim>
Point<dim> BoatSurface<dim>::HullNormal(const Point<dim> point) const
{


    double L = 2.5; // wigley
    double B = L/10.0; // wigley
    double T = B/1.6; // wigley
    //double L = 2.5; // ellipse
    //double B = 0.5; // ellipse
    //double T = 1; //ellipse

    double x = point(0);
    double y = point(1);
    double z = point(2);

    double sign_y;
    
    if (fabs(y)/y > 0)
       sign_y = 1;
    else if (fabs(y)/y < 0)
       sign_y = -1;
    else
       sign_y = 1;
    
    Point<dim> nHull;

    if (std::abs(x) > L/2.0 || z < -T)
    	{
    	return Point<dim>(0.0,1.0,0.0);
    	}
    else
    	{
    	if (z > 0.0)
    		{
    		nHull = Point<dim>(4.0*B/(L*L)*x,
					 sign_y,
					 0.0); // wigley
                //nHull = Point<dim>(-2*B*x/L/L/(sqrt(1.0-pow(x/(L/2.0),2.0))),
                //                        -sign_y,
                //                        0.0); // ellipse
    		nHull = nHull/nHull.distance(Point<dim>(0,0,0));

                //nHull = Point<dim>(0.0,1.0,0.0); //nothing 

   		return nHull;
    		}
        //else if (z == -T) //ellipse
        //        {         // ellipse
        //        //return Point<dim>(0,0,1); // ellipse
        //        return Point<dim>(0.0,1.0,0.0);
        //        }         // ellipse
    	else
    		{
    		nHull = Point<dim>(4.0*B/(L*L)*(1.0-pow(z/T,2.0))*x,
					 sign_y,
			          B/(T*T)*(1.0-pow(2.0*x/L,2.0))*z); // wigley
                //nHull = Point<dim>(-2*B*x/L/L/(sqrt(1.0-pow(x/(L/2.0),2.0))),
                //                         -sign_y,
                //                         0.0); // ellipse
    		nHull = nHull/nHull.distance(Point<dim>(0,0,0));
                //return Point<dim>(0.0,1.0,0.0); //nothing
    		return nHull;
    		}
    	}

}

template <int dim>
double BoatSurface<dim>::HullMeanCurvature(const Point<dim> point) const
{


    double L = 2.5; //wigley
    double B = L/10.0; //wigley
    double T = B/1.6; //wigley

    //double L = 2.5; // ellipse
    //double B = 0.5; // ellipse
    //double T = 1; //ellipse

    double x = point(0);
    double y = point(1);
    double z = point(2);
    
    double sign_y;
    
    if (fabs(y)/y > 0)
       sign_y = 1;
    else if (std::abs(y)/y < 0)
       sign_y = -1;
    else
       sign_y = 0;   
       
    if (std::abs(x) > L/2.0 || z < -T)
    	{
    	return 0.0;
    	}
    else
    	{
    	if (z > 0.0)
    		{
                double ds_dx = B/2.0*(-pow(2.0/L,2.0)*2*x);
                double dds_ddx = B/2.0*(-pow(2.0/L,2.0)*2.0);
                double ds_dz = 0.0;
                double dds_ddz = 0.0;
                double dds_dxdz = 0.0;
   		return sign_y*0.5*((1+ds_dx*ds_dx)*dds_ddz -2.0*ds_dx*ds_dz*dds_dxdz  + (1+ds_dz*ds_dz)*dds_ddx)/
                              pow(1+ds_dx*ds_dx+ds_dz*ds_dz,1.5);
                }
    	else
    		{
                double ds_dx = B/2.0*(-pow(2.0/L,2.0)*2*x)*(1.0-pow(z/T,2.0));
                double dds_ddx = B/2.0*(-pow(2.0/L,2.0)*2.0)*(1.0-pow(z/T,2.0));
                double ds_dz = B/2.0*(1.0-pow(2.0*x/L,2.0))*(-pow(1/T,2.0)*2.0*z);
                double dds_ddz = B/2.0*(1.0-pow(2.0*x/L,2.0))*(-pow(1/T,2.0)*2.0*z);
                double dds_dxdz = B/2.0*(-pow(2.0/L,2.0)*2*x)*(-pow(1/T,2.0)*2.0*z);
                
                return sign_y*0.5*((1+ds_dx*ds_dx)*dds_ddz -2.0*ds_dx*ds_dz*dds_dxdz  + (1+ds_dz*ds_dz)*dds_ddx)/
                              pow(1+ds_dx*ds_dx+ds_dz*ds_dz,1.5);
                }
    	}

}


template <int dim>
Point<dim> BoatSurface<dim>::get_new_point_on_line
(const typename Triangulation< dim-1,dim >::line_iterator &line) const
{
  Point<dim> newp = StraightBoundary<dim-1,dim>::get_new_point_on_line(line);
  newp(1) = HullFunction(newp);
  return newp;
}


template <int dim>
Point<dim> BoatSurface<dim>::get_new_point_on_quad
(const typename Triangulation< dim-1,dim >::quad_iterator &quad) const
{

  Point<dim> newp = StraightBoundary<dim-1,dim>::get_new_point_on_quad(quad);
  newp(1) = HullFunction(newp);
  return newp;
}

template class BoatSurface<3>;
