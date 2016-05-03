#ifndef lalpace_kernel_h
#define lalpace_kernel_h
#include <deal.II/base/point.h>

#include <cmath>

namespace LaplaceKernel
{
  template <int dim>
  double single_layer(const Point<dim> &R)
  {
    switch(dim)
      {
      case 2:
	return (-std::log(R.norm()) / (2*numbers::PI) );
	
      case 3:
	return (1./( R.norm()*4*numbers::PI ) );
	
      default:
	Assert(false, ExcInternalError());
	return 0.;
      }
  }
        


  template <int dim>
  Point<dim> double_layer(const Point<dim> &R)
  {
    switch(dim)
      {
      case 2:
	return R / ( -2*numbers::PI * R.square());
      case 3:
	return R / ( -4*numbers::PI * R.square() * R.norm() );
	
      default:
	Assert(false, ExcInternalError());
	return Point<dim>();
      }
  }
  
    template <int dim>
  void kernels(const Point<dim> &R, Point<dim> &D, double &d)
  {
    double r = R.norm();
    double r2 = r*r;
    switch(dim)
      {
      case 2:
	d = -std::log(r) / (2*numbers::PI);
	D = R / ( -2*numbers::PI * r2);
    break;
      case 3:
	d = (1./( r*4*numbers::PI ) );
	D = R / ( -4*numbers::PI * r2*r );
    break;
      default:
	Assert(false, ExcInternalError());
      }
  }
  
}
#endif
