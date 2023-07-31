#include "surface_smoothing.h"

#include <deal.II/base/logstream.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/numerics/matrix_tools.h>

// all include files you need here

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/solution_transfer.h>

#include "occ_utilities.h"

using namespace dealii;
using namespace std;
using namespace OpenCascade;

SurfaceSmoothing::SurfaceSmoothing(Vector<double> &euler_vector,
                                   Vector<double> &curvature_vector,
                                   const DoFHandler<2,3> &dh,
                                   const Mapping<2,3> &mapping) :
  euler_vector(euler_vector),
  curvature_vector(curvature_vector),
  dh(dh),
  mapping(mapping)
{
  update_reference();
}

void SurfaceSmoothing::update_reference()
{
  Assert(dh.n_dofs() == euler_vector.size(),
         ExcDimensionMismatch(dh.n_dofs(), euler_vector.size()));

  Assert(dh.n_dofs() == curvature_vector.size(),
         ExcDimensionMismatch(dh.n_dofs(), curvature_vector.size()));

  reference_identity.reinit(dh.n_dofs());
  solution.reinit(dh.n_dofs());
  rhs.reinit(dh.n_dofs());

  matrix.clear();
  mass_matrix.clear();

  constraints.clear();
  DoFTools::make_hanging_node_constraints (dh,constraints);
  constraints.close();

  DynamicSparsityPattern csp (dh.n_dofs(), dh.n_dofs());
  DoFTools::make_sparsity_pattern (dh, csp, constraints);
  sparsity.copy_from (csp);

  matrix.reinit(sparsity);
  mass_matrix.reinit(sparsity);

  vector<Point<3> > ref_support_points(dh.n_dofs());
  // Compute the reference support
  // points.
  DoFTools::map_dofs_to_support_points<2,3>(StaticMappingQ1<2,3>::mapping,
                                            dh, ref_support_points);

  // build the reference_identity
  // vector
  for (unsigned int i=0; i<dh.n_dofs(); ++i)
    reference_identity(i) = ref_support_points[i](i%3);

  std::vector< bool > comp_sel(3, true);
  boundary_dofs.resize(dh.n_dofs());

  DoFTools::extract_boundary_dofs(dh, comp_sel, boundary_dofs);
}

void SurfaceSmoothing::smooth()
{
  Assert(euler_vector.size() == reference_identity.size(),
         ExcDimensionMismatch(euler_vector.size(), reference_identity.size()));

  fix_boundary_values();
//cout<<"S1"<<endl;
  assemble_system();
//cout<<"S2"<<endl;
  solve_system();
//cout<<"S3"<<endl;
  euler_vector.sadd(0., 1., solution);
  euler_vector.sadd(1.,-1., reference_identity);
//cout<<"S4"<<endl;
}

void SurfaceSmoothing::fix_boundary_values()
{
  boundary_values.clear();

  for (unsigned int i=0; i<dh.n_dofs(); ++i)
    if (boundary_dofs[i] == true)
      {
        boundary_values[i] = reference_identity(i)+euler_vector(i);
        curvature_vector(i) = 0;
      }
}

void SurfaceSmoothing::assemble_system()
{
  assemble_system(curvature_vector);
}


void SurfaceSmoothing::assemble_system(const Vector<double> &curvature)
{
  solution = 0;
  rhs = 0;
  matrix = 0;
  //cout<<"SS1"<<endl;
  const FiniteElement<2,3> &fe = dh.get_fe();
  QGauss<2> quad(fe.degree*2+1);

  FEValues<2,3> fe_v(mapping, fe, quad,
                     update_values |
                     update_gradients |
                     update_JxW_values);
  //cout<<"SS2"<<endl;
  const unsigned int n_q_points = fe_v.n_quadrature_points;
  const unsigned int dofs_per_cell   = fe.dofs_per_cell;
  std::vector<unsigned int> local_dof_indices (dofs_per_cell);
  std::vector<Vector<double> > local_curvature (n_q_points, Vector<double>(3));

  FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
  FullMatrix<double>   local_mass_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       local_rhs (dofs_per_cell);
  //cout<<"SS3"<<endl;
  DoFHandler<2,3>::active_cell_iterator
  cell = dh.begin_active(),
  endc = dh.end();
  //cout<<"SS4"<<endl;
  for (; cell!=endc; ++cell)
    {
      fe_v.reinit (cell);
      local_matrix = 0;
      local_mass_matrix = 0;
      local_rhs = 0;
      fe_v.get_function_values(curvature,
                               local_curvature);

      unsigned int comp_i, comp_j;

      for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          comp_i = fe.system_to_component_index(i).first;
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            {
              comp_j = fe.system_to_component_index(j).first;
              if (comp_i == comp_j)

                for (unsigned int q=0; q<n_q_points; ++q)
                  {
                    local_matrix(i,j) += fe_v.shape_grad(i,q)*
                                         fe_v.shape_grad(j,q)*
                                         fe_v.JxW(q);
                    local_mass_matrix(i,j) += fe_v.shape_value(i,q)*
                                              fe_v.shape_value(j,q)*
                                              fe_v.JxW(q);

                  }
            }
          for (unsigned int q=0; q<n_q_points; ++q)
            {
              local_rhs(i) += fe_v.shape_value(i,q)*
                              local_curvature[q](comp_i)*
                              fe_v.JxW(q);
            }
        }
      //cout<<"SS5"<<endl;
      cell->get_dof_indices (local_dof_indices);

      constraints.distribute_local_to_global
      (local_matrix, local_rhs, local_dof_indices, matrix, rhs);

      constraints.distribute_local_to_global
      (local_mass_matrix, local_dof_indices, mass_matrix);
      //cout<<"SS6"<<endl;
    }
  MatrixTools::apply_boundary_values(boundary_values, matrix, solution, rhs);
  //cout<<"SS7"<<endl;
}

void SurfaceSmoothing::solve_system()
{
  SparseDirectUMFPACK inverse;
  inverse.initialize(matrix);
  inverse.vmult(solution, rhs);
  constraints.distribute(solution);
}


void SurfaceSmoothing::compute_curvatures(Vector<double> &curvatures)
{
  Assert(curvatures.size() == dh.n_dofs(),
         ExcDimensionMismatch(curvatures.size(), dh.n_dofs()));
  Vector <double> positions(dh.n_dofs());
  for (unsigned int i=0; i<dh.n_dofs(); ++i)
    {
      positions(i) = reference_identity(i)+euler_vector(i);
    }
  assemble_system();
  matrix.vmult(curvatures,positions);
  SparseDirectUMFPACK inverse_mass;
  inverse_mass.initialize(mass_matrix);
  inverse_mass.solve(curvatures);
  constraints.distribute(curvatures);
}

void SurfaceSmoothing::apply_curvatures(const Vector<double> &curvatures,
                                        const vector<bool> &boundary_dofs)
{
  AssertThrow(curvatures.size() == dh.n_dofs(),
              ExcDimensionMismatch(curvatures.size(), dh.n_dofs()));
  //fix_boundary_values();

  boundary_values.clear();

  for (unsigned int i=0; i<dh.n_dofs(); ++i)
    if (boundary_dofs[i] == true)
      {
        boundary_values[i] = reference_identity(i)+euler_vector(i);
        curvature_vector(i) = 0;
      }

  assemble_system(curvatures);
  solve_system();
  euler_vector.sadd(0., 1., solution);
  euler_vector.sadd(1.,-1., reference_identity);
}

