#ifndef surface_smoothing_h
#define surface_smoothing_h

#include <string>

#include <deal.II/base/point.h>
#include <deal.II/grid/tria.h>

#include <deal.II/fe/mapping.h>

#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/matrix_tools.h>

using namespace dealii;

class SurfaceSmoothing
{
public:
  //! Smooth all dofs in
  // euler_vector, by maintaining
  // the boundary dofs at their
  // position.
  //
  // The euler_vector is
  // interpreted as a collection
  // of triples associated with
  // the displacement of the
  // nodes of the mesh, like the
  // one which is used in
  // MappingQEulerian. Mapping
  // should be parent of
  // MappingQEulerian, in some
  // sense. It is used to compute
  // all integrals, which means
  // it is required to return
  // points in real space.
  //
  // Some vectors are stored
  // internally to allow for
  // faster computations between
  // unchanged reference
  // meshes. If you change the
  // mesh, you should call
  // update_reference() before
  // attempting to call smooth().
  //
  // You should also call
  // update_reference() the first
  // time you use this class,
  // otherwise exceptions will be
  // thrown.
  //
  // By default, the smoothing is
  // performed only on the euler
  // vector. If you want the
  // smoothing on the global
  // position of the nodes, set
  // on_euler to false.
  SurfaceSmoothing(Vector<double> &euler_vector,
                   Vector<double> &curvature_vector,
                   const DoFHandler<2,3> &dh,
                   const Mapping<2,3> &mapping);

  //! Whenever the underlying dh
  // changes, you should call
  // this function to update the
  // internal vectors.
  void update_reference();

  //! All boundaries are fixed to
  // the current value of the
  // stored identity plus the
  // given euler_vector.
  void fix_boundary_values();
  //! Assemble the Laplace
  // Beltrami matrix, rhs and
  // solution, ready to be
  // solved. Boundary values are
  // applied, so be sure to call
  // fix_boundary_values, before
  // you call this function.
  void assemble_system();

  void assemble_system(const Vector<double> &curvature);

  //! Solve the system.
  void solve_system();

  //! Perform the full smoothing
  // cycle of the given mesh.
  void smooth();

  //! Compute curvatures at the
  // mesh nodes (to be used before
  // refinements)
  void compute_curvatures(Vector<double> &curvatures);
  //! Apply curvatures at the
  // mesh nodes to obtain geometrically
  // consistent mesh (to be used after
  // refinements)
  void apply_curvatures(const Vector<double> &curvatures, const std::vector<bool> &boundary_dofs);



private:
  Vector<double> &euler_vector;
  Vector<double> &curvature_vector;
  const DoFHandler<2,3> &dh;
  const Mapping<2,3> &mapping;

  Vector<double> reference_identity;
  std::vector<bool> boundary_dofs;
  std::map<unsigned int, double> boundary_values;

  Vector<double> solution;
  Vector<double> rhs;

  ConstraintMatrix constraints;

  SparsityPattern sparsity;
  SparseMatrix<double> matrix;
  SparseMatrix<double> mass_matrix;
};


#endif
