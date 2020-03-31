#pragma once
#include "angem/VTK_ID.hpp"
#include "mesh/Mesh.hpp"
#include "angem/Tensor2.hpp"

namespace discretization {

using angem::VTK_ID;
using Point = angem::Point<3,double>;
template<VTK_ID id> constexpr size_t N_ELEMENT_VERTICES;
template<VTK_ID id> constexpr size_t ELEMENT_DIM;

/**
 * This is an abstract class that combines methods for standard finite elements.
 * It can be used to compute the values of shape functions, their gradients,
 * and JxW values.
 */
template<VTK_ID vtk_id>
class FeValues
{
 public:
  /**
   * Constructor.
   * This is an abstract class so the constuctor is protected.
   */
  FeValues(const mesh::Mesh & grid);
  /**
   * Update the needed internal quantities in the cell.
   * Use this method before using the fem quantities in the new cell during the loop.
   */
  void update(const mesh::Cell & cell);
  /**
   * Update the needed internal quantities in the cell in the gived integration point.
   * The point coordinate is the real (not reference) coordinates of the integration point.
   */
  void update(const mesh::Cell & cell, const angem::Point<3,double> & point);
  /**
   * Update the needed internal quantities in the face element.
   * Use this method before using the fem quantities in the new face during the loop.
   */
  void update(const mesh::Face & face);
  /**
   * Update the needed quantities in the face in the gived integration point.
   * The point coordinate is the real (not reference) coordinates of the integration point.
   */
  void update(const mesh::Face & cell, const angem::Point<3,double> & point);
  /**
   * Computes a vector of integration points in the master element
   */
  std::vector<Point> get_master_integration_points() const;
  /**
   * Computes a vector of integration point weights in the master element
   */
  std::vector<double> get_master_integration_weights() const;

  /**
   * Return the value of the shape function indexed by shape_index in the
   * qpoint integration point in the current element.
   * NOTE: must call update() before calling this function.
   */
  double value(const size_t shape_index, const size_t qpoint) const;
  /**
   * Return the gradient of the shape function indexed by shape_index in the
   * qpoint integration point in the current element.
   * NOTE: must call update() before calling this function.
   */
  Point grad(const size_t shape_index, const size_t qpoint) const;
  /**
   * Return the JxW value at the integration point indexed by qpoint.
   * JxW is a product of the determinant of the transformation jacobian by
   * the weight of the integration point.
   * NOTE: must call update() before calling this function.
   */
  double JxW(const size_t qpoint) const;

  /**
   * Number of integration points in the quadrature rule.
   */
  size_t n_integration_points() const { return 1; }

 protected:
  /**
   * Given element vertices and integration points, update shape values,
   * gradients, and the determinant of the trasformation jacobian.
   * This function is for cells (3D elements).
   */
  void update_();

  /**
   * Given element vertices and integration points, update shape values,
   * gradients, and the determinant of the trasformation jacobian.
   * Same as above but for faces (2d elemenet).
   */
  void update_face_();
  /**
   * Given a vector of the grid vertex indices, fill out the interal
   * _vertex_coord array.
   */
  void update_vertex_coord_(const std::vector<size_t> & vertex_indices);
  /**
   * Evalues the values of the ith shape function in point in the reference element
   * \param[in] point  : coordinates in the reference element
   * \param[in] vertex : the local number of the shape function
   * returns the value of the given shape function in a given point
   */
  double eval_(const Point & point, const size_t vertex) const;
  /**
   * Evalues the derivatives of the ith shape function in point in the reference element
   * \param[in] point  : coordinates in the reference element
   * \param[in] vertex : the local number of the shape function
   * returns the derivative of the given shape function in a given point
   */
   Point eval_derivative_(const Point & point, const size_t vertex) const;

  /**
   * Compute transformation jacobian matrix dx / dx_ref
   * Input:
   * \param[in] ref_grad : shape function grad in current q-point in master element
   * Returns transposed matrix dx / dx_ref
   */
   angem::Tensor2<3, double> compute_jacobian_(const std::array<Point,N_ELEMENT_VERTICES<vtk_id>> & ref_grad) const;

  // compute shape values in q-points in master element
  void update_shape_values_();
  // compute shape grads in q-points in master element
  void update_shape_grads_(const std::array<Point,N_ELEMENT_VERTICES<vtk_id>> & ref_grad,
                           const angem::Tensor2<3, double> & du_dx,
                           const size_t integration_point);
  // compute the gradients of shape fucntions in qpoint in master element
  std::array<Point,N_ELEMENT_VERTICES<vtk_id>> compute_ref_gradient_(const size_t qpoint) const;

  /**
   * Map from real coordinates xyz to local coordinates in the master element uvw
   */
  Point map_real_to_local_(const Point & p) const;

  /**
   * Compute face jacobian, its determinant, and invert it.
   * Need a separate method for it since the jacobian will be 2x2 instead of 3x3.
   */
  // std::tuple<angem::Tensor2<3, double>>
  std::pair<angem::Tensor2<3, double>,double>
  compute_detJ_and_invert_face_jacobian_(const std::array<Point,N_ELEMENT_VERTICES<vtk_id>> & ref_grad) const;

  const mesh::Mesh & _grid;
  std::vector<Point> _qpoints;  // gauss point coordinates in master element
  std::vector<std::array<double,N_ELEMENT_VERTICES<vtk_id>>> _shape_values;
  std::vector<std::array<Point,N_ELEMENT_VERTICES<vtk_id>>> _shape_grads;
  std::vector<double> _weights;
  std::vector<double> _determinants;
  std::array<Point,N_ELEMENT_VERTICES<vtk_id>> _vertex_coord;
};

template<VTK_ID vtk_id>
FeValues<vtk_id>::FeValues(const mesh::Mesh & grid)
    : _grid(grid)
{}


template<VTK_ID vtk_id>
void FeValues<vtk_id>::update(const mesh::Cell & cell)
{
  static_assert(ELEMENT_DIM<vtk_id> == 3, "This function only exists for 3D elements");
  update_vertex_coord_(cell.vertices());
  _qpoints = get_master_integration_points();
  _weights = get_master_integration_weights();
  update_();
}

template<VTK_ID vtk_id>
void FeValues<vtk_id>::update(const mesh::Cell & cell, const angem::Point<3,double> & point)
{
  static_assert(ELEMENT_DIM<vtk_id> == 3, "This function only exists for 3D elements");
  update_vertex_coord_(cell.vertices());
  _weights = {1.0};
  _qpoints = {map_real_to_local_(point)};
  update_();
}

template<VTK_ID vtk_id>
void FeValues<vtk_id>::update(const mesh::Face & face, const angem::Point<3,double> & point)
{
  static_assert(ELEMENT_DIM<vtk_id> == 2, "This function only exists for 2D elements");

  update_vertex_coord_(face.vertices());
  _weights = {1.0};
  _qpoints = {map_real_to_local_(point)};
  update_face_();
}

template<VTK_ID vtk_id>
void FeValues<vtk_id>::update(const mesh::Face & face)
{
  static_assert(ELEMENT_DIM<vtk_id> == 2, "This function only exists for 2D elements");

  update_vertex_coord_(face.vertices());
  _qpoints = get_master_integration_points();
  _weights = get_master_integration_weights();
  update_face_();
}

template<VTK_ID vtk_id>
angem::Tensor2<3, double> FeValues<vtk_id>::compute_jacobian_(const std::array<Point,N_ELEMENT_VERTICES<vtk_id>> & ref_grad) const
{
  // compute jacobian
  // see Becker E., Carey G., Oden J. Finite elements. An Introduction Volume 1 1981
  // Eq. 5.3.6
  // NOTE: this is a transposed jacobian matrix
  angem::Tensor2<3, double> J;
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j)
      for (std::size_t v=0; v<N_ELEMENT_VERTICES<vtk_id>; ++v)
        J( j, i ) += ref_grad[v][j] * _vertex_coord[v][i];
  return J;
}

template<VTK_ID vtk_id>
void FeValues<vtk_id>::update_()
{
  // compute shape values in q-points in master element
  update_shape_values_();

  _determinants.assign(_qpoints.size(), 0.0);
  _shape_grads.assign(_qpoints.size(), std::array<Point,N_ELEMENT_VERTICES<vtk_id>>());
  // loop integration points
  for (std::size_t q=0; q<_qpoints.size(); ++q)
  {
    // compute shape grad in in master element
    const auto ref_grad = compute_ref_gradient_(q);
    // compute transformation jacobian
    angem::Tensor2<3, double> dx_du = compute_jacobian_(ref_grad);
    // compute the determinant of transformation jacobian
    _determinants[q] = det(dx_du);
    // must be positive
    if ( _determinants[q] <= 0 )
    {
      std::cout << "error: " << std::endl;
      std::cout << "vtk_id = " << vtk_id << std::endl;
      std::cout << dx_du << std::endl;
      throw std::runtime_error("Transformation Jacobian is not invertible");
    }
    // invert the jacobian to compute shape function gradients
    angem::Tensor2<3, double> du_dx = invert(dx_du);

    // compute the true shape function gradients
    update_shape_grads_(ref_grad, du_dx, q);
  }
}

template<VTK_ID vtk_id>
void FeValues<vtk_id>::update_face_()
{
  // compute shape values in q-points in master element
  update_shape_values_();

  _determinants.assign(_qpoints.size(), 0.0);
  _shape_grads.assign(_qpoints.size(), std::array<Point,N_ELEMENT_VERTICES<vtk_id>>());

  // loop integration points
  for (std::size_t q=0; q<_qpoints.size(); ++q)
  {
    // compute shape grad in in master element
    const auto ref_grad = compute_ref_gradient_(q);
    const auto [du_dx, detJ] = compute_detJ_and_invert_face_jacobian_(ref_grad);
    _determinants[q] = detJ;
    // compute the true shape function gradients
    update_shape_grads_(ref_grad, du_dx, q);
  }
}

template<VTK_ID vtk_id>
void FeValues<vtk_id>::update_shape_grads_(const std::array<Point,N_ELEMENT_VERTICES<vtk_id>> & ref_grad,
                                           const angem::Tensor2<3, double> & du_dx, const size_t q)
{
    // compute the true shape function gradients
    for (size_t vertex = 0; vertex < N_ELEMENT_VERTICES<vtk_id>; ++vertex)
      for (size_t i=0; i<3; ++i)
        for (size_t j = 0; j < 3; ++j)
          // d phi_vert / dx_i = (d phi_vert / d u_j) * (d u_j / d x_i)
          _shape_grads[q][vertex][i] += ref_grad[vertex][j] * du_dx(j, i);
}

template<VTK_ID vtk_id>
void FeValues<vtk_id>::update_vertex_coord_(const std::vector<size_t> & vertex_indices)
{
  for (size_t iv=0; iv<vertex_indices.size(); ++iv)
    _vertex_coord[iv] = _grid.vertex(vertex_indices[iv]);
}

template<VTK_ID vtk_id>
double FeValues<vtk_id>::value(const size_t shape_index, const size_t qpoint) const
{
  assert( qpoint < _qpoints.size() && "qpoint too large" );
  assert( shape_index < N_ELEMENT_VERTICES<vtk_id> && "shape_index too large");
  return _shape_values[qpoint][shape_index];
}

template<VTK_ID vtk_id>
Point FeValues<vtk_id>::grad(const size_t shape_index, const size_t qpoint) const
{
  assert( qpoint < _qpoints.size() && "qpoint too large" );
  assert( shape_index < N_ELEMENT_VERTICES<vtk_id> && "shape_index too large");
  return _shape_grads[qpoint][shape_index];
}

template<VTK_ID vtk_id>
double FeValues<vtk_id>::JxW(const size_t qpoint) const
{
  assert( qpoint < _qpoints.size() && "qpoint too large" );
  return _determinants[qpoint] * _weights[qpoint];
}


template<VTK_ID vtk_id>
void FeValues<vtk_id>::update_shape_values_()
{
  _shape_values.resize( _qpoints.size() );
  for (std::size_t q=0; q<_qpoints.size(); ++q)
    for (std::size_t v=0; v<N_ELEMENT_VERTICES<vtk_id>; ++v)
      _shape_values[q][v] = eval_(_qpoints[q], v);
}

template<VTK_ID vtk_id>
std::array<Point,N_ELEMENT_VERTICES<vtk_id>> FeValues<vtk_id>::compute_ref_gradient_(const size_t q) const
{
  std::array<Point,N_ELEMENT_VERTICES<vtk_id>> ref_grad;
  for (std::size_t v=0; v<N_ELEMENT_VERTICES<vtk_id>; ++v)
    ref_grad[v] = eval_derivative_(_qpoints[q], v);
  return ref_grad;
}

template<VTK_ID vtk_id>
std::pair<angem::Tensor2<3, double>,double>
FeValues<vtk_id>::
compute_detJ_and_invert_face_jacobian_(const std::array<Point,N_ELEMENT_VERTICES<vtk_id>> & ref_grad) const
{
  // angem::Tensor2<3,double>
  angem::Plane<double> plane (_vertex_coord[0], _vertex_coord[1], _vertex_coord[2]);
  // const auto & basis = plane.get_basis();
  std::vector<Point> loc_coord(_vertex_coord.size());
  for (size_t i=0; i<_vertex_coord.size(); ++i)
    loc_coord[i] = plane.local_coordinates(_vertex_coord[i]);

  angem::Tensor2<3, double> J;
  for (std::size_t i=0; i<3; ++i)
    for (std::size_t j=0; j<3; ++j)
      for (std::size_t v=0; v<N_ELEMENT_VERTICES<vtk_id>; ++v)
        J( j, i ) += ref_grad[v][j] * loc_coord[v][i];

  // since the shape is 2d, we need to cast this down in 2D
  // since the third row and column of J are zero
  const angem::Tensor2<2, double> J2 = {J(0, 0), J(0, 1),
                                        J(1, 0), J(1, 1)};
  const angem::Tensor2<2, double> J2_inv = invert(J2);

  // cast it back to 3D since that's what the code uses to compute shape gradients
  const angem::Tensor2<3, double> J_inv = {J2_inv(0,0), J2_inv(0,1), 0.0,
                                           J2_inv(1,0), J2_inv(1,1), 0.0,
                                           0.0,         0.0,         0.0};
  const double detJ = det(J2);

  return {J_inv, detJ};
}

}  // end namespace discretization

#include "FeValuesTriangle.hpp"
#include "FeValuesTetra.hpp"