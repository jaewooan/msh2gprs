#ifdef WITH_EIGEN
#include "DFEMElement.hpp"
#include "EdgeComparison.hpp"
#include "gmsh_interface/GmshInterface.hpp"
#include "gmsh_interface/FeValues.hpp"
#include "VTKWriter.hpp"
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>

namespace discretization {

using api = gprs_data::GmshInterface;
using FeValues = gprs_data::FeValues;
using Point = angem::Point<3,double>;

DFEMElement::DFEMElement(const mesh::Cell & cell,
                         const double msrsb_tol)
    : _cell(cell), _msrsb_tol(msrsb_tol)
{
  build_();
}

void DFEMElement::build_()
{
  api::build_triangulation(_cell, _element_grid);
  std::cout << "done building discr" << std::endl;
  std::cout << "_element_grid.n_active_cells " << _element_grid.n_active_cells() << std::endl;

  // build system matrix
  build_jacobian_();
  std::cout << "done building jac" << std::endl;

  // Direct method
  build_support_edges_();
  compute_shape_functions_();
  // save_support_edges_();

  //
  // // MSRSB method
  // // // gmsh::write("cell.msh");
  // build_support_boundaries_();
  // initialize_shape_functions_();
  // initial_guess_();
  // // initial_guess2_();
  // _system_matrix.makeCompressed();
  // std::cout << "compression done" << std::endl;
  // debug_save_shape_functions_("shape_functions-initial-guess.vtk");
  // //
  // run_msrsb_();


  // postprocessing
  debug_save_shape_functions_("shape_functions-final.vtk");
  // save_support_boundaries_();
  // find_integration_points_();
  // compute_fe_quantities_();
}

void DFEMElement::build_support_boundaries_()
{
  const std::vector<const mesh::Face*> parent_faces = _cell.faces();
  const std::vector<size_t> parent_vertices = _cell.sorted_vertices();
  _support_boundaries.resize(parent_vertices.size());
  for (auto face = _element_grid.begin_active_faces(); face != _element_grid.end_active_faces(); ++face)
  {
    if (face->marker() > 0 && face->neighbors().size() == 1)
    {
      const mesh::Face * parent_face = parent_faces[face->marker() - 1];
      std::vector<size_t> non_adj;  // parent vertices not in parent face
      for (size_t iv_parent=0; iv_parent<parent_vertices.size(); ++iv_parent)
      {
        // index in global domain grid
        const size_t parent_vertex = parent_vertices[iv_parent];
        if ( !parent_face->has_vertex(parent_vertex) )
          non_adj.push_back(iv_parent);
      }
      for (const size_t vertex : face->vertices())
        for (const size_t parent_vertex : non_adj)
          _support_boundaries[parent_vertex].insert(vertex);
    }
  }
}

void DFEMElement::save_support_edges_()
{
  const std::string fname = "support_edges.vtk";
  std::cout << "saving " << fname << std::endl;

  std::ofstream out;
  out.open(fname.c_str());
  IO::VTKWriter::write_geometry(_element_grid, out);
  const size_t nv = _element_grid.n_vertices();
  IO::VTKWriter::enter_section_point_data(nv, out);

  const size_t n_parent_vertices = _cell.vertices().size();
  for (std::size_t j=0; j<n_parent_vertices; ++j)
  {
    std::vector<double> output(nv, -1);
    for (size_t iv=0; iv<_support_edge_vertices[j].size(); ++iv)
    {
      const size_t v = _support_edge_vertices[j][iv];
      const double value = _support_edge_values[j][iv];
      output[v] = value;
    }
    IO::VTKWriter::add_data(output, "support-"+std::to_string(j), out);
  }
  // save bnd markers
  std::vector<double> output(nv, 0);
  for (auto face = _element_grid.begin_active_faces(); face != _element_grid.end_active_faces(); ++face)
  {
    if (face->marker() > 0 && face->neighbors().size() == 1)
    {
      for ( const size_t v : face->vertices() )
        output[v] = face->marker();
    }
  }
  IO::VTKWriter::add_data(output, "bnd-marker", out);
  out.close();
}

void DFEMElement::save_support_boundaries_()
{
  const std::string fname = "support_bnd.vtk";
  std::cout << "saving " << fname << std::endl;

  std::ofstream out;
  out.open(fname.c_str());
  IO::VTKWriter::write_geometry(_element_grid, out);
  const size_t nv = _element_grid.n_vertices();
  IO::VTKWriter::enter_section_point_data(nv, out);

  const size_t n_parent_vertices = _cell.vertices().size();
  for (std::size_t j=0; j<n_parent_vertices; ++j)
  {
    std::vector<double> output(nv, 0);
    for (const size_t v : _support_boundaries[j])
      output[v] = 1;
    IO::VTKWriter::add_data(output, "support-"+std::to_string(j), out);
  }
  // save bnd markers
  std::vector<double> output(nv, 0);
  for (auto face = _element_grid.begin_active_faces(); face != _element_grid.end_active_faces(); ++face)
  {
    if (face->marker() > 0 && face->neighbors().size() == 1)
    {
      for ( const size_t v : face->vertices() )
        output[v] = face->marker();
    }
  }
  IO::VTKWriter::add_data(output, "bnd-marker", out);
  out.close();
}

void DFEMElement::run_msrsb_()
{
  JacobiPreconditioner preconditioner(_system_matrix);
  std::vector<Eigen::VectorXd> solutions;
  for (size_t i=0; i < _cell.vertices().size(); ++i)
    solutions.push_back(Eigen::VectorXd::Zero( _element_grid.n_vertices() ));

  double dphi = 10  * _msrsb_tol;
  size_t iter = 0;
  const size_t max_iter = 10 * _element_grid.n_vertices();
  do
  {
    dphi = jacobi_iteration_(solutions, preconditioner) / solutions[0].size();

    if (!(iter % 50))  // li'l progress prinout
      std::cout << "iter = " << iter << " dphi = " << dphi << std::endl;

    if (iter++ >= max_iter)
    {
      std::cout << "cell index = " << _cell.index() << std::endl;
      api::save_gmsh_grid("unconverged.msh");
      debug_save_shape_functions_("shape_functions-unconverged.vtk");
      throw std::runtime_error("msrsb did not converge after " + std::to_string(iter) +
                               " iterations (current error = " + std::to_string(dphi) + ")");
    }

  } while (dphi > _msrsb_tol);
  std::cout << "MSRSB converged after " << iter << " iterations" << std::endl;
}

double DFEMElement::jacobi_iteration_(std::vector<Eigen::VectorXd> & solutions,
                                      const JacobiPreconditioner & prec)
{
  for (size_t parent_node = 0; parent_node < _cell.vertices().size(); parent_node++)
    prec.solve(_system_matrix, _basis_functions[parent_node], solutions[parent_node]);

  for (size_t vertex=0; vertex < _element_grid.n_vertices(); vertex++)
  {
    enforce_zero_on_boundary_(vertex, solutions);
    enforce_partition_of_unity_(vertex, solutions);
  }

  double error = 0;
  for (size_t parent_node = 0; parent_node < _cell.vertices().size(); parent_node++)
  {
    _basis_functions[parent_node] += solutions[parent_node];
    const double value = solutions[parent_node].norm();
    error = std::max(value, error);
  }
  return error;
}

void DFEMElement::enforce_partition_of_unity_(const size_t fine_vertex,
                                              std::vector<Eigen::VectorXd> & solutions)
{
  const size_t parent_n_vert = _cell.vertices().size();
  double sum_shape_functions = 0;
  for (size_t parent_vertex = 0; parent_vertex < parent_n_vert; parent_vertex++)
  {
    const double shape = _basis_functions[parent_vertex][fine_vertex];
    const double new_value = shape + solutions[parent_vertex][fine_vertex];
    sum_shape_functions += new_value;
  }
  for (size_t parent_vertex = 0; parent_vertex < parent_n_vert; parent_vertex++)
  {
    const double shape = _basis_functions[parent_vertex][fine_vertex];
    const double new_value = (shape + solutions[parent_vertex][fine_vertex]) / sum_shape_functions;
    solutions[parent_vertex][fine_vertex] = new_value - shape;
  }
}

void DFEMElement::enforce_zero_on_boundary_(const size_t fine_vertex,
                                            std::vector<Eigen::VectorXd> & solutions)
{
  const size_t parent_n_vert = _cell.vertices().size();
  for (size_t parent_vertex = 0; parent_vertex < parent_n_vert; parent_vertex++)
  {
    if (in_support_boundary_(fine_vertex, parent_vertex))
      solutions[parent_vertex][fine_vertex] = 0;
  }
}

void DFEMElement::build_jacobian_()
{
  // build element jacobian for the homogeneous laplace equation
  _system_matrix = Eigen::SparseMatrix<double,Eigen::RowMajor>(_element_grid.n_vertices(),
                                                               _element_grid.n_vertices());
  // since we only build tetrahedral element mesh
  const int element_type = api::get_gmsh_element_id(angem::VTK_ID::TetrahedronID);
  FeValues fe_values( element_type );

  Eigen::MatrixXd cell_matrix(fe_values.n_vertices(), fe_values.n_vertices());
  for (auto cell = _element_grid.begin_active_cells(); cell != _element_grid.end_active_cells(); ++cell)
  {
    cell_matrix.setZero();
    fe_values.update(cell->index());

    const std::vector<size_t> & cell_vertices = cell->vertices();
    const size_t nv = cell_vertices.size();

    // This assembles a local matrix for the laplace equation
    for (size_t q = 0; q < fe_values.n_q_points(); ++q)
    {
      for (size_t i = 0; i < nv; ++i)
        for (size_t j = 0; j < nv; ++j)
          cell_matrix(i, j) += -(fe_values.grad(i, q) * // grad phi_i(x_q)
                                 fe_values.grad(j, q) * // grad phi_j(x_q)
                                 fe_values.JxW(q));     // dV
    }

    /* distribute local to global */
    for (size_t i = 0; i < nv; ++i)
      for (size_t j = 0; j < nv; ++j)
      {
        const size_t idof = cell_vertices[i];
        const size_t jdof = cell_vertices[j];
        _system_matrix.coeffRef(idof, jdof) += cell_matrix(i, j);
      }
  }
}

void DFEMElement::initialize_shape_functions_()
{
  const size_t parent_n_vert = _cell.vertices().size();
  for (std::size_t i=0; i<parent_n_vert; ++i)
    _basis_functions.push_back(Eigen::VectorXd::Zero(_element_grid.n_vertices()));
}

void DFEMElement::initial_guess_()
{
  const auto parent_nodes = _cell.polyhedron()->get_points();
  for (size_t v=0; v < _element_grid.n_vertices(); ++v)
  {
    double min_dist = std::numeric_limits<double>::max();
    size_t closest = 0;
    for (size_t j = 0; j < parent_nodes.size(); ++j)
    {
      const auto & pn = parent_nodes[j];
      const double dist = pn.distance(_element_grid.vertex(v));
      if (dist < min_dist)
      {
        min_dist = dist;
        closest = j;
      }
    }
    _basis_functions[closest][v] = 1.0;
  }
}

void DFEMElement::initial_guess2_()
{
  const auto parent_nodes = _cell.polyhedron()->get_points();
  for (size_t v=0; v < _element_grid.n_vertices(); ++v)
  {
    const Point & x = _element_grid.vertex(v);
    for (size_t i = 0; i < parent_nodes.size(); ++i)
    {
      const Point & xi = parent_nodes[i];
      double num = 1, denom = 1;
      for (size_t j = 0; j < parent_nodes.size(); ++j)
        if (i != j)
        {
          const Point & xj = parent_nodes[j];
          num *= x.distance(xj);
          denom *= xi.distance(xj);
        }
      _basis_functions[i][v] = num / denom;
    }
  }
  for (size_t v=0; v < _element_grid.n_vertices(); ++v)
    for (size_t i = 0; i < parent_nodes.size(); ++i)
    {
      const double val = _basis_functions[i][v];
      _basis_functions[i][v] = std::max(0.0, std::min( val, 1.0 ) );
    }
}

void DFEMElement::debug_save_shape_functions_(const std::string fname)
{
  std::cout << "saving " << fname << std::endl;
  std::ofstream out;
  out.open(fname.c_str());

  IO::VTKWriter::write_geometry(_element_grid, out);
  const size_t nv = _element_grid.n_vertices();
  IO::VTKWriter::enter_section_point_data(nv, out);

  const size_t n_parent_vertices = _cell.vertices().size();
  for (std::size_t j=0; j<n_parent_vertices; ++j)
  {
    std::vector<double> output(nv, 0.0);
    for (size_t i = 0; i < nv; ++i)
      output[i] = _basis_functions[j][i];
    IO::VTKWriter::add_data(output, "basis-"+std::to_string(j), out);
  }
  out.close();
}

bool DFEMElement::in_support_boundary_(const size_t fine_vertex, const size_t parent_vertex) const
{
  if ( _support_boundaries[parent_vertex].find(fine_vertex) == _support_boundaries[parent_vertex].end() )
    return false;
  else return true;
}

bool DFEMElement::in_global_support_boundary_(const size_t fine_vertex) const
{
  for (size_t parent_vertex = 0; parent_vertex < _cell.vertices().size(); parent_vertex++)
    if (in_support_boundary_(fine_vertex, parent_vertex))
      return true;
  return false;
}


void DFEMElement::find_integration_points_()
{
  /* Split a parent cell into pyramids and find their centers */
  const auto polyhedron = _cell.polyhedron();
  std::vector<Point> parent_vertices = polyhedron->get_points();
  parent_vertices.push_back(_cell.center());
  for (const auto & face : polyhedron->get_faces())
    _cell_gauss_points.push_back(create_pyramid_and_compute_center_(face, parent_vertices));
}

Point DFEMElement::create_pyramid_and_compute_center_(const std::vector<size_t> & face,
                                                      const std::vector<Point> & vertices) const
{
  const size_t vertex_center = vertices.size() - 1;
  std::vector<std::vector<size_t>> pyramid_faces;
  for (size_t iv=0; iv<face.size(); ++iv)
  {
    size_t v1, v2;
    v1 = face[iv];
    if ( iv < face.size() - 1 )
      v2 = face[iv+1];
    else
      v2 = face[0];
    pyramid_faces.push_back( {v1, v2, vertex_center} );
  }

  angem::Polyhedron<double> pyramid(vertices, pyramid_faces);
  return pyramid.center();
}

void DFEMElement::compute_fe_quantities_()
{
  // allocate vector of point data
  const size_t n_parents = _basis_functions.size();
  _cell_data.points.resize(_cell_gauss_points.size());
  // 4 is a gmsh tetra
  FeValues fe_values( 4, _element_grid.n_cells() );
  for (auto cell = _element_grid.begin_active_cells(); cell != _element_grid.end_active_cells(); ++cell)
    for (size_t q=0; q<_cell_gauss_points.size(); ++q)
      if ( cell->polyhedron()->point_inside(_cell_gauss_points[q]) )
      {
        std::vector<angem::Point<3,double>> points = {_cell_gauss_points[q]};
        const std::vector<size_t> & cell_verts = cell->vertices();
        fe_values.update(cell->index(), points);
        const size_t nv = cell->vertices().size();
        const size_t q_loc = 0;  // only one integration point
        FEPointData data;
        data.values.resize( _basis_functions.size() );
        data.grads.resize( _basis_functions.size() );
        for (size_t parent_vertex=0; parent_vertex<n_parents; ++parent_vertex)
        {
          data.values[parent_vertex] = 0;
          for (size_t vertex=0; vertex<nv; ++vertex)
          {
            data.values[parent_vertex] += fe_values.value( vertex, q_loc ) *
                                          _basis_functions[parent_vertex][cell_verts[vertex]];
            data.grads[parent_vertex] += fe_values.grad(vertex, q_loc) *
                                         _basis_functions[parent_vertex][cell_verts[vertex]];
          }
        }
        data.weight = 1.0 / static_cast<double>( _cell_gauss_points.size() );
        _cell_data.points.push_back( std::move(data) );
      }
}

const FiniteElementData & DFEMElement::get_cell_data() const
{
  return _cell_data;
}

void DFEMElement::build_support_edges_()
{
  /* Build structure that store face markers for each vertex in grid */
  std::vector<std::vector<size_t>> vertex_markers( _element_grid.n_vertices() );
  for (auto face = _element_grid.begin_active_faces(); face != _element_grid.end_active_faces(); ++face)
    if (face->neighbors().size() == 1 && face->marker() > 0)
    {
      const size_t ipf = face->marker();  // i parent face
      for (size_t v : face->vertices())
        if (std::find( vertex_markers[v].begin(), vertex_markers[v].end(), ipf ) ==
            vertex_markers[v].end())
          vertex_markers[v].push_back(ipf);
    }

  // map parent vertex to parent faces
  std::vector<std::list<size_t>> parent_vertex_markers( _cell.n_vertices() );
  const std::vector<size_t> parent_vertices = _cell.sorted_vertices();
  const std::vector<const mesh::Face*> parent_faces = _cell.faces();
  for (size_t ipf=0; ipf<parent_faces.size(); ++ipf)
  {
    const auto parent_face = parent_faces[ipf];

    for (size_t iv_parent=0; iv_parent<parent_vertices.size(); ++iv_parent)
    {
      const size_t parent_vertex = parent_vertices[iv_parent];
      if (parent_face->has_vertex(parent_vertex))
        parent_vertex_markers[iv_parent].push_back( ipf + 1);
    }
  }

  const auto pair_markers_to_edge = edgecmp::EdgeComparison::get_edges( parent_vertex_markers );

  _support_edge_vertices.resize( parent_vertices.size() );
  _support_edge_values.resize( parent_vertices.size() );
  const auto parent_nodes = _cell.polyhedron()->get_points();
  for (size_t v=0; v<_element_grid.n_vertices(); ++v)
  {
    if ( vertex_markers[v].size() > 1 )
    {
      const auto & markers = vertex_markers[v];
      for (size_t mi=0; mi<markers.size(); ++mi)
        for (size_t mj=mi+1; mj<markers.size(); ++mj)
        {
          const size_t m1 = markers[mi];
          const size_t m2 = markers[mj];
          assert( pair_markers_to_edge.contains(m1, m2) );
          const auto & parent_edge = pair_markers_to_edge.get_data(m1, m2);
          const size_t vp1 = parent_edge.either();
          const size_t vp2 = parent_edge.other(vp1);
          const Point p1 = parent_nodes[vp1];
          const Point p2 = parent_nodes[vp2];
          const double dp = p1.distance(p2);
          const double d1 = _element_grid.vertex(v).distance(p1);
          const double d2 = _element_grid.vertex(v).distance(p2);
          _support_edge_vertices[vp1].push_back(v);
          _support_edge_values[vp1].push_back( (dp - d1 ) / dp );
          _support_edge_vertices[vp2].push_back(v);
          _support_edge_values[vp2].push_back( (dp - d2) / dp );

          // also set zero on edges that are not either parent
          for (size_t vp=0; vp<parent_vertices.size(); ++vp)
            if (vp != vp1 && vp != vp2)
            {
              _support_edge_vertices[vp].push_back(v);
              _support_edge_values[vp].push_back(0.0);
            }

        }
    }
    else if (vertex_markers[v].size() == 1)
    {
      const size_t iface = vertex_markers[v].front();
      for (size_t pv = 0; pv < _cell.n_vertices(); pv++)
      {
        if ( std::find( parent_vertex_markers[pv].begin(), parent_vertex_markers[pv].end(),
                        iface) == parent_vertex_markers[pv].end() )
        {
          _support_edge_vertices[pv].push_back(v);
          _support_edge_values[pv].push_back(0.0);
        }
      }
    }
  }
}

void DFEMElement::compute_shape_functions_()
{
  const size_t n_parent_vertices = _cell.vertices().size();

  // allocate shape functions
  for (std::size_t i=0; i<n_parent_vertices; ++i)
    _basis_functions.push_back(Eigen::VectorXd::Zero(_element_grid.n_vertices()));

  for (size_t pv=0; pv<n_parent_vertices; ++pv)
  {
    // assemble problem with BC's for the parent vertex
    // copy system matrix and create a rhs vector
    Eigen::SparseMatrix<double,Eigen::RowMajor> mat = _system_matrix;
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero( _element_grid.n_vertices() );
    // impose boundary conditions
    impose_boundary_conditions_(mat, rhs, pv);
    mat.makeCompressed();

    // TODO: I cannot currently use conjugate gradient since the matrix is non-symmetric
    // due to dirichlet conditions. I need to symmetrize the matrix by modifying the
    // columns that correspond to dirichlet vertices
    // Eigen::ConjugateGradient<Eigen::SparseMatrix<double,Eigen::RowMajor>, Eigen::Lower|Eigen::Upper> solver;
    // solver.analyzePattern(mat);
    // solver.factorize(mat);
    // solver.setMaxIterations(200);
    // solver.setTolerance(1e-10);
    // solver.compute(mat);

    // Eigen::SparseLU<Eigen::SparseMatrix<double,Eigen::RowMajor>> solver;
    // Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver(mat);
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver(mat);
    // solver.factorize(mat);
    solver.analyzePattern(mat);
    solver.factorize(mat);

    _basis_functions[pv] = solver.solve(rhs);
    // if ( solver.info() ==  Eigen::ComputationInfo::NumericalIssue)
    //   std::cout << "numerical issue" << std::endl;
    // else if( solver.info() ==  Eigen::ComputationInfo::NoConvergence)
    //   std::cout << "no convergence" << std::endl;
    if (solver.info() != Eigen::ComputationInfo::Success)
    {
      throw std::runtime_error( "Conjugate gratient did not converge" );
    }
    std::cout << "done solving " << pv << std::endl;
  }

  for (size_t i=0; i<_element_grid.n_vertices(); ++i)
  {
    double sum = 0;
    for (size_t j=0; j<n_parent_vertices; ++j)
      sum += _basis_functions[j][i];
    std::cout << "sum = " << sum << std::endl;
  }

}

void DFEMElement::impose_boundary_conditions_(Eigen::SparseMatrix<double,Eigen::RowMajor> & mat,
                                              Eigen::VectorXd & rhs,
                                              const size_t ipv)
{
  for (size_t iv=0; iv<_support_edge_vertices[ipv].size(); ++iv)
  {
    const size_t i = _support_edge_vertices[ipv][iv];
    for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(mat, i); it; ++it)
      it.valueRef() = (it.row() == it.col()) ? 1.0 : 0.0;
    rhs[i] = _support_edge_values[ipv][iv];
  }
}

}  // end namespace discretization

#endif  // with_eigen
