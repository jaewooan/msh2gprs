#include "DFEMElement.hpp"
#include "gmsh_interface/GmshInterface.hpp"
#include "gmsh_interface/FeValues.hpp"

namespace discretization {

using api = gprs_data::GmshInterface;
using FeValues = gprs_data::FeValues;

DFEMElement::DFEMElement(const mesh::Cell & cell)
    : _cell(cell)
{
  build_();
}

void DFEMElement::build_()
{
  api::build_triangulation(_cell);
  build_jacobian_();
  initial_guess_();
}

void DFEMElement::build_jacobian_()
{
  // 1. ask gmsh to provide gaussian points, shape functions, and jacobians
  // for our tetras
  // 2. build element jacobian for the homogeneous laplace equation
  api::get_elements(_element_types, _element_tags, _node_tags, /* dim = */ 3);
  numberNodesEndElements_(_element_types, _element_tags, _node_tags);
  _system_matrix = Eigen::SparseMatrix<double,Eigen::RowMajor>(_node_numbering.size(),
                                                              _node_numbering.size());

  for (std::size_t itype = 0; itype < _element_types.size(); ++itype)
  {
    const size_t type = _element_types[itype];
    FeValues fe_values(type, _element_tags[itype].size());

    for (std::size_t itag=0; itag<_element_tags[itype].size(); ++itag)
    {
      const size_t tag = _element_tags[itype][itag];
      const size_t cell_number = _cell_numbering[tag];
      fe_values.update(cell_number);
      Eigen::MatrixXd cell_matrix = Eigen::MatrixXd::Zero(fe_values.n_vertices(),
                                                          fe_values.n_vertices());

      for (size_t q = 0; q < fe_values.n_q_points(); ++q)
      {
        for (size_t i = 0; i < fe_values.n_vertices(); ++i)
          for (size_t j = 0; j < fe_values.n_vertices(); ++j)
            cell_matrix(i, j) += (fe_values.grad(i, q) * // grad phi_i(x_q)
                                  fe_values.grad(j, q) * // grad phi_j(x_q)
                                  fe_values.JxW(q));           // dx
      }

      for (size_t i = 0; i < fe_values.n_vertices(); ++i)
        for (size_t j = 0; j < fe_values.n_vertices(); ++j)
        {
          const size_t inode = _node_tags[itype][fe_values.n_vertices()*itag + i];
          const size_t jnode = _node_tags[itype][fe_values.n_vertices()*itag + j];
          const size_t idof = _node_numbering[inode];
          const size_t jdof = _node_numbering[jnode];
          _system_matrix.coeffRef(idof, jdof) += cell_matrix(i, j);
        }
    }
  }
  _system_matrix.makeCompressed();
}

void DFEMElement::
numberNodesEndElements_(std::vector<int> &element_types,
                        std::vector<std::vector<std::size_t> > & element_tags,
                        const std::vector<std::vector<std::size_t> > &node_tags)
{
  _node_numbering.clear();
  _cell_numbering.clear();
  size_t dof_cell = 0, dof_node = 0;
  for (std::size_t itype = 0; itype < element_types.size(); ++itype)
  {
    const size_t nv = api::get_n_vertices(element_types[itype]);

    for (std::size_t itag=0; itag<element_tags[itype].size(); ++itag)
    {
      _cell_numbering[itag] = dof_cell++;
      for (std::size_t inode=0; inode<nv; ++inode)
      {
        const size_t node = node_tags[itype][nv * itag + inode];
        if (_node_numbering.find(node) == _node_numbering.end())
          _node_numbering[node] = dof_node++;
      }
    }
  }
}

void DFEMElement::initial_guess_()
{
  const size_t parent_n_vert = _cell.vertices().size();
  for (std::size_t i=0; i<parent_n_vert; ++i)
    _basis_functions.emplace_back(_node_numbering.size());

  // std::cout << "vertices" << std::endl;

  // get nodes that form the vertices of the original elements
  std::vector<size_t> node_tags;
  std::vector<double> coord, parametric_coord;
gmsh::model::mesh::getNodes(node_tags, coord, parametric_coord,
                          0, 0, /* incl boundary */ true,
                          /* return parametric = */ false);
  for (auto & it : node_tags)
      std::cout << it << std::endl;
  std::cout << "coord" << std::endl;
  for (auto v : coord)
    std::cout << v << std::endl;
  std::cout << "cell v" << std::endl;
  std::cout << _cell.polyhedron()->get_points()[0] << std::endl;

}


}  // end namespace discretization
