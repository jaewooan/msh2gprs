#pragma once

#include "Mesh.hpp"

namespace mesh {

/*
** This class implements splitting mesh faces for Gomechanics DFM fractures.
** It does not actually modify the mesh since AD-GPRS needs only the information
** about the master faces.
 */
class FaceSplitter {
 public:
  FaceSplitter(const Mesh & grid);
  // tell splitter which faces to split before calling split_faces method
  void mark_for_split(const size_t face_index);
  // split faces marked for splitting with mark_for_split
  // returns SurfaceMesh of master DFM faces
  // cleans marked_for_split array upon completion
  SurfaceMesh<double> split_faces();

  const std::vector<angem::Point<3,double>> & get_all_vertices() const {return _vertex_coord;}
  const std::vector<std::vector<size_t>> & get_cell_vertices() const {return _cell_vertices;}

  const std::unordered_map<size_t, std::vector<size_t>> & get_child_vertices() const
  { return _parent_to_child_vertices;}

 private:
  void split_vertex_(const std::size_t               vertex_index,
                     const std::vector<std::size_t> &splitted_face_indices);
  // two elements are in the same group if they are neighbors and
  // the neighboring face is not in vertex_faces array
  std::vector<std::vector<std::size_t>>
  group_cells_based_on_split_faces_(const std::vector<std::size_t> & affected_cells,
                                    const std::vector<std::size_t> & splitted_face_indices) const;
  std::vector<std::vector<std::size_t>>
  group_cells_based_on_split_faces_2(const std::vector<std::size_t> & affected_cells,
                                     const std::vector<std::size_t> & splitted_face_indices) const;

  void create_fracture_face_grid_(SurfaceMesh<double> & face_mesh,
                                  std::unordered_map<size_t,size_t> & map_face_surface_element,
                                  std::unordered_map<size_t,size_t> & map_frac_vertex_vertex) const;
  void find_vertices_to_split_(const SurfaceMesh<double> & face_mesh);
  void find_vertices_to_split_2(const SurfaceMesh<double> & face_mesh);

  std::unordered_set<size_t> find_fracture_tip_vertices_(const SurfaceMesh<double> & face_mesh);
  std::unordered_set<size_t> find_global_boundary_();

  // ----------------------- Variables ----------------------- //
  const Mesh & _grid;
  // vector of grid vertex coordinates + those that happened after split
  std::vector<angem::Point<3,double>> _vertex_coord;
  // vector of grid cell vertex coordinates after splitting
  std::vector<std::vector<size_t>> _cell_vertices;
  std::vector<size_t> _marked_for_split;
  // map surfacemesh vertex -> 3d mesh vertex
  std::unordered_map<size_t,size_t> _surface_vertex_to_global;
  // map 2d-element -> 3d face index
  std::unordered_map<size_t,size_t> _surface_to_face;
  std::unordered_map<size_t, std::vector<size_t>> _vertices_to_split;
  std::unordered_map<size_t, std::vector<size_t>> _parent_to_child_vertices;
};

}  // end namespace mesh
