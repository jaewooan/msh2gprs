#pragma once

#include <angem/Point.hpp>
#include <angem/Polyhedron.hpp>
#include <angem/Polygon.hpp>
#include <mesh_methods.hpp>
#include <ShapeID.hpp>
#include <Face.hpp>
#include <cell_iterator.hpp>
#include <face_iterator.hpp>

#include <algorithm> // std::sort


/*
 * This class implements a structure for unstructure grid storage
 * It features constant lookup and insertion times
 * for faces, cells, and their neighbors.
 * Due to this feature though, the order of faces is not guaranteed
 * so that unordered_map[hash]->value should be used to store
 * data related to faces
 */

namespace mesh
{


using Point = angem::Point<3,double>;
// using Polygon = angem::Polygon<double>;
using Polyhedron = angem::Polyhedron<double>;
using Polygon = angem::Polygon<double>;
// using Face = std::vector<std::size_t>;
using FaceiVertices = std::vector<std::size_t>;
// using hash_type = uint256_t;
// using hash_type = vli::integer<128>;
// using FaceMap = std::unordered_map<hash_type, Face>;


// maximum 6-vert polygons as faces (hashing)
// this allows hashing about (2^256)^(1/6) ≈ 7x10¹² vertices
class Mesh
{
 public:
  Mesh();
  // this method does not allow duplicates in vertices
  // infers type by number of vertices
  void insert(const Polyhedron & poly,
              const int          marker = -1);
  // same as above but without searching vertices
  void insert_cell(const std::vector<std::size_t> & ivertices,
                   const int                        vtk_id,
                   const int                        marker = -1);
  // insert marker into map_physical_faces
  void insert(const Polygon & poly,
              const int       marker);
  // insert element without searching vertices
  void insert_face(const std::vector<std::size_t> & ivertices,
                   const int                        vtk_id,
                   const int                        marker);

  // iterators
  // cell iterators
  cell_iterator create_cell_iterator(const std::size_t icell)
  {return cell_iterator(icell, vertices, cells, map_faces,
                        shape_ids, cell_markers);}
  cell_iterator begin_cells(){return create_cell_iterator(0);}
  cell_iterator end_cells()  {return create_cell_iterator(cells.size());}

  // face iterators
  face_iterator create_face_iterator(const FaceMap::iterator & it)
  {return face_iterator(it, vertices);}
  face_iterator begin_faces(){return create_face_iterator(map_faces.begin());}
  face_iterator end_faces()  {return create_face_iterator(map_faces.end());}

  // GETTERS
  // get vector of neighbor cell indices
  std::vector<std::size_t> get_neighbors( const std::size_t icell ) const;
  // vector of indices of cells neighboring a face
  const std::vector<std::size_t> & get_neighbors( const FaceiVertices & face ) const;
  // get vector of vectors of indices representing faces of a cell
  std::vector<std::vector<std::size_t>> get_faces( const std::size_t ielement ) const;
  // true if vector of cells is empty
  bool empty() const {return cells.empty();}
  // get number of cells
  inline std::size_t n_cells() const {return cells.size();}
  // get number of vertices
  inline std::size_t n_vertices() const {return vertices.size();}
  // get number of faces
  inline std::size_t n_faces() const {return map_faces.size();}
  // get cell center coordinates
  Point get_center(const std::size_t icell) const;
  std::unique_ptr<Polyhedron> get_polyhedron(const std::size_t icell) const;

  // MANIPULATION
  // delete an element from the mesh
  void delete_element(const std::size_t ielement);
  // merges jcell into icell if they have a common face
  std::size_t merge_cells(const std::size_t icell,
                          const std::size_t jcell);
  // tell grid which faces to split before calling split_faces method
  void mark_for_split(const face_iterator & face);
  // split faces marked for splitting with mark_for_split
  void split_faces();

  // ATTRIBUTES
  angem::PointSet<3,double>             vertices;
  std::vector<std::vector<std::size_t>> cells;  // vertex indices
  // map face -> neighbor elements
  std::unordered_map<hash_type, Face> map_faces;
  std::vector<int> shape_ids;
  std::vector<int> cell_markers;

 private:
  /* split a vertex
   * retults in adding new vertices (pushed to the back of vertices set)
   * modifies elements of cells array
   */
  void split_vertex(const std::size_t                            ivertex,
                    const std::vector<face_iterator>             & vertex_faces,
                    std::unordered_map<std::size_t, std::size_t> & map_old_new_cells,
                    std::vector<std::vector<std::size_t>>        & new_cells);

  // two elements are in the same group if they are neighbors and
  // the neighboring face is not in vertex_faces array
  std::vector<std::vector<std::size_t>>
  group_cells_based_on_split_faces(const std::unordered_set<std::size_t> & affected_cells,
                                   const std::vector<face_iterator>      & vertex_faces) const;

  // get global indices of polygon face vertices
  std::vector<std::vector<std::size_t>> get_faces(const Polyhedron & poly) const;
  std::vector<hash_type> marked_for_split;

};


}
