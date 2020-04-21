#pragma once

#include "Face.hpp"
#include "Cell.hpp"
#include "active_cell_iterator.hpp"
#include "active_cell_const_iterator.hpp"
#include "active_face_const_iterator.hpp"
#include "SurfaceMesh.hpp"  // mesh::SurfaceMesh
#include <angem/Point.hpp> // angem::Point
#include <angem/Polyhedron.hpp> // angem::Polyhedron
#include <angem/Polygon.hpp> // angem::Polygon
#include <algorithm> // std::sort
#include <unordered_set> // unordered_set


namespace mesh
{

using Point = angem::Point<3,double>;
using Polyhedron = angem::Polyhedron<double>;
using Polygon = angem::Polygon<double>;
using FaceiVertices = std::vector<std::size_t>;

// buddies
class Subdivision;
namespace io { class VTKReader;}

struct FaceTmpData
{
  std::vector<size_t> vertices;
  size_t parent = constants::invalid_index;
  int marker = constants::default_face_marker;
  int vtk_id = angem::VTK_ID::GeneralPolyhedronID;
};


/* This class implements a structure for unstructure grid storage
 * It features constant lookup and insertion times
 * for faces, cells, and their neighbors.
 * Due to this feature though, the order of faces is not guaranteed
 * so that unordered_map[hash]->value should be used to store
 * data related to faces */
class Mesh
{
 public:
  /**
   * Default constructor.
   * Creates an empty grid object with zero cells, faces, and vertices.
   */
  Mesh();
  /**
   * Copy assignment operator
   * 
   * @param  {Mesh} other : a mesh object to be copied
   * @return {Mesh}       : a reference to this.
   */
  Mesh & operator=(const Mesh & other);
  /**
   * Insert a standard vtk polyhedron cell into grid.
   * This will always insert a new cell into the grid, even if a cell with the 
   * same vertices and faces exists.
   * @param  {std::vector<std::size_t>} ivertices : vector of grid cell vertices
   * @param  {int} vtk_id                         : vtk_id of the cell polyhedron
   * @param  {int} marker                         : cell marker (aka domain id)
   * @return {std::size_t}                        : index of the inserted cell
   */
  std::size_t insert_cell(const std::vector<std::size_t> & ivertices,
                          const int                        vtk_id,
                          const int                        marker = constants::default_face_marker);
  /**
   * Insert a face into the grid. If a face with the same vertices exists,
   * this function will rewrite its vtk_id, optionally marker, and optionally parent.
   * A new face will not be created if a face with the same vertices is found.
   * 
   * @param  {std::vector<std::size_t>} ivertices : vector of face vertex indices
   * @param  {int} vtk_id                         : vtk id of the face shape
   * @param  {int} marker                         : face marker
   * @param  {std::size_t} face_parent            : optional, face parent
   * @return {std::size_t}                        : inserted face index
   */
  std::size_t insert_face(const std::vector<std::size_t> & ivertices,
                          const int                        vtk_id,
                          const int                        marker = constants::default_face_marker,
                          const std::size_t                face_parent = constants::invalid_index);

  /**
   * Insert a new vertex into the grid.
   * Returns the index of the new vertex.
   * There is no checks performed whether there are duplicates.
   */
  std::size_t insert_vertex(const angem::Point<3,double> & coord);

  // ITERATORS
  //  create cell iterator for the first active cell
  active_cell_iterator begin_active_cells();
  //  create cell iterator for the first active cell
  active_cell_const_iterator begin_active_cells() const;
  // Returns an iterator referring to the past-the-end active cell
  inline active_cell_iterator end_active_cells() {return active_cell_iterator(nullptr);}
  // Returns an iterator referring to the past-the-end active cell
  inline active_cell_const_iterator end_active_cells() const {return active_cell_const_iterator(nullptr);}

  /* -------------------------- RAW cell iterators, use with caution ---------- */
  //  create cell iterator for the first cell
  //  NOTE: RAW cell iterator, use with caution
  inline std::vector<Cell>::iterator begin_cells() {return m_cells.begin();}
  // create cell const_iterator for the first cell
  //  NOTE: RAW cell iterator, use with caution
  inline std::vector<Cell>::const_iterator begin_cells() const {return m_cells.begin();}
  // Returns an iterator referring to the past-the-end element in the cell container
  //  NOTE: RAW cell iterator, use with caution
  inline std::vector<Cell>::iterator end_cells()   {return m_cells.end();}
  // end iterator for cells. Must only be used as the range indicator
  //  NOTE: RAW cell iterator, use with caution
  inline std::vector<Cell>::const_iterator end_cells() const {return m_cells.end();}

  /* -------------------------- Active cell iterators ---------------------------- */
  // create an active face iterator
  active_face_const_iterator begin_active_faces() const;
  // create an active face end-iterator
  inline active_face_const_iterator end_active_faces() const {return active_face_const_iterator(nullptr, m_faces);}

  /* -------------------------- RAW face iterators, use with caution ---------- */
  // create a face iterator
  //  NOTE: RAW cell iterator, use with caution
  inline std::vector<Face>::iterator begin_faces(){return m_faces.begin();}
  // create a face const_iterator
  //  NOTE: RAW cell iterator, use with caution
  inline std::vector<Face>::const_iterator begin_faces() const {return m_faces.begin();}
  // create an end iterator for faces
  //  NOTE: RAW cell iterator, use with caution
  inline std::vector<Face>::iterator end_faces(){return m_faces.end();}
  // create an end const_iterator for faces
  //  NOTE: RAW cell iterator, use with caution
  inline std::vector<Face>::const_iterator end_faces() const {return m_faces.end();}

  /* -------------------------- Vertex iterators ------------------------------- */
  // create a vertex iterator pointing to the first vertex
  inline std::vector<Point>::iterator begin_vertices() { return m_vertices.begin(); }
  // create a vertex const_iterator pointing to the first vertex
  inline std::vector<Point>::const_iterator begin_vertices() const { return m_vertices.begin(); }
  // create a vertex iterator pointing to the past last vertex
  inline std::vector<Point>::iterator end_vertices() { return m_vertices.end(); }
  // create a vertex const_iterator pointing to the past last vertex
  inline std::vector<Point>::const_iterator end_vertices() const { return m_vertices.end(); }

  /* -------------------------- Access operators (Getters) ---------------------- */
  // get reference to a cell object by index
  inline Cell & cell(const size_t cell_index) { return m_cells[cell_index]; }
  // get const reference to a cell object by index
  inline const Cell & cell(const size_t cell_index) const { return m_cells[cell_index]; }
  // get reference to a face object by index
  inline Face & face(const size_t face_index) { return m_faces[face_index]; }
  // get const reference to a face object by index
  inline const Face & face(const size_t face_index) const { return m_faces[face_index]; }
  // get reference to a vertex object by index
  inline Point& vertex(const size_t vertex_index) {return m_vertices[vertex_index];}
  // get const reference to a vertex object by index
  inline const Point& vertex(const size_t vertex_index) const {return m_vertices[vertex_index];}
  // get vector of all the grid vertex node coordinates
  inline std::vector<angem::Point<3,double>> & vertices() {return m_vertices;}
  // get const vector of all the grid vertex node coordinates
  inline const std::vector<angem::Point<3,double>> & vertices() const {return m_vertices;}
  // get vector of cells
  inline std::vector<Cell> & cells() {return m_cells;}
  // get vector of cells
  inline const std::vector<Cell> & cells() const {return m_cells;}
  // true if vector of cells is empty
  bool empty() const {return m_cells.empty();}
  // get number of cells
  inline std::size_t n_cells() const {return m_cells.size();}
  // get number of active cells
  inline std::size_t n_active_cells() const { return n_cells() - m_n_split_cells - m_n_cells_with_hanging_nodes; }
  // get number of vertices
  inline std::size_t n_vertices() const {return m_vertices.size();}
  // get number of faces
  inline std::size_t n_faces() const {return m_faces.size();}

  // SEARCHING //
  // returns face index or invalid index
  size_t find_face(const std::vector<size_t> & face_vertices) const;

  // MANIPULATION //
  // tell grid which faces to split before calling split_faces method
  void mark_for_split(const std::size_t face_index);
  // split faces marked for splitting with mark_for_split
  // returns SurfaceMesh of master DFM faces
  // cleans marked_for_split array upon completion
  SurfaceMesh<double> split_faces();

  /* Split a cell by cutting it with a plane. New cell indices are appended
   * at the back, so it is safe to split multiple cells in a row.
   * Note: cell is copied since inserting new cells invalidates the pointers. */
  void split_cell(Cell cell, const angem::Plane<double> & plane,
                  const int splitting_face_marker = constants::marker_splitting_plane);
  /* split a vertex
   * retults in adding new vertices (pushed to the back of vertices set) */
  void split_vertex(const std::size_t               vertex_index,
                    const std::vector<std::size_t> &splitted_face_indices);
  // coarsen cells split built by split_cell method and restore active cells
  void coarsen_cells();

 protected:
  // find edge neighboring cells
  std::vector<size_t> neighbors_indices_(const vertex_pair & edge) const;

  // two elements are in the same group if they are neighbors and
  // the neighboring face is not in vertex_faces array
  std::vector<std::vector<std::size_t>>
  group_cells_based_on_split_faces(const std::vector<std::size_t> & affected_cells,
                                   const std::vector<std::size_t> & splitted_face_indices) const;

  /* private insert cell class that does all the cell insertion work */
  std::size_t insert_cell_(const std::vector<std::size_t> & ivertices,
                           const std::vector<size_t> take_faces,
                           const std::vector<FaceTmpData> &big_face_vector,
                           const int                        vtk_id,
                           const int                        marker);

  /* Insert an arbitrary polyhedron cell into grid .
   * A wrapper on the above function to minimize bookkeeping. */
  std::size_t insert_cell_(const std::vector<size_t> take_faces,
                           const std::vector<FaceTmpData> &big_face_vector,
                           const int                        marker = constants::default_cell_marker);
  std::size_t insert_face_(const FaceTmpData & f);

  void insert_hanging_node_(const Cell parent,
                            const vertex_pair edge,
                            const size_t inserted_vertex);
  /* if any of the cell faces contain the vertices in edge, split that face in two */
  void split_face_in_cell_(const Cell parent, const vertex_pair new_edge);
  /* get a vector of polygon global vertex indices given a vector with
   * local polygon vertex indices and a mapping vector. */
  std::vector<std::size_t>
  build_global_face_indices_(const std::vector<size_t> & polygon_local_indices,
                             const std::vector<size_t> & local_to_global) const;
  /* find out whether the face is a parent or child cell based on the
   * split information. returns pair
   * 1) true if the face and parent match identically (or new face)
   * 2) parent index. if new face then invalid_index */
  std::pair<bool,std::size_t> determine_face_parent_(const std::vector<size_t> & face_vertices,
                                                     const Cell                & parent_cell,
                                                     const std::vector<size_t> & splitting_face_vertices) const;

  std::map<vertex_pair,size_t> find_affected_edges_(const std::vector<size_t> &new_vertices,
                                                    const Cell & cell) const;

  int face_vtk_id_(const size_t n_vertices) const
  {
    switch (n_vertices)
    {
      case 3:
        return 5;  //  triangle
      case 4:
        return 9;  //  vtk_quad
      default:
        return 7;  //  vtk_polygon
    }
  }

  std::vector<size_t> sort_copy_(const std::vector<size_t>& vec) const
  {
    std::vector<size_t> copy(vec.size());
    partial_sort_copy(vec.begin(), vec.end(), copy.begin(), copy.end());
    return copy;
  }

 private:
  // ATTRIBUTES
  std::vector<angem::Point<3,double>>   m_vertices;      // vector of vertex coordinates
  std::vector<Cell>                     m_cells;
  std::vector<Face>                     m_faces;
  public:
  std::vector<std::vector<std::size_t>> m_vertex_cells;  // vertex neighboring cells
  std::vector<std::vector<std::size_t>> m_vertex_faces;  // vertex neighboring faces
  private:
  // Used as a tmp container when splitting faces for dfm
  std::vector<std::size_t> m_faces_marked_for_split;
  // to count active cells
  size_t m_n_split_cells, m_n_cells_with_hanging_nodes;
  // for keeping track of cell vertices during cell splitting
  angem::PointSet<3, double> m_vertices_from_cell_splitting;
  // for keeping track of cell vertices during cell splitting
  std::vector<size_t> m_vertices_from_cell_splitting_indices;

  friend class Subdivision;
  friend class io::VTKReader;
};


}
