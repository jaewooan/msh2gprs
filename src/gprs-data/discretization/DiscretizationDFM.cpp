#include "DiscretizationDFM.hpp"
#include "angem/Tensor2.hpp"

namespace discretization
{

using Point = angem::Point<3,double>;
using Tensor = angem::Tensor2<3,double>;

DiscretizationDFM::
DiscretizationDFM(const mesh::Mesh                              & grid,
                  const std::set<int>                           & dfm_markers,
                  const std::unordered_map<std::size_t, PhysicalFace> & dfm_faces,
                  const std::vector<std::vector<double>>        & props,
                  const std::vector<std::string>                & keys,
                  const size_t                                    shift_matrix,
                  const size_t                                    shift_dfm)
: DiscretizationBase(grid, dfm_markers, props, keys),
  dfm_faces(dfm_faces), shift_matrix(shift_matrix),
  shift_dfm(shift_dfm)
{
  // check that ranges don't intersect
  assert( (shift_matrix + grid.n_cells() <= shift_dfm) ||
          (shift_matrix >= shift_dfm + dfm_faces.size()) );
}


// size_t DiscretizationDFM::matrix_cv(const size_t cell)
// {

// }


void DiscretizationDFM::build()
{
  build_cell_data();
  build_connections();

  for (auto & con : con_data)
    switch (con.type) {
      case ConnectionType::fracture_fracture :
        build_fracture_fracture(con);
        break;
      case ConnectionType::matrix_fracture :
        build_matrix_fracture(con);
        break;
      default:
        throw std::runtime_error("invalid connection type");
    }
}


void DiscretizationDFM::build_connections()
{
  // build fracture-fracture-connections
  // map edges to dfm faces
  const auto edge_face_connections = map_edge_to_faces();
  // we won't need all of those since there are edges connected to
  // only one face
  con_data.reserve(2*dfm_faces.size() + edge_face_connections.size());
  for (const auto & edge_faces : edge_face_connections)
    if (edge_faces.size() > 1)
    {
      auto &con = con_data.emplace_back();
      con.elements = edge_faces;
      con.type = ConnectionType::fracture_fracture;
    }

  // build fracture_matrix connections
  for (auto face = grid.begin_faces(); face != grid.end_faces(); ++face)
    if (is_fracture(face.index()))
    {
      auto &con = con_data.emplace_back();
      con.elements = face.neighbors();
      con.type = ConnectionType::matrix_fracture;
    }
}


void DiscretizationDFM::build_cell_data()
{
  // could be less though, just to be save
  mapping.resize(grid.n_cells());
  cv_data.reserve(dfm_faces.size() * 3);
  std::unordered_set<size_t> bounding_cells;

  for (auto face = grid.begin_faces(); face != grid.end_faces(); ++face)
    if (face.neighbors().size() == 2)
      if (is_fracture(face.index()))
      {
        const auto it = dfm_faces.find(face.index());
        assert(it == dfm_faces.end());
        const auto & face_props = it->second;

        auto & data = cv_data[face_props.cv_index];
        data.type = ControlVolumeType::face;
        data.master = face.index();
        data.center = face.center();
        data.volume = face_props.aperture * face.area();
        data.permeability = Tensor::make_unit_tensor();
        data.permeability *= (face_props.conductivity / face_props.aperture);
        data.porosity = 1.0;

        // take custom properties as weighted average from bounding cells
        const auto cells = face.neighbors();
        const auto cell1 = grid.create_const_cell_iterator(cells[0]);
        const auto cell2 = grid.create_const_cell_iterator(cells[1]);
        data.custom.resize(custom_keys.size());
        for (size_t j = 0; j < custom_keys.size(); ++j)
          data.custom[j] +=
              (props[cell1.index()][custom_keys[j]] * cell1.volume() +
               props[cell2.index()][custom_keys[j]] * cell2.volume()) /
              ( cell1.volume() + cell2.volume() ) ;

        bounding_cells.insert(cells[0]);
        bounding_cells.insert(cells[1]);
      }

  size_t cv = dfm_faces.size();
  for (const std::size_t i : bounding_cells)
  {
    const auto cell = grid.create_const_cell_iterator(i);
    auto & data = cv_data[cv];
    data.type = ControlVolumeType::cell;
    data.master = i;
    data.center = cell.center();
    data.volume = cell.volume();
    data.permeability = get_permeability(i);

    data.custom.resize(custom_keys.size());
    for (size_t j = 0; j < custom_keys.size(); ++j)
      data.custom[j] = props[i][custom_keys[j]];

    cv++;
  }
}


hash_algorithms::ConnectionMap<std::vector<size_t>>
DiscretizationDFM::map_edge_to_faces()
{
  hash_algorithms::ConnectionMap<std::vector<size_t>> edge_face_connections;
  for (auto face = grid.begin_faces(); face != grid.end_faces(); ++face)
    if (face.neighbors().size() == 2)
      if (is_fracture(face.index()))
      {
        const auto it = dfm_faces.find(face.index());
        assert(it == dfm_faces.end());
        const auto & props = it->second;
        auto perm = Tensor::make_unit_tensor();
        perm *= it->second.conductivity / it->second.aperture;

        for (const auto & edge : face.edges())
        {
          size_t index;
          if (edge_face_connections.has(edge.first, edge.second))
            index = edge_face_connections.index(edge.first, edge.second);
          else
            index = edge_face_connections.insert(edge.first, edge.second);

          auto & cvs = edge_face_connections.get_data(index);

          // get frac properties
          cvs.push_back(face.index());
        }
      }
  return edge_face_connections;
}


void DiscretizationDFM::build_fracture_fracture(ConnectionData & con)
{

}


void DiscretizationDFM::build_matrix_fracture(ConnectionData & con)
{

}

} // end namespace
