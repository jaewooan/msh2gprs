#include "DiscretizationDFM.hpp"
#include "angem/Tensor2.hpp"
#include <cmath>  // std::isnan

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


void DiscretizationDFM::build()
{
  std::cout << "build cv data" << std::endl;
  build_cell_data();
  std::cout << "build connection list" << std::endl;
  // build connection lists (no data)
  build_fracture_fracture_connections();
  build_fracture_matrix_connections();
  // build discretization

  std::cout << "build face discretization" << std::endl;
  for (auto & con : con_data)
    switch (con.type) {
      case ConnectionType::fracture_fracture :
        {
          std::cout << "F-F" << std::endl;
          build_fracture_fracture(con);
          break;
        }
      case ConnectionType::matrix_fracture :
        {
          std::cout << "M-F" << std::endl;
          build_matrix_fracture(con);
          break;
        }
      default:
        throw std::runtime_error("invalid connection type");
    }
}


void DiscretizationDFM::build_cell_data()
{
  // could be less though, just to be save
  mapping.resize(grid.n_cells());
  cv_data.resize(dfm_faces.size() * 3);
  std::unordered_set<size_t> bounding_cells;

  for (auto face = grid.begin_faces(); face != grid.end_faces(); ++face)
    if (face.neighbors().size() == 2)
    {
      if (is_fracture(face.marker()))
      {
        const auto it = dfm_faces.find(face.index());
        assert(it != dfm_faces.end());
        const auto & face_props = it->second;

        auto & data = cv_data[face_props.cv_index];
        data.type = ControlVolumeType::face;
        std::cout << "dfm face " << face_props.cv_index << std::endl;
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
    }

  size_t cv = dfm_faces.size();
  for (const std::size_t i : bounding_cells)
  {
    std::cout << i << std::endl;
    const auto cell = grid.create_const_cell_iterator(i);
    auto & data = cv_data[cv];
    // auto & data = cv_data.emplace_back();
    mapping[i] = cv;

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

  cv_data.resize(cv);
}


hash_algorithms::ConnectionMap<std::vector<size_t>>
DiscretizationDFM::map_edge_to_faces()
{
  hash_algorithms::ConnectionMap<std::vector<size_t>> edge_face_connections;
  for (auto face = grid.begin_faces(); face != grid.end_faces(); ++face)
    if (face.neighbors().size() == 2)
      if (is_fracture(face.marker()))
      {
        const auto it = dfm_faces.find(face.index());
        assert(it != dfm_faces.end());
        const auto & props = it->second;
        auto perm = Tensor::make_unit_tensor();
        perm *= it->second.conductivity / it->second.aperture;

        for (const auto & edge : face.edges())
        {
          size_t index;
          if (edge_face_connections.has(edge.first, edge.second))
          {
            index = edge_face_connections.index(edge.first, edge.second);
          }
          else
          {
            index = edge_face_connections.insert(edge.first, edge.second);
          }

          auto & cvs = edge_face_connections.get_data(index);

          // get frac properties
          cvs.push_back(face.index());
        }
      }
  return edge_face_connections;
}


void DiscretizationDFM::build_fracture_fracture(ConnectionData & con)
{
  std::cout << "not implemented" << std::endl;
  abort();
}


void DiscretizationDFM::build_matrix_fracture(ConnectionData & con)
{
  // cause I built them that way
  const auto & cv_frac = cv_data[con.elements[0]];
  const auto & cv_cell = cv_data[con.elements[1]];
  assert(cv_frac.type == ControlVolumeType::face);
  assert(cv_cell.type == ControlVolumeType::cell);
  assert(con.type == ConnectionType::matrix_fracture);
  std::cout << "got data" << std::endl;

  // project cell permeability
  const auto f = cv_cell.center - con.center;
  std::cout << "accessed data" << std::endl << std::flush;
  std::cout << "f = " << f << std::endl;
  const double K_cell = (cv_cell.permeability * (f/f.norm())).norm();
  // cause that's how Mo did it ¯\_(ツ)_/¯
  const double K_frac = cv_frac.permeability(0, 0);
  std::cout << "K_frac = " << K_frac << std::endl;

  const double T_cell = con.area * K_cell / f.norm();
  const double T_face = con.area * K_frac;

  // connection transmissibility
  double T = 0;
  if ( !std::isnan(1. / (T_cell + T_face) ) )
    T = T_cell*T_face / (T_cell + T_face);

  con.coefficients.resize(2);
  con.coefficients[0] = -T;
  con.coefficients[1] =  T;
}


void DiscretizationDFM::build_fracture_matrix_connections()
{
  for (auto face = grid.begin_faces(); face != grid.end_faces(); ++face)
  {
    if (is_fracture(face.marker()))
    {
      const auto & neighbors = face.neighbors();
      // find frac cv index
      std::cout << face.index() << std::endl;
      const auto it = dfm_faces.find(face.index());
      assert(it != dfm_faces.end());
      const std::size_t cv_frac = it->second.cv_index;

      const auto & face_props = it->second;
      // connection fracture-cell1
      {
        con_data.emplace_back();
        auto &con = con_data.back();
        con.elements.push_back(cv_frac);
        con.elements.push_back(mapping[neighbors[0]]);
        con.type = ConnectionType::matrix_fracture;
        con.area = face.area();
        con.normal = face.normal();
        con.center = face.center();
      }

      //  connection fracture-cell2
      {
        con_data.emplace_back();
        auto &con = con_data.back();
        con.type = ConnectionType::matrix_fracture;
        con.elements.push_back(cv_frac);
        con.elements.push_back(mapping[neighbors[1]]);
        con.area = face.area();
        con.normal = face.normal();
        con.center = face.center();
      }
    }
  }
}


void DiscretizationDFM::build_fracture_fracture_connections()
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
      std::cout << "shitty shit" << std::endl;
      auto &con = con_data.emplace_back();
      con.elements = edge_faces;
      con.type = ConnectionType::fracture_fracture;
    }

  std::cout << "bye bye" << std::endl;
  exit(0);
  
}

} // end namespace
