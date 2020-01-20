#include "DiscretizationEDFM.hpp"
#include "DiscretizationDFM.hpp"
#include "DiscretizationTPFA.hpp"
#include "angem/Projections.hpp"

namespace discretization
{

using Point = angem::Point<3,double>;
using Tensor = angem::Tensor2<3, double>;
using mesh::Face;

DiscretizationEDFM::
DiscretizationEDFM(const DoFNumbering & split_dof_numbering,
                   const DoFNumbering & combined_dof_numbering,
                   gprs_data::SimData & data,
                   std::vector<ControlVolumeData> & cv_data,
                   std::vector<ConnectionData> & connection_data,
                   const std::vector<int> & edfm_markers,
                   const EDFMMethod method)
    : m_split_dofs(split_dof_numbering),
      m_edfm_markers( edfm_markers.begin(), edfm_markers.end() ),
      m_method(method),
      DiscretizationBase(combined_dof_numbering, data, cv_data, connection_data)
{}

void DiscretizationEDFM::build()
{
  DiscretizationDFM discr_dfm(m_split_dofs, m_data, m_split_cv, m_split_con);
  discr_dfm.build();
  if (m_method != EDFMMethod::compartmental)
  {
    identify_edfm_faces_();
    if (m_edfm_faces.empty())
    {
      m_cv_data = std::move(m_split_cv);
      m_con_data = std::move(m_split_con);
      return;
    }

    build_control_volume_data_();
    build_connection_data_();

    if (m_method == EDFMMethod::projection)
      build_pedfm_();

    // copy to condata
    m_con_data.reserve(m_con_map.size());
    for (auto it = m_con_map.begin(); it != m_con_map.end(); ++it)
      m_con_data.push_back(*it);
  }
  else
  {
    m_cv_data = std::move(m_split_cv);
    m_con_data = std::move(m_split_con);
  }
}

void DiscretizationEDFM::build_control_volume_data_()
{
  for (auto & cv : m_cv_data)
  {
    cv.volume = 0;
    cv.center = {0.0, 0.0, 0.0};
    cv.aperture = 0;
  }

  // first compute parent volumes since some props are weighted by them
  m_dof_mapping.resize(m_split_cv.size());
  for (size_t i = 0; i < m_split_cv.size(); i++)
  {
    const auto &cv = m_split_cv[i];
    size_t parent_dof;
    if (cv.type == ControlVolumeType::cell)
    {
      parent_dof = m_dofs.cell_dof(cv.master);
      m_cv_data[parent_dof].master = m_grid.cell(cv.master).ultimate_parent().index();
    }
    else // if (cv.type == ControlVolumeType::face)
    {
      parent_dof = m_dofs.face_dof(cv.master);
      m_cv_data[parent_dof].master = cv.master;
    }

    m_dof_mapping[i] = parent_dof;
    m_cv_data[parent_dof].volume += cv.volume;
  }

  for (size_t i = 0; i < m_split_cv.size(); i++)
    {
      const auto &cv = m_split_cv[i];
      const size_t parent_dof = m_dof_mapping[i];
      auto &parent_cv = m_cv_data[parent_dof];
      parent_cv.type = cv.type;

      const double volume_fraction = cv.volume / parent_cv.volume;
      parent_cv.aperture += cv.aperture * volume_fraction;
      parent_cv.center += cv.center * volume_fraction;
      parent_cv.porosity += cv.porosity * volume_fraction;
      parent_cv.permeability = cv.permeability; // assume they are the same
      parent_cv.custom = cv.custom;             // assume they are the same
    }

  for(size_t i = 0; i < m_cv_data.size(); i++)
  {
    const auto & cv = m_cv_data[i];
    // std::cout << "cv.index= " << i << std::endl;
    if ( cv.type == ControlVolumeType::cell && cv.volume == 0 )
      throw std::runtime_error("discr edfm zero cell volume " +
                               std::to_string(i) + " (cell " +
                               std::to_string(cv.master) + ")");
    else if ( cv.type == ControlVolumeType::face && cv.volume == 0 )
      throw std::runtime_error("discr edfm zero face volume");
  }

}

void DiscretizationEDFM::build_connection_data_()
{
  std::vector<size_t> m_m_rebuild = create_connections_();
  // we only need to build F-M connections for both dfms and edfms
  // since the frac dofs do not get merged
  for (const size_t icon : m_m_rebuild)
  {
    auto &con = m_con_map.get_data(icon);
    // std::cout << "rebuild = :" ;
    // for (auto i : con.elements) std::cout << m_cv_data[i].master << " ";
    // for (auto i : con.elements) std::cout << i << " ";
    // std::cout << std::endl;

    con.center /= con.area;
    assert( con.elements[0] < m_cv_data.size() );
    assert( con.elements[1] < m_cv_data.size() );
    DiscretizationTPFA::build_mo(con, m_cv_data[con.elements[0]], m_cv_data[con.elements[1]]);
  }

  for (auto con = m_con_map.begin(); con != m_con_map.end(); ++con)
    if (con->type == ConnectionType::matrix_fracture)
    {
      if (find_edfm_elements_(*con).empty())
        build_matrix_dfm_(*con);
      else build_matrix_edfm_(*con);
    }
}

void DiscretizationEDFM::build_matrix_dfm_(ConnectionData &con)
{
  assert(con.elements.size() == 2);
  const auto & cell = m_cv_data[con.elements[1]];
  const auto & frac = m_cv_data[con.elements[0]];
  assert ( cell.type == ControlVolumeType::cell );
  assert ( frac.type == ControlVolumeType::face );
  DiscretizationDFM::build_matrix_fracture(con, frac, cell);
}

void DiscretizationEDFM::build_matrix_edfm_(ConnectionData &con)
{
  assert(con.elements.size() == 2);
  const auto & cell = m_cv_data[con.elements[1]];
  const auto & frac = m_cv_data[con.elements[0]];
  assert ( cell.type == ControlVolumeType::cell );
  assert ( frac.type == ControlVolumeType::face );
  const double k_d = con.normal * (cell.permeability * con.normal);  // directional permeability
  const double T = k_d * con.area / con.distances[1]; // average distance from cell to fracture
  con.coefficients = {-T, T};
}

void DiscretizationEDFM::identify_edfm_faces_()
{
  for (auto face = m_grid.begin_active_faces(); face != m_grid.end_active_faces(); ++face)
  {
    if (m_edfm_markers.find(face->marker()) != m_edfm_markers.end())
      m_edfm_faces.insert(m_split_dofs.face_dof(face->index()));
  }
}

std::vector<size_t> DiscretizationEDFM::find_edfm_elements_(const ConnectionData & con)
{
  std::vector<size_t> result;
  for (std::size_t i=0; i<con.elements.size(); ++i)
    if (m_edfm_faces.find(con.elements[i]) != m_edfm_faces.end())
      result.push_back(i);
  return result;
}

std::vector<size_t> DiscretizationEDFM::create_connections_()
{
  std::unordered_set<size_t> rebuild_connections;
  for (const auto &con : m_split_con)
  {
    assert( con.elements[0] < m_dof_mapping.size() );
    assert( con.elements[1] < m_dof_mapping.size() );
    const size_t dof1 = m_dof_mapping[con.elements[0]];
    const size_t dof2 = m_dof_mapping[con.elements[1]];

    if (dof1 != dof2)
    {
      if (!m_con_map.contains(dof1, dof2))
      {
        // std::cout << "insert " << m_cv_data[dof1].master << " " <<m_cv_data[dof2].master
        //     << " (" << dof1 << " " << dof2 << ") "<< std::endl;
        m_con_map.insert(dof1, dof2);
      }
      if (con.type == ConnectionType::matrix_matrix)
      {
        if ((dof1 == 458 or dof1 == 731) and (dof2 == 458 or dof2 == 731))
        {
          std::cout << "here" << std::endl;
          std::cout << dof1 << " " << dof2 << std::endl;
        }
      }

      auto &new_con = m_con_map.get_data(dof1, dof2);
      new_con.elements = {dof1, dof2};
      new_con.type = con.type;

      if ((con.type == ConnectionType::fracture_fracture))
      {
        new_con = con;
      }
      else if (con.type == ConnectionType::matrix_fracture)
      {
        const auto &cell = m_split_cv[con.elements[1]];
        const auto &parent_cell = m_cv_data[dof2];
        const auto &face = m_split_cv[con.elements[0]];
        new_con.area += 2 * con.area;
        new_con.normal = con.normal;
        const double dist = (cell.center - face.center).dot(con.normal);
        const double cell_volume_ratio = cell.volume / parent_cell.volume;
        new_con.distances.resize(2);
        new_con.distances[0] += 0.0; // from fracture to connection
        new_con.distances[1] += std::fabs(dist) * cell_volume_ratio; // from cell to connection
      }
      else // if (con.type == matrix_matrix)
      {
        bool copy_connection = true;
        for (const size_t cv : con.elements)
          if (m_grid.cell(m_split_cv[cv].master).ultimate_parent() != m_grid.cell(m_split_cv[cv].master))
            copy_connection = false;
        if (copy_connection)
        {
          new_con = con;
          new_con.elements = {dof1, dof2};
        }
        else
        {
          rebuild_connections.insert(m_con_map.index(dof1, dof2));
          new_con.area += con.area;
          new_con.normal = con.normal;
          new_con.center += con.center * con.area;
        }
      }  // end M-M
    }    // end insertion
  }      //end con loop
  return std::vector<size_t>(rebuild_connections.begin(), rebuild_connections.end());
}

void DiscretizationEDFM::build_pedfm_()
{
  PureConnectionMap cleared_connections;
  for (auto face = m_grid.begin_active_faces(); face != m_grid.end_active_faces(); ++face)
  {
    if (m_edfm_markers.find(face->marker()) != m_edfm_markers.end())
    {
      const auto & frac_cell = face->neighbors()[0]->ultimate_parent();
      // std::cout << "searching nighbors for frac cell " << frac_cell.index() << std::endl;

      assert(!pedfm_select_faces_(*face).empty());
      for (const mesh::Face* face2 : pedfm_select_faces_(*face))
      {
        const size_t iother_cell = pedfm_find_other_cell_(*face, *face2);

        // determing all participating dofs
        const size_t cell_dof1 = m_dofs.cell_dof(frac_cell.index());
        const size_t cell_dof2 = m_dofs.cell_dof(iother_cell);
        const size_t face_dof = m_dofs.face_dof(face->index());

        // if already cleared, then skip
        if (cleared_connections.contains(cell_dof1, cell_dof2)) continue;
        // this happens when several edfm's in a cell
        if (iother_cell == frac_cell.index()) continue;
        // have connection between two parent cells
        // std::cout << "\tfrac " << frac_cell.index() << " neib " << iother_cell << std::endl;

        // {
        //   std::cout << "cell dof1 = " << cell_dof1 << std::endl;
        //   std::cout << "cell dof2 = " << cell_dof2 << std::endl;
        //   std::cout << "face_dof = " << face_dof << std::endl;
        // }

        assert( m_con_map.contains( cell_dof1, cell_dof2 ) );

        // projection of frac face onto M-M connecting face
        const auto projection_points = angem::project(face->polygon(), face2->polygon());
        if (projection_points.size() < 3) continue;  // no real projection

        const size_t ifrac_mat_con = m_con_map.insert(face_dof, cell_dof2);
        auto & frac_mat_con = m_con_map.get_data(ifrac_mat_con);
        frac_mat_con.elements = { face_dof, cell_dof2 };
        // note that we invoke it after the insertion insertion killls references
        auto & mat_mat_con = m_con_map.get_data(cell_dof1, cell_dof2);

        const angem::Polygon<double> projection(projection_points);
        frac_mat_con.normal = projection.normal();
        frac_mat_con.area = projection.area();
        build_pedfm_(mat_mat_con, frac_mat_con);
        // if small trans then kill connection
        if (std::fabs(mat_mat_con.coefficients[0]) < 1e-6 * std::fabs(frac_mat_con.coefficients[0]))
        {
          std::cout << "clear connection "<< frac_cell.index() << " " << iother_cell << std::endl;
          cleared_connections.insert(cell_dof1, cell_dof2);
        }

      }
    }
  }

  // clear connections
  for (auto it = cleared_connections.begin(); it!=cleared_connections.end(); ++it)
  {
    const auto pair_elements = it.elements();
    m_con_map.remove(pair_elements.first, pair_elements.second);
  }
}

std::vector<const mesh::Face*> DiscretizationEDFM::pedfm_select_faces_(const mesh::Face & frac_face) const
{
  //  select largest cell neighbor
  const auto & neighbors = frac_face.neighbors();
  const mesh::Cell* smallest_neighbor;
  if (neighbors[0]->volume() > neighbors[1]->volume())
    smallest_neighbor = neighbors[1];
  else smallest_neighbor = neighbors[0];

  std::vector<const Face*> result;
  for (auto face : smallest_neighbor->faces())
  {
    if (*face == frac_face) continue; // same face
    if (face->neighbors().size() < 2) continue;  // skip boundary
    if (std::fabs(face->normal().dot(frac_face.normal())) < 1e-6) continue;  // face ⟂ frac
    if (m_dofs.is_active_face(face->index())) continue; // skip frac faces

    // if (*face == frac_face)
    // {
    //   std::cout << "sam face" << std::endl; continue;
    // }
    // if (face->neighbors().size() < 2)
    // {
    //   std::cout << face->index() << " boundary face" << std::endl;
    //   continue; // skip boundary
    // }
    // if (std::fabs(face->normal().dot(frac_face.normal())) < 1e-6)
    // {
    //   std::cout << face->index() << " perpendicular" << std::endl;
    //   continue; // face ⟂ frac
    // }
    // if (m_dofs.is_active_face(face->index()))
    // {
    //   std::cout << face->index() << " dfm face" << std::endl;
    //   continue; // skip frac faces
    // }

    result.push_back(&*face);
  }
  return result;
}

size_t DiscretizationEDFM::pedfm_find_other_cell_(const Face & frac, const Face & other) const
{
  const auto & frac_cell = frac.neighbors()[0]->ultimate_parent();
  for (const auto & cell : other.neighbors())
  {
    const auto & parent = cell->ultimate_parent();
    if (frac_cell != parent)
      return parent.index();
  }
  return frac_cell.index();
}

void DiscretizationEDFM::build_pedfm_(ConnectionData & mm_con,
                                      ConnectionData & fm_con)
{
  assert(mm_con.elements.size() == 2);
  assert(fm_con.elements.size() == 2);
  auto pair_cells = find_fracture_cv_and_nonfracture_cv_(mm_con,fm_con);
  assert( pair_cells.first < m_cv_data.size() );
  assert( pair_cells.second < m_cv_data.size() );
  assert( fm_con.elements[0] < m_cv_data.size() );

  const auto & cell1 = m_cv_data[pair_cells.first];
  const auto & cell2 = m_cv_data[pair_cells.second];
  const auto & frac = m_cv_data[fm_con.elements[0]];

  // first modify old M_M connection
  {
    const Point &c1 = cell1.center;
    const Point &c2 = cell2.center;
    const Point &cf = mm_con.center;
    const Point &n = mm_con.normal;
    const double t = (cf - c1).dot(n) / (c2 - c1).dot(n);
    const Point cp = c1 + t * (c2 - c1);
    const Tensor &K1 = cell1.permeability;
    const Tensor &K2 = cell2.permeability;
    const double Kp1 = (K1 * (c1 - cp).normalize()).norm();
    const double Kp2 = (K2 * (c2 - cp).normalize()).norm();
    const double subtracted_area = fm_con.area;
    const double dT1 = subtracted_area * Kp1 / (c1 - cp).norm();
    const double dT2 = subtracted_area * Kp2 / (c2 - cp).norm();
    const double dT = dT1 * dT2 / (dT1 + dT2);
    const double new_mm_trans =
        std::min(std::fabs(mm_con.coefficients[0]) - dT, 0.0);
    mm_con.coefficients = {-new_mm_trans, new_mm_trans};
  }

  // now compute non-neighboring F-M trans
  {
    const Point &c1 = frac.center;
    const Point &c2 = cell2.center;
    const Point &cf = fm_con.center;
    const Point &n = fm_con.normal;
    // projection point
    const double t =  (cf - c1).dot(n) / (c2 - c1).dot(n);
    const Point cp = c1 + t*(c2 - c1);
    // project permeability
    const Tensor & K2 = cell2.permeability;
    const double & Kp1 = frac.permeability(0, 0);;
    const double Kp2 = (K2 * (c2 - cp).normalize()).norm();
    // cell-face transmissibility
    const double face_area = fm_con.area;
    const double T1 = face_area * Kp1 / (c1 - cp).norm();
    const double T2 = face_area * Kp2 / (c2 - cp).norm();
    // face transmissibility
    const double T = T1*T2 / ( T1 + T2 );
    fm_con.coefficients = {-T, T};
  }

}

std::pair<size_t,size_t> DiscretizationEDFM::find_fracture_cv_and_nonfracture_cv_(const ConnectionData & mm_con, const ConnectionData & fm_con) const
{
  if ( fm_con.elements[1] == mm_con.elements[0] )
    return {mm_con.elements[0], mm_con.elements[1]};
  else return {mm_con.elements[1], mm_con.elements[0]};
}

}  // end namespace discretization
