#pragma once

#include "DiscretizationBase.hpp"
#include "ConnectionMap.hpp"

namespace discretization
{

class DiscretizationDFM : public DiscretizationBase
{
 public:
  DiscretizationDFM(const mesh::Mesh                              & grid,
                    const std::set<int>                           & dfm_markers,
                    const std::unordered_map<std::size_t, PhysicalFace> & dfm_faces,
                    const std::vector<std::vector<double>>        & props,
                    const std::vector<std::string>                & keys,
                    const size_t                                    shift_matrix,
                    const size_t                                    shift_dfm);

  virtual void build() override;

 protected:
  // build data like volumes, depth, poro, etc.
  virtual void build_cell_data() override;
  // build F-F and M-F connections
  void build_connections();
  // build map edge to neighboring faces
  hash_algorithms::ConnectionMap<std::vector<size_t>>
  map_edge_to_faces();
  // build matrix-fracture connection
  void build_matrix_fracture(ConnectionData & con);
  // build fracture-fracture connection
  void build_fracture_fracture(ConnectionData & con);

  // storage for the properties of dfm fractures
  const std::unordered_map<std::size_t, PhysicalFace> & dfm_faces;
  //  numbering shift of matrix CVs
  const size_t shift_matrix;
  // numbering shift of dfm CVs
  const size_t shift_dfm;
  // map cell index -> control volume
  std::vector<std::size_t> mapping;
};


} // end namespace
