#pragma once

#include "mesh/Mesh.hpp"
#include "PhysicalFace.hpp"
#include "angem/Tensor2.hpp"

namespace discretization
{

using PhysicalFace = gprs_data::PhysicalFace;

/* This is an abstract base class for
 * all discretization classes out there. */
class DiscretizationBase
{
 public:
  DiscretizationBase(const mesh::Mesh                                    & grid,
                     const std::set<int>                                 & dfm_markers,
                     const std::unordered_map<std::size_t, PhysicalFace> & dfm_faces,
                     const std::vector<std::vector<double>>              & props,
                     const std::vector<std::string>                      & keys);
  virtual void build() = 0;

 protected:
  // is a face a dfm fracture
  bool is_fracture (const int marker);
  angem::Tensor2<3,double> get_permeability(const std::size_t cell) const;
  void infer_perm_assignment();

  // ATTRIBUTES
  // input
  const mesh::Mesh & grid;
  const std::set<int> & dfm_markers;
  const std::unordered_map<std::size_t, PhysicalFace> & dfm_faces;
  const std::vector<std::vector<double>> & props;
  const std::vector<std::string> & keys;

  //  computed
  std::array<int, 9> perm_keys = {-1,-1,-1,-1,-1,-1,-1,-1,-1};
};

}
