#pragma once

#include "mesh/Mesh.hpp"
#include "PhysicalFace.hpp"
#include "angem/Tensor2.hpp"

namespace discretization
{

using PhysicalFace = gprs_data::PhysicalFace;

struct FaceData
{
  std::vector<double> coefficients;
  std::vector<size_t> elements;
};


struct CellData
{
  double volume;
  double porosity;
  double depth;
  std::vector<double> custom;
};

/* This is an abstract base class for
 * all discretization classes out there. */
class DiscretizationBase
{
 public:
  DiscretizationBase(const mesh::Mesh                                    & grid,
                     const std::set<int>                                 & dfm_markers,
                     const std::vector<std::vector<double>>              & props,
                     const std::vector<std::string>                      & keys);
  virtual void build() = 0;
  virtual void build_cell_data();

 protected:
  // is a face a dfm fracture
  bool is_fracture (const int marker);
  angem::Tensor2<3,double> get_permeability(const std::size_t cell) const;
  double get_porosity(const std::size_t cell) const;
  // calculate the positions of the perm keys (in keys array)
  void infer_perm_assignment();
  // calculate the positions of the porosity keys (in keys array)
  void infer_poro_assignment();
  // calculate keys that are not reserved (poro, perm)
  void infer_custom_keys();

  // ATTRIBUTES
  // input
  const mesh::Mesh & grid;
  const std::set<int> & dfm_markers;
  const std::vector<std::vector<double>> & props;
  const std::vector<std::string> & keys;

  //  computed
  std::array<int, 9> perm_keys = {-1,-1,-1,-1,-1,-1,-1,-1,-1};
  int poro_key = -1;
  std::vector<size_t> custom_keys;
  std::vector<FaceData> face_data;
  std::vector<CellData> cell_data;
};

}
