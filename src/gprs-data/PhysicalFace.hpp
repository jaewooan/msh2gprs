#pragma once
#include "angem/Point.hpp"
#include <vector>

namespace gprs_data
{

struct PhysicalFace
{
  int ntype;
  std::size_t nface;
  int nmark;
  int ifracture;
  std::size_t nfluid;
  bool coupled;
  angem::Point<3,double> condition;
  std::vector<std::size_t> neighbor_cells;
  double aperture;
  double conductivity;
};


}
