#pragma once

#include <vector>

namespace discretization {

struct FEPointData
{
  std::vector<double> values;
  std::vector<angem::Point<3,double>> grads;
  double weight;
};

struct FiniteElementData
{
  std::vector<FEPointData> points;
};

}  // end namespace discretization