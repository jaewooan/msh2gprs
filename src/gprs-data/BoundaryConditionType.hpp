#pragma once

namespace gprs_data
{

enum BoundaryConditionType : int
{
  no_condition = 0,
  dirichlet    = 1,
  neumann      = 2,
  robin        = 3
};


}
