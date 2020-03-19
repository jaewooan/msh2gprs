#pragma once

#include "discretization/flow/DoFNumbering.hpp"
#include "mesh/Mesh.hpp"

namespace gprs_data {

using discretization::DoFNumbering;

/**
 * This class implements the numbering of active grid cells.
 * Currently, mesh::Mesh does not keep track of the active cells
 * numbers, insead using cell indices that may be not consequtive if
 * some of the cells are refined.
 * For this purpose I implemented this simple class
 */
class GridCellNumberingManager
{
 public:
  GridCellNumberingManager(const mesh::Mesh & grid);
  DoFNumbering * get_cell_numbering();

 private:
  const mesh::Mesh & _grid;
};

}  // end namespace gprs_data
