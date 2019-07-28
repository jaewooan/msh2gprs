#pragma once
#include "angem/Point.hpp"
#include "BoundaryConditionType.hpp"
#include <vector>

namespace gprs_data
{

struct FractureFace
{
  std::size_t face_index;                // face index
  int marker;                       // face marker
  int fracture_index;               // dfm fracture index
  std::size_t cv_index;             // control volume index
  std::pair<std::size_t,std::size_t> neighbor_cells;
  double aperture;
  double conductivity;
};


/* That shit that requires mesh splitting */
struct MechanicalFractureFace
{
  bool coupled;
  std::pair<std::size_t,std::size_t> neighbor_cells;
};


/* Used for geomechanics boundary conditions */
struct BoundaryFace
{
  BoundaryConditionType ntype;      // dirichlet or neumann
  std::size_t face_index;           // face index
  int marker;                       // face marker
  angem::Point<3,double> condition; // displacement or force
  // can still be more than 
  std::vector<std::size_t> neighbor_cells;
};

}
