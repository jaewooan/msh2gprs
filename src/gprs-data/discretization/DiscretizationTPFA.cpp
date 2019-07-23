#include "DiscretizationTPFA.hpp"
#include "angem/Point.hpp"
#include "angem/Tensor2.hpp"

namespace discretization
{

using Point = angem::Point<3,double>;
using Tensor = angem::Tensor2<3,double>;

DiscretizationTPFA::
DiscretizationTPFA(const mesh::Mesh                                    & grid,
                   const std::set<int>                                 & dfm_markers,
                   const std::unordered_map<std::size_t, PhysicalFace> & dfm_faces,
                   const std::vector<std::vector<double>>              & props,
                   const std::vector<std::string>                      & keys)
    :
    DiscretizationBase(grid, dfm_markers, dfm_faces, props, keys)
{}


void DiscretizationTPFA::build()
{
  for (auto face = grid.begin_faces(); face != grid.end_faces(); ++face)
  {
    build(face);
  }
}


void DiscretizationTPFA::build(const mesh::const_face_iterator & face)
{
  // Point nfmf, center0, center1, center_face;
  const auto cells = face.neighbors();
  if (cells.size() == 2)
  {
    const auto cell0 = grid.create_const_cell_iterator(cells[0]);
    const auto cell1 = grid.create_const_cell_iterator(cells[1]);

    const Point center_face = face.center();
    const Point d1 = cell0.center() - center_face;
    const Point d2 = cell1.center() - center_face;

    if (!is_fracture(face.marker()))
    {
      const Point n = face.normal();
      const Tensor K = get_permeability(cell0.index());

    //   int ret1 = transposed_tensor_prod(Ktype, cell0->RealArray(K), nrmf, f1);
    }

  }
  else return;
}



}
