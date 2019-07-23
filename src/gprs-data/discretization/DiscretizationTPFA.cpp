#include "DiscretizationTPFA.hpp"
#include "angem/Point.hpp"

namespace discretization
{

using Point = angem::Point<3,double>;

DiscretizationTPFA::DiscretizationTPFA(const mesh::Mesh & grid)
    :
    DiscretizationBase(grid)
{}


void DiscretizationTPFA::init()
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

    // if (not fracture)
    {
    //   f->UnitNormal(nrmf);

    //   int ret1 = transposed_tensor_prod(Ktype, cell0->RealArray(K), nrmf, f1);
    }

  }
  else return;
}



}
