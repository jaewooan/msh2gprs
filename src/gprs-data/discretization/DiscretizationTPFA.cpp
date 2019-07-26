#include "DiscretizationTPFA.hpp"
#include "angem/Point.hpp"
#include "angem/Tensor2.hpp"

namespace discretization
{

using Point = angem::Point<3,double>;
using Tensor = angem::Tensor2<3,double>;

DiscretizationTPFA::
DiscretizationTPFA(const mesh::Mesh                       & grid,
                   const std::set<int>                    & dfm_markers,
                   const std::vector<std::vector<double>> & props,
                   const std::vector<std::string>         & keys,
                   const size_t                             numbering_shift)
    :
    DiscretizationBase(grid, dfm_markers, props, keys),
    shift(numbering_shift)
{}


void DiscretizationTPFA::build()
{
  face_data.resize(grid.n_faces());
  size_t iface = 0;
  for (auto face = grid.begin_faces(); face != grid.end_faces(); ++face)
  {
    if (face.neighbors().size() == 2)
      if (!is_fracture(face.marker()))
      {
        build(face, face_data[iface]);
        iface++;
      }
  }

  DiscretizationBase::build_cell_data();
}


void DiscretizationTPFA::build(const mesh::const_face_iterator & face,
                               FaceData          & data)
{
  // Point nfmf, center0, center1, center_face;
  const auto cells = face.neighbors();
  const auto cell1 = grid.create_const_cell_iterator(cells[0]);
  const auto cell2 = grid.create_const_cell_iterator(cells[1]);

  const Point center_face = face.center();
  const Point d1 = cell1.center() - center_face;
  const Point d2 = cell2.center() - center_face;

  double dd1, dd2;   // projection of cell-face distance d1,d2 on normal
  Point Kn1, Kn2;    // perm projections onto the face normal
  double Knn1, Knn2; // // projection of conormal along normal

  const Point face_normal = face.normal();
  const Tensor K1 = get_permeability(cell1.index());
  const Tensor K2 = get_permeability(cell2.index());

  // project permeability on face normal
  Kn1 = K1*face_normal;
  Kn2 = K2*face_normal;

  // projection of distance along normal
  dd1 = fabs(face_normal.dot(d1));
  dd2 = fabs(face_normal.dot(d2));

  // projection of conormal along normal
  Knn1 = face_normal.dot(Kn1);
  Knn2 = face_normal.dot(Kn2);

  double T; // transmissibility
  double D; // inverse sum dd1 + dd2
  static const double eps = 1e-8;
  if (dd1 > eps && dd2 > eps)
  {
    T = Knn1*Knn2 / (dd1*Knn2 + dd2*Knn1 + 1.0e-35);
    D = 1.0 / (dd1 + dd2 + 1.0e-35);
  }
  else if (dd1 > eps)
  {
    T = Knn1 / dd1;
    D = 1.0 / dd1;
  }
  else if (dd2 > eps)
  {
    T = Knn2 / dd2;
    D = 1.0 / dd2;
  }
  else
  {
    T = 0;
    D = 0;
  }

  const double face_area = face.area();
  T *= face_area;
  D *= face_area;

  data.elements.resize(2);
  data.elements[0] = shift + cell1.index();
  data.elements[1] = shift + cell2.index();

  data.coefficients.resize(2);
  data.coefficients[0] = -T;
  data.coefficients[1] =  T;
}


}
