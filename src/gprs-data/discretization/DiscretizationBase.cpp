#include "DiscretizationBase.hpp"

namespace discretization
{

DiscretizationBase::
DiscretizationBase(const mesh::Mesh                                    & grid,
                   const std::set<int>                                 & dfm_markers,
                   const std::unordered_map<std::size_t, PhysicalFace> & dfm_faces,
                   const std::vector<std::vector<double>>              & props,
                   const std::vector<std::string>                      & keys)
    : grid(grid),
      dfm_markers(dfm_markers),
      dfm_faces(dfm_faces),
      props(props),
      keys(keys)
{
  infer_perm_assignment();
}


bool DiscretizationBase::is_fracture(const int marker)
{
  const auto it = dfm_markers.find(marker);
  if (it != dfm_markers.end())
    return true;
  else return false;
}


angem::Tensor2<3,double>
DiscretizationBase::get_permeability(const std::size_t cell) const
{
  assert(cell < props.size());
  angem::Tensor2<3,double> K;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      K(i, j) = (perm_keys[3*i+j] >= 0) ? props[cell][perm_keys[3*i+j]] : 0;
  return K;
}


void DiscretizationBase::infer_perm_assignment()
{
  bool found_perm_x = false;
  bool found_perm_y = false;
  bool found_perm_z = false;
  bool found_perm = false;
  for (std::size_t i = 0; i < keys.size(); i++)
  {
    const auto & key = keys[i];
    if (key == "PERMX") perm_keys[0] = i;
    if (key == "PERMY") perm_keys[1*3 + 1] = i;
    if (key == "PERMZ") perm_keys[2*3 + 2] = i;
    if (key == "PERM")
    {
      perm_keys[0] = i;
      perm_keys[1*3 + 1] = i;
      perm_keys[2*3 + 2] = i;
    }
  }

  throw std::invalid_argument("permebility is undefined");
  return;
}

}
