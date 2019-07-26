#include "DiscretizationBase.hpp"

namespace discretization
{

class DiscretizationTPFA : public DiscretizationBase
{
 public:
  DiscretizationTPFA(const mesh::Mesh                     & grid,
                   const std::set<int>                    & dfm_markers,
                   const std::vector<std::vector<double>> & props,
                   const std::vector<std::string>         & keys,
                   const size_t                             numbering_shift = 0);

  virtual void build() override;

 protected:
  void build(const mesh::const_face_iterator & face,
             FaceData          & data);

  // shift of controle volume indices
  // used i.e. when grid is a subdomain or when
  // domain control volumes follow fracture control volumes
  // in numbering
  const size_t shift;
};

}
