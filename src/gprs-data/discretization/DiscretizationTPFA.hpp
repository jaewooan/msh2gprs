#include "DiscretizationBase.hpp"

namespace discretization
{

class DiscretizationTPFA : public DiscretizationBase
{
 public:
  DiscretizationTPFA(const mesh::Mesh                                    & grid,
                   const std::set<int>                                 & dfm_markers,
                   const std::unordered_map<std::size_t, PhysicalFace> & dfm_faces,
                   const std::vector<std::vector<double>>              & props,
                   const std::vector<std::string>                      & keys);
  virtual void build() override;

 protected:
  void build(const mesh::const_face_iterator & face);
};

}
