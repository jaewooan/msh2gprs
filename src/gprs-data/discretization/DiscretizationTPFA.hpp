#include "DiscretizationBase.hpp"

namespace discretization
{

class DiscretizationTPFA : public DiscretizationBase
{
 public:
  DiscretizationTPFA(const mesh::Mesh & grid);
  virtual void init();

 protected:
  void build(const mesh::const_face_iterator & face);
};

}
