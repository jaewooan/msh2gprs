#pragma once

#include "mesh/Mesh.hpp"

namespace discretization
{

/* This is an abstract base class for
 * all discretization classes out there. */
class DiscretizationBase
{
 public:
  DiscretizationBase(const mesh::Mesh & grid);
  virtual void init() = 0;

 protected:
  const mesh::Mesh & grid;
};

}
