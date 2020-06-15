#pragma once

#ifdef WITH_EIGEN
#include "TributaryRegion2dBase.hpp"

namespace discretization {

/** This class divides the faces of the parent PFEM cell into a numbsr of tributary regions
 * equal to the number of vertices in the face.
 */
class TributaryRegion2dVertices : public TributaryRegion2dBase
{
 public:
  TributaryRegion2dVertices(PolyhedralElementBase & element);
  virtual ~TributaryRegion2dVertices() = default;

 protected:
  void build_tributary_shapes_face_(const size_t iface, const angem::Polygon<double> & face_poly);
};

}  // end namespace discretization


#endif
