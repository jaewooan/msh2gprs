#include <const_face_iterator.hpp>
#include <mesh_methods.hpp>

namespace mesh
{

const_face_iterator::
const_face_iterator(FaceMap::const_iterator   it,
                    const angem::PointSet<3,double> & vertices)
    :
    face_it(it),
    p_mesh_vertices(&vertices)
{}


const_face_iterator::
const_face_iterator(const const_face_iterator & other)
    :
    face_it(other.face_it),
    p_mesh_vertices(other.p_mesh_vertices)
{}


const_face_iterator &
const_face_iterator::operator=(const const_face_iterator & other)
{
  face_it = other.face_it;
  p_mesh_vertices = other.p_mesh_vertices;
  return (*this);
}


bool const_face_iterator::operator==(const const_face_iterator & other) const
{
  if (face_it != other.face_it)
    return false;

  return true;
}


bool const_face_iterator::operator!=(const const_face_iterator & other) const
{
  return !(*this == other);
}


const_face_iterator & const_face_iterator::operator++()
{
  face_it++;
  return (*this);
}


int const_face_iterator::marker() const
{
  return face_it->second.marker;
}


std::vector<Point> const_face_iterator::vertices() const
{
  std::vector<std::size_t> ivertices(invert_hash(face_it->first));
  return get_vertex_coordinates(p_mesh_vertices, ivertices);
}


std::vector<std::size_t> const_face_iterator::vertex_indices() const
{
  return face_it->second.ordered_indices;
}


Point const_face_iterator::normal() const
{
  const auto poly = polygon();
  return poly.plane.normal();
}


angem::Polygon<double> const_face_iterator::polygon() const
{
  try
  {
    return angem::Polygon<double>(vertices());
  }
  catch (std::runtime_error & e)
  {
    throw std::runtime_error("invalid face " + std::to_string(index()));
  }
}


double const_face_iterator::area() const
{
  return polygon().area();
}


angem::Point<3,double> const_face_iterator::center() const
{
  const auto verts = vertex_indices();
  angem::Point<3,double> c = {0, 0, 0};
  for (const size_t vert : verts)
    c += (*p_mesh_vertices)[vert];

  for (int comp = 0; comp < 3; comp++)
    c[comp] /= verts.size();
  return c;
}


std::vector<vertex_pair> const_face_iterator::edges() const
{
  const auto & verts = face_it->second.ordered_indices;
  std::vector<vertex_pair> pairs(verts.size());
  for (std::size_t i=0; i<verts.size(); ++i)
  {
    std::size_t i1, i2;
    if (i < verts.size() - 1)
    {
      i1 = verts[i]; i2 = verts[i+1];
    }
    else
    {
      i1 = verts[ i ]; i2 = verts[ 0 ];
    }
    pairs[i] = std::make_pair(i1, i2);
  }
  return pairs;
}

}  // end namespace
