// Copyright (c) 2014, 2017 INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// Author(s)     : Stephen Kiazyk,
//                 Mael Rouxel-Labbé
//
#ifndef CGAL_POLYGON_MESH_PROCESSING_LOCATE_H
#define CGAL_POLYGON_MESH_PROCESSING_LOCATE_H

#include <CGAL/license/Polygon_mesh_processing/locate.h>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/array.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/Default.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/Random.h>

#include <boost/foreach.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/variant.hpp>

#include <iostream>
#include <utility>

// For the complete file:
// \details If `tm` is the input graph and given the pair (`f`, `bc`)
//          such that `bc` is `(w0, w1, w2)`, the correspondance with the weights in `bc`
//          and the vertices of the face `f` is the following:
//          - `w0 = source(halfedge(f, tm), tm)`
//          - `w1 = target(halfedge(f, tm), tm)`
//          - `w2 = target(next(halfedge(f, tm), tm), tm)`

// The following overloads are offered:
// - Face_location locate(vertex_descriptor, TriangleMesh)
// - Face_location locate(halfedge_descriptor, FT, TriangleMesh)
// - Face_location locate(face_descriptor, Point, TriangleMesh, NamedParameters)
// - Face_location locate(face_descriptor, Point, TriangleMesh)
// - Face_location locate(Face_location, Point, TriangleMesh)
// - Face_location locate(Point, AABBtree, TriangleMesh, NamedParameters)
// - Face_location locate(Point, AABBtree, TriangleMesh)
// - Face_location locate(Point, TriangleMesh, NamedParameters)
// - Face_location locate(Point, TriangleMesh)

// @todo check that all functions compile (probably need to use sfinae a few more times)

namespace CGAL{

namespace Polygon_mesh_processing {

namespace internal {

// to not carry 5-lines-long types everywhere
template <typename TriangleMesh,
          typename NamedParameters = Default>
struct Locate_types
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor      face_descriptor;
  typedef boost::variant<vertex_descriptor,
                         halfedge_descriptor,
                         face_descriptor>                                  descriptor_variant;

  typedef typename boost::mpl::if_<
                     boost::is_same<
                       NamedParameters, Default>,
                       typename boost::property_map<TriangleMesh,
                                                    CGAL::vertex_point_t>::type,
                       typename GetVertexPointMap<TriangleMesh,
                                                  NamedParameters>::type
                     >::type                                               VertexPointMap;
  typedef typename boost::property_traits<VertexPointMap>::value_type      Point;

  typedef typename CGAL::Kernel_traits<Point>::type                        Kernel;
  typedef typename Kernel::FT                                              FT;

  typedef CGAL::cpp11::array<FT, 3>                                        Barycentric_coordinates;
  typedef std::pair<face_descriptor, Barycentric_coordinates>              Face_location;
};

} // namespace internal

// forward declarations
template <typename TriangleMesh>
typename internal::Locate_types<TriangleMesh>::descriptor_variant
get_descriptor_from_location(const typename internal::Locate_types<TriangleMesh>::Face_location& loc,
                             const TriangleMesh& tm);

template <typename TriangleMesh>
typename internal::Locate_types<TriangleMesh>::Face_location
locate(typename boost::graph_traits<TriangleMesh>::halfedge_descriptor he,
       typename internal::Locate_types<TriangleMesh>::FT t,
       const TriangleMesh& tm);
// end of forward declarations

namespace internal {

template<typename TriangleMesh, typename Container>
void
incident_faces(const typename internal::Locate_types<TriangleMesh>::Face_location& location,
               Container& c,
               const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor      face_descriptor;

  typedef boost::variant<vertex_descriptor,
                         halfedge_descriptor,
                         face_descriptor>                                  descriptor_variant;

  descriptor_variant dv = get_descriptor_from_location(location, tm);

  if(dv.which() == 0)
  {
    vertex_descriptor vd = boost::get<vertex_descriptor>(dv);
    //@todo replace by CGAL_FOREACH when it's available
    BOOST_FOREACH(face_descriptor fd, CGAL::faces_around_target(halfedge(vd, tm), tm))
    {
      c.insert(fd);
    }
  }
  else if(dv.which() == 1)
  {
    halfedge_descriptor hd = boost::get<halfedge_descriptor>(dv);
    c.insert(face(hd, tm));
    c.insert(face(opposite(hd, tm), tm));
  }
  else // dv.which() == 2
  {
    CGAL_assertion(dv.which() == 2);
    face_descriptor fd = boost::get<face_descriptor>(dv);
    c.insert(fd);
  }
}

// Snapping coordinates for robustness
//@todo should this snap to closest border without any tolerance ?
template<typename TriangleMesh>
void
snap_coordinates_to_border(typename Locate_types<TriangleMesh>::Barycentric_coordinates& coords,
                           const typename Locate_types<TriangleMesh>::FT tolerance =
                             std::numeric_limits<typename Locate_types<TriangleMesh>::FT>::epsilon())
{
  typedef typename internal::Locate_types<TriangleMesh>::FT            FT;

  std::cout << "Pre-snapping: " << coords[0] << " "
                                << coords[1] << " "
                                << coords[2] << std::endl;
  std::cout << "Sum: " << coords[0] + coords[1] + coords[2] << std::endl;

  std::cout << "tolerance: " << tolerance << std::endl;

  // To still keep a sum roughly equals to 1, keep in memory the small changes
  FT residue = 0.;

  for(int i=0; i<3; ++i)
  {
    if(std::abs(coords[i]) <= tolerance)
    {
      residue += coords[i];
      coords[i] = 0.;
    }
    else if(std::abs(1 - coords[i]) <= tolerance)
    {
      residue -= 1. - coords[i];
      coords[i] = 1.;
    }
  }

  // Dump the residue into one of the barycentric values that is neither 0 nor 1
  for(int i=0; i<3; ++i)
  {
    if(coords[i] != 0. && coords[i] != 1.)
    {
      coords[i] += residue;
      break;
    }
  }

  std::cout << "Post-snapping: " << coords[0] << " "
                                 << coords[1] << " "
                                 << coords[2] << std::endl;
  std::cout << "Sum: " << coords[0] + coords[1] + coords[2] << std::endl;
}

template<typename TriangleMesh>
void
snap_location_to_border(typename Locate_types<TriangleMesh>::Face_location& loc,
                        const typename Locate_types<TriangleMesh>::FT tolerance =
                          std::numeric_limits<typename Locate_types<TriangleMesh>::FT>::epsilon())
{
  return snap_coordinates_to_border<TriangleMesh>(loc.second, tolerance);
}

template<typename TriangleMesh>
boost::optional<typename boost::graph_traits<TriangleMesh>::halfedge_descriptor>
common_halfedge(const typename boost::graph_traits<TriangleMesh>::face_descriptor first_fd,
                const typename boost::graph_traits<TriangleMesh>::face_descriptor second_fd,
                const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor  halfedge_descriptor;

  halfedge_descriptor hd = halfedge(first_fd, tm);
  std::size_t count = 0;

  while(count++ < 3)
  {
    if(face(opposite(hd, tm), tm) == second_fd)
      return hd;
    hd = next(hd, tm);
  }
  return boost::none;
}

// given two locations, return the common location with the largest dimension
#if 0
template<typename TriangleMesh>
boost::optional<typename boost::graph_traits<TriangleMesh>::descriptor_variant>
find_common_entity(typename Locate_types<TriangleMesh>::Face_location& first_loc,
                   typename Locate_types<TriangleMesh>::Face_location& second_loc,
                   const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor      face_descriptor;

  typedef typename boost::graph_traits<TriangleMesh>::descriptor_variant   descriptor_variant;

  if(first_loc.first != second_loc.first)
    return;

  descriptor_variant first_dv = get_descriptor_from_location(first_loc, tm);
  descriptor_variant second_dv = get_descriptor_from_location(second_loc, tm);

  // @ todo
}
#endif

} // namespace internal

// Point to barycentric coordinates
template <typename K>
CGAL::array<typename K::FT, 3>
barycentric_coordinates(const typename K::Point_2& p0, const typename K::Point_2& p1,
                        const typename K::Point_2& p2, const typename K::Point_2& query)
{
  typedef typename K::FT                        FT;
  typedef typename K::Vector_2                  Vector_2;

  K k;

  typename K::Compute_scalar_product_2 csp2 = k.compute_scalar_product_2_object();
  typename K::Construct_vector_2 cv2 = k.construct_vector_2_object();

  Vector_2 v0 = cv2(p0, p1);
  Vector_2 v1 = cv2(p0, p2);
  Vector_2 v2 = cv2(p0, query);

  FT d00 = csp2(v0, v0);
  FT d01 = csp2(v0, v1);
  FT d11 = csp2(v1, v1);
  FT d20 = csp2(v2, v0);
  FT d21 = csp2(v2, v1);

  FT denom = d00 * d11 - d01 * d01;

  FT v = (d11 * d20 - d01 * d21) / denom;
  FT w = (d00 * d21 - d01 * d20) / denom;

  return CGAL::make_array(FT(1.0) - v - w, v, w);
}

template <typename K>
CGAL::array<typename K::FT, 3>
barycentric_coordinates(const typename K::Point_3& p0, const typename K::Point_3& p1,
                        const typename K::Point_3& p2, const typename K::Point_3& query)
{
  typedef typename K::FT                        FT;
  typedef typename K::Vector_3                  Vector_3;

  K k;

  typename K::Compute_scalar_product_3 csp2 = k.compute_scalar_product_3_object();
  typename K::Construct_vector_3 cv2 = k.construct_vector_3_object();

  Vector_3 v0 = cv2(p0, p1);
  Vector_3 v1 = cv2(p0, p2);
  Vector_3 v2 = cv2(p0, query);

  FT d00 = csp2(v0, v0);
  FT d01 = csp2(v0, v1);
  FT d11 = csp2(v1, v1);
  FT d20 = csp2(v2, v0);
  FT d21 = csp2(v2, v1);

  FT denom = d00 * d11 - d01 * d01;

  FT v = (d11 * d20 - d01 * d21) / denom;
  FT w = (d00 * d21 - d01 * d20) / denom;

  return CGAL::make_array(FT(1.0) - v - w, v, w);
}

// return the number of 'next' one has to apply 'hd' to get source(hd, mesh) = vd,
// starting from hd = halfedge(fd, tm)
template <typename TriangleMesh>
int vertex_index_in_face(typename boost::graph_traits<TriangleMesh>::vertex_descriptor vd,
                         typename boost::graph_traits<TriangleMesh>::face_descriptor fd,
                         const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  halfedge_descriptor start = halfedge(fd, tm);
  halfedge_descriptor current = start;
  int count = 0;

  while (source(current, tm) != vd && count < 3)
  {
    current = next(current, tm);
    ++count;
  }

  // Just to make sure that 'vd' actually is a vertex of the face 'fd'
  CGAL_postcondition(count < 3);

  return count;
}

// return the number of 'next' one has to apply to get hd, starting from
// halfedge(face(hd, tm), tm)
template <typename TriangleMesh>
int halfedge_index_in_face(typename boost::graph_traits<TriangleMesh>::halfedge_descriptor he,
                           const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor     face_descriptor;

  CGAL_precondition(he != boost::graph_traits<TriangleMesh>::null_halfedge());
  CGAL_precondition(!is_border(he, tm));

  face_descriptor f = face(he, tm);
  halfedge_descriptor start = halfedge(f, tm);
  halfedge_descriptor current = start;
  int count = 0;

  while (current != he)
  {
    current = next(current, tm);
    ++count;
  }

  return count;
}

// Random locations
template <typename TriangleMesh>
typename internal::Locate_types<TriangleMesh>::Face_location
random_location_on_face(typename boost::graph_traits<TriangleMesh>::face_descriptor fd,
                        const TriangleMesh& /*tm*/,
                        CGAL::Random& rnd)
{
  typedef typename internal::Locate_types<TriangleMesh>::FT               FT;

  FT u = rnd.uniform_real(FT(0.0), FT(1.0));
  FT v = rnd.uniform_real(FT(0.0), FT(FT(1.0) - u));
  return std::make_pair(fd, boost::make_tuple(u, v, FT(FT(1.0) - u - v)));
}

template <typename TriangleMesh>
typename internal::Locate_types<TriangleMesh>::Face_location
random_location_on_halfedge(typename boost::graph_traits<TriangleMesh>::halfedge_descriptor hd,
                            const TriangleMesh& tm,
                            CGAL::Random& rnd)
{
  typedef typename internal::Locate_types<TriangleMesh>::FT               FT;

  FT t = rnd.uniform_real(FT(0.0), FT(1.0));
  return locate(hd, t, tm);
}

// vertex, halfedge, or face descriptor from a Face_location
template <typename TriangleMesh>
typename internal::Locate_types<TriangleMesh>::descriptor_variant
get_descriptor_from_location(const typename internal::Locate_types<TriangleMesh>::Face_location& loc,
                             const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor     face_descriptor;

  const face_descriptor fd = loc.first;

  // the first barycentric coordinate corresponds to source(halfedge(fd, tm), tm)
  halfedge_descriptor hd = prev(halfedge(fd, tm), tm);

  // check if the point is a vertex
  for(int i=0; i<3; ++i)
  {
    if(loc.second[i] == 1.) // coordinate at target(hd, tm)
      return target(hd, tm);
    hd = next(hd, tm);
  }
  CGAL_assertion(hd == prev(halfedge(fd, tm), tm));

  // check if the point is on an edge
  for(int i=0; i<3; ++i)
  {
    if(loc.second[i] == 0) // coordinate at target(hd, tm)
      return prev(hd, tm);
    hd = next(hd, tm);
  }

  return fd;
}

// Location to Point
template <typename TriangleMesh, typename NamedParameters>
typename internal::Locate_types<TriangleMesh, NamedParameters>::Point
location_to_point(const typename internal::Locate_types<TriangleMesh>::Face_location& loc,
                  const TriangleMesh& tm,
                  const NamedParameters& np)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::type VertexPointMap;
  typedef typename VertexPointMap::value_type                             Point;

  using boost::choose_param;
  using boost::get_param;

  VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                                      get_const_property_map(vertex_point, tm));

  halfedge_descriptor hd = halfedge(loc.first, tm);
  const Point p0 = get(vpmap, source(hd, tm));
  const Point p1 = get(vpmap, target(hd, tm));
  const Point p2 = get(vpmap, target(next(hd, tm), tm));

  return CGAL::barycenter(p0, loc.second[0], p1, loc.second[1], p2, loc.second[2]);
}

template <typename TriangleMesh>
typename property_map_value<TriangleMesh, CGAL::vertex_point_t>::type
location_to_point(const typename internal::Locate_types<TriangleMesh>::Face_location& loc,
                  const TriangleMesh& tm)
{
  return location_to_point(loc, tm, parameters::all_default());
}

/// \brief Given a `Face_location`, that is an ordered pair composed of a
///        `boost::graph_traits<TriangleMesh>::face_descriptor` and an array
///        of barycentric coordinates, and a second face adjacent to the first,
///        returns whether the location is on the vertex `vd` or not.
///
/// \details If `tm` is the input graph and given the pair (`f`, `bc`)
///          such that `bc` is `(w0, w1, w2)`, the correspondance with the weights in `bc`
///          and the vertices of the face `f` is the following:
///          - `w0 = source(halfedge(f, tm), tm)`
///          - `w1 = target(halfedge(f, tm), tm)`
///          - `w2 = target(next(halfedge(f, tm), tm), tm)`
///
/// \param loc the location
/// \param vd the vertex
/// \param tm the triangle mesh
///
template <typename TriangleMesh>
bool
is_on_vertex(const typename internal::Locate_types<TriangleMesh>::Face_location& loc,
             const typename boost::graph_traits<TriangleMesh>::vertex_descriptor vd,
             const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor     vertex_descriptor;
  typedef typename internal::Locate_types<TriangleMesh>::descriptor_variant descriptor_variant;

  descriptor_variant dv = get_descriptor_from_location(loc, tm);

  if(dv.which() == 0)
  {
    vertex_descriptor loc_vd = boost::get<vertex_descriptor>(dv);
    return (vd == loc_vd);
  }

  return false;
}

/// \brief Given a `Face_location`, that is an ordered pair composed of a
///        `boost::graph_traits<TriangleMesh>::face_descriptor` and an array
///        of barycentric coordinates, and a second face adjacent to the first,
///        returns whether the location is on the halfedge `hd` or not.
///
/// \details If `tm` is the input graph and given the pair (`f`, `bc`)
///          such that `bc` is `(w0, w1, w2)`, the correspondance with the weights in `bc`
///          and the vertices of the face `f` is the following:
///          - `w0 = source(halfedge(f, tm), tm)`
///          - `w1 = target(halfedge(f, tm), tm)`
///          - `w2 = target(next(halfedge(f, tm), tm), tm)`
///
/// \param loc the location
/// \param hd the halfedge
/// \param tm the triangle mesh
///
template <typename TriangleMesh>
bool
is_on_halfedge(const typename internal::Locate_types<TriangleMesh>::Face_location& loc,
               const typename boost::graph_traits<TriangleMesh>::halfedge_descriptor hd,
               const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor     vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor   halfedge_descriptor;
  typedef typename internal::Locate_types<TriangleMesh>::descriptor_variant descriptor_variant;

  descriptor_variant dv = get_descriptor_from_location(loc, tm);

  if(dv.which() == 0)
  {
    vertex_descriptor loc_vd = boost::get<vertex_descriptor>(dv);
    return (loc_vd == source(hd, tm) || loc_vd == target(hd, tm));
  }
  else if(dv.which() == 1)
  {
    halfedge_descriptor loc_hd = boost::get<halfedge_descriptor>(dv);
    return (loc_hd == hd);
  }

  return false;
}

/// \brief Given a `Face_location`, that is an ordered pair composed of a
///        `boost::graph_traits<TriangleMesh>::face_descriptor` and an array
///        of barycentric coordinates, and a second face adjacent to the first,
///        returns whether the location is in the face (border included) or not.
///
/// \details If `tm` is the input graph and given the pair (`f`, `bc`)
///          such that `bc` is `(w0, w1, w2)`, the correspondance with the weights in `bc`
///          and the vertices of the face `f` is the following:
///          - `w0 = source(halfedge(f, tm), tm)`
///          - `w1 = target(halfedge(f, tm), tm)`
///          - `w2 = target(next(halfedge(f, tm), tm), tm)`
///
/// \param loc the location
/// \param tm the triangle mesh
///
template <typename TriangleMesh>
bool
is_in_face(const typename internal::Locate_types<TriangleMesh>::Face_location& loc,
           const TriangleMesh& /*tm*/)
{
  typedef typename internal::Locate_types<TriangleMesh>::Face_location    Face_location;
  typedef typename Face_location::second_type                             Barycentric_coordinates;

  const Barycentric_coordinates& baryc = loc.second;

  for(int i=0; i<3; ++i)
    if(baryc[i] < 0. || baryc[i] > 1.)
      return false;

  return true;
}

/// \brief Given a `Face_location`, that is an ordered pair composed of a
///        `boost::graph_traits<TriangleMesh>::face_descriptor` and an array
///        of barycentric coordinates, and a second face adjacent to the first,
///        returns whether the location is on the border of the face or not.
///
/// \details If `tm` is the input graph and given the pair (`f`, `bc`)
///          such that `bc` is `(w0, w1, w2)`, the correspondance with the weights in `bc`
///          and the vertices of the face `f` is the following:
///          - `w0 = source(halfedge(f, tm), tm)`
///          - `w1 = target(halfedge(f, tm), tm)`
///          - `w2 = target(next(halfedge(f, tm), tm), tm)`
///
/// \param loc the location
/// \param tm the triangle mesh
///
template <typename TriangleMesh>
bool
is_on_face_border(const typename internal::Locate_types<TriangleMesh>::Face_location& loc,
                  const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor     face_descriptor;

  typedef typename internal::Locate_types<TriangleMesh>::Face_location    Face_location;
  typedef typename Face_location::second_type                             Barycentric_coordinates;

  const face_descriptor fd = loc.first;
  const Barycentric_coordinates& baryc = loc.second;

  // the first barycentric coordinate corresponds to source(halfedge(fd, tm), tm)
  halfedge_descriptor hd = prev(halfedge(fd, tm), tm);

  // check if the point is a vertex
  for(int i=0; i<3; ++i)
  {
    if(baryc[i] == 1.) // coordinate at target(hd, tm)
      return true;
    hd = next(hd, tm);
  }
  CGAL_assertion(hd == prev(halfedge(fd, tm), tm));

  // check if the point is on an edge
  for(int i=0; i<3; ++i)
  {
    if(baryc[i] == 0) // coordinate at target(hd, tm)
      return true;
    hd = next(hd, tm);
  }

  // point is strictly within the face, so it's not on the border
  return false;
}

/// \brief Given a `Face_location`, that is an ordered pair composed of a
///        `boost::graph_traits<TriangleMesh>::face_descriptor` and an array
///        of barycentric coordinates, and a second face adjacent to the first,
///        returns whether the location is on the border of the mesh or not.
///
/// \details If `tm` is the input graph and given the pair (`f`, `bc`)
///          such that `bc` is `(w0, w1, w2)`, the correspondance with the weights in `bc`
///          and the vertices of the face `f` is the following:
///          - `w0 = source(halfedge(f, tm), tm)`
///          - `w1 = target(halfedge(f, tm), tm)`
///          - `w2 = target(next(halfedge(f, tm), tm), tm)`
///
/// \param loc the location
/// \param tm the triangle mesh
///
template <typename TriangleMesh>
bool
is_on_mesh_border(const typename internal::Locate_types<TriangleMesh>::Face_location& loc,
                  const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor     face_descriptor;

  typedef typename internal::Locate_types<TriangleMesh>::Face_location    Face_location;
  typedef typename Face_location::second_type                             Barycentric_coordinates;

  const face_descriptor fd = loc.first;
  const Barycentric_coordinates& baryc = loc.second;

  // the first barycentric coordinate corresponds to source(halfedge(fd, tm), tm)
  halfedge_descriptor hd = prev(halfedge(fd, tm), tm);

  // check if the point is a vertex
  for(int i=0; i<3; ++i)
  {
    if(baryc[i] == 1.) // coordinate at target(hd, tm)
      return CGAL::is_border(target(hd, tm), tm);
    hd = next(hd, tm);
  }
  CGAL_assertion(hd == prev(halfedge(fd, tm), tm));

  // check if the point is on an edge
  for(int i=0; i<3; ++i)
  {
    if(baryc[i] == 0) // coordinate at target(hd, tm)
      return CGAL::is_border(prev(hd, tm), tm);
    hd = next(hd, tm);
  }

  // point is strictly within the face, so it's not on the border
  return false;
}

/// \name Surface Face Location Constructions
/// @{

/// \brief Returns the location of the given vertex as a `Face_location`, that is
///        an ordered pair specifying a face containing the location and the
///        barycentric coordinates of that location in that face.
///
/// \details If `tm` is the input graph and given the pair (`f`, `bc`)
///          such that `bc` is `(w0, w1, w2)`, the correspondance with the weights in `bc`
///          and the vertices of the face `f` is the following:
///          - `w0 = source(halfedge(f, tm), tm)`
///          - `w1 = target(halfedge(f, tm), tm)`
///          - `w2 = target(next(halfedge(f, tm), tm), tm)`
///
/// \param v a vertex of the TriangleMesh `tm`
/// \param tm the triangle mesh to which `v` belongs
///
template <typename TriangleMesh>
typename internal::Locate_types<TriangleMesh>::Face_location
locate(typename boost::graph_traits<TriangleMesh>::vertex_descriptor v,
       const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor     face_descriptor;

  typedef typename internal::Locate_types<TriangleMesh>::FT               FT;

  halfedge_descriptor he = next(halfedge(v, tm), tm);
  face_descriptor fd = face(he, tm);

  FT coords[3] = { FT(0.0), FT(0.0), FT(0.0) };
  std::size_t edge_local_index = halfedge_index_in_face(he, tm);
  coords[edge_local_index] = FT(1.0);

  return std::make_pair(fd, CGAL::make_array(coords[0], coords[1], coords[2]));
}

/// \brief Returns a location along the given edge as a `Face_location`, that is
///        an ordered pair specifying a face containing the location and the
///        barycentric coordinates of that location in that face.
///
/// \details If `tm` is the input graph and given the pair (`f`, `bc`)
///          such that `bc` is `(w0, w1, w2)`, the correspondance with the weights in `bc`
///          and the vertices of the face `f` is the following:
///          - `w0 = source(halfedge(f, tm), tm)`
///          - `w1 = target(halfedge(f, tm), tm)`
///          - `w2 = target(next(halfedge(f, tm), tm), tm)`
///
/// \param he a halfedge of the input face graph
/// \param t  the parametric distance of the desired point along `he`
/// \param tm the triangle mesh to which `he` belongs
///
template <typename TriangleMesh>
typename internal::Locate_types<TriangleMesh>::Face_location
locate(typename boost::graph_traits<TriangleMesh>::halfedge_descriptor he,
       typename internal::Locate_types<TriangleMesh>::FT t,
       const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor     face_descriptor;
  typedef typename internal::Locate_types<TriangleMesh>::FT               FT;

  face_descriptor fd = face(he, tm);
  std::size_t edge_local_index = halfedge_index_in_face(he, tm);

  const FT one_minus_t(FT(1.0) - t);
  FT coords[3] = { FT(0.0), FT(0.0), FT(0.0) };
  coords[edge_local_index] = one_minus_t;
  coords[(edge_local_index + 1) % 3] = t;

  return std::make_pair(fd, CGAL::make_array(coords[0], coords[1], coords[2]));
}

/// \brief Returns a location along in a face as a `Face_location`, that is
///        an ordered pair specifying a face containing the location and the
///        barycentric coordinates of that location in that face.
///
/// \details If `tm` is the input graph and given the pair (`f`, `bc`)
///          such that `bc` is `(w0, w1, w2)`, the correspondance with the weights in `bc`
///          and the vertices of the face `f` is the following:
///          - `w0 = source(halfedge(f, tm), tm)`
///          - `w1 = target(halfedge(f, tm), tm)`
///          - `w2 = target(next(halfedge(f, tm), tm), tm)`
///
/// \param he a halfedge of the input face graph
/// \param t  the parametric distance of the desired point along `he`
/// \param tm the triangle mesh to which `he` belongs
///
template <typename TriangleMesh, typename NamedParameters>
typename internal::Locate_types<TriangleMesh>::Face_location
locate(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
       const typename internal::Locate_types<TriangleMesh, NamedParameters>::Point& query,
       const TriangleMesh& tm,
       const NamedParameters& np)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor   vertex_descriptor;

  typedef typename internal::Locate_types<TriangleMesh>::Kernel           K;
  typedef typename K::FT                                                  FT;

  // VertexPointMap
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::type VertexPointMap;
  typedef typename boost::property_traits<VertexPointMap>::value_type     Point;

  using boost::choose_param;
  using boost::get_param;

  VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                                      get_const_property_map(vertex_point, tm));

  vertex_descriptor vd0 = source(halfedge(f, tm), tm);
  vertex_descriptor vd1 = target(halfedge(f, tm), tm);
  vertex_descriptor vd2 = target(next(halfedge(f, tm), tm), tm);

  const Point& p0 = get(vpmap, vd0);
  const Point& p1 = get(vpmap, vd1);
  const Point& p2 = get(vpmap, vd2);

  CGAL::cpp11::array<FT, 3> coords = barycentric_coordinates<K>(p0, p1, p2, query);

  if(coords[0] < 0. || coords[0] > 1. ||
     coords[1] < 0. || coords[1] > 1. ||
     coords[2] < 0. || coords[2] > 1.)
  {
    std::cerr << "Warning: point " << query << " is not in the face " << f << std::endl;
    std::cerr << "Coordinates: " << coords[0] << " " << coords[1] << " " << coords[2] << std::endl;

    // Try to to snap the coordinates, hoping the problem is just a -10e-17ish epsilon
    // pushing the coordinates over the edge
    // @tmp keep or not ?
    internal::snap_coordinates_to_border<TriangleMesh>(coords);
  }

  return std::make_pair(f, coords);
}

template <typename TriangleMesh>
typename internal::Locate_types<TriangleMesh>::Face_location
locate(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
       const typename property_map_value<TriangleMesh, CGAL::vertex_point_t>::type& query,
       const TriangleMesh& tm)
{
  return locate(f, query, tm, parameters::all_default());
}

/// \brief Given a `Face_location`, that is an ordered pair composed of a
///        `boost::graph_traits<TriangleMesh>::face_descriptor` and an array
///        of barycentric coordinates, and a second face adjacent to the first,
///        return the `Face_location` of the point in the second face.
///
/// \details If `tm` is the input graph and given the pair (`f`, `bc`)
///          such that `bc` is `(w0, w1, w2)`, the correspondance with the weights in `bc`
///          and the vertices of the face `f` is the following:
///          - `w0 = source(halfedge(f, tm), tm)`
///          - `w1 = target(halfedge(f, tm), tm)`
///          - `w2 = target(next(halfedge(f, tm), tm), tm)`
///
/// \param loc the first location
/// \param f the second face, adjacent to loc.first
/// \param tm the triangle mesh to which `he` belongs
///
/// \pre loc is on a subface common to `loc.first` and `f`
///
template <typename TriangleMesh>
typename internal::Locate_types<TriangleMesh>::Face_location
locate(const typename internal::Locate_types<TriangleMesh>::Face_location& loc,
       const typename boost::graph_traits<TriangleMesh>::face_descriptor f,
       const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor   vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor     face_descriptor;

  typedef boost::variant<vertex_descriptor,
                         halfedge_descriptor,
                         face_descriptor>                                descriptor_variant;

  typedef typename internal::Locate_types<TriangleMesh>::Face_location   Face_location;

  Face_location loc_in_f = std::make_pair(f, CGAL::make_array(0.,0.,0.));
  descriptor_variant dv = get_descriptor_from_location(loc, tm);

  // This function is meaningless if the point is not on the border of the face
  CGAL_precondition(dv.which() != 2);

  if(dv.which() == 0) // we're on a vertex
  {
    vertex_descriptor vd = boost::get<vertex_descriptor>(dv);
    int index_of_vd = vertex_index_in_face(vd, f, tm);
    loc_in_f.second[index_of_vd] = 1.;
    // note that the barycentric coordinates were initialized at 0,
    // so the second and third coordinates are already set up properly
  }
  else // dv.which() == 1, we're on a halfedge (and importantly, not on a vertex!)
  {
    halfedge_descriptor hd = boost::get<halfedge_descriptor>(dv);
    halfedge_descriptor opp_hd = opposite(hd, tm);
    CGAL_assertion(face(hd, tm) == loc.first);
    CGAL_assertion(face(opp_hd, tm) == f);

    int index_of_hd = halfedge_index_in_face(hd, tm);
    int index_of_opp_hd = halfedge_index_in_face(opp_hd, tm);

    // - Coordinates will be non-null at indices `index_of_hd`
    //   and `index_of_hd + 1` in loc.first.
    // - Coordinates will be non-null at indices `index_of_opp_hd`
    //   and `index_of_opp_hd + 1` in f.
    // - The halfedges `hd` and `opp_hd` have opposite directions.
    loc_in_f.second[index_of_opp_hd] = loc.second[(index_of_hd + 1)%3];
    loc_in_f.second[(index_of_opp_hd + 1)%3] = loc.second[index_of_hd];
    // note that the barycentric coordinates were initialized at 0,
    // so the third coordinate is already set up properly
  }

  // @todo clean up
//  std::cout << "Ini loc: " << loc.first << " b: " << loc.second[0] << " " << loc.second[1] << " " << loc.second[2] << std::endl;
//  std::cout << "Out loc: " << loc_in_f.first << " b: " << loc_in_f.second[0] << " " << loc_in_f.second[1] << " " << loc_in_f.second[2] << std::endl;

  CGAL_postcondition(loc_in_f.first == f);
  return loc_in_f;
}

// Finding a common face to a location and a point
// - the first location must be known
// - the second must be a point in a face incident to get_descriptor_from_location(known_location)
template <typename TriangleMesh>
bool
locate_in_common_face(const typename internal::Locate_types<TriangleMesh>::Point& query, //@todo namedparameters?
                      typename internal::Locate_types<TriangleMesh>::Face_location& known_location,
                      typename internal::Locate_types<TriangleMesh>::Face_location& query_location,
                      const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor      face_descriptor;

  typedef boost::variant<vertex_descriptor, halfedge_descriptor, face_descriptor> descriptor_variant;

  std::cout << "locate in common face:" << std::endl;
  std::cout << "query: " << query << std::endl;
  std::cout << "known location: " << known_location.first
            << " bc: " << known_location.second[0] << " "
                      << known_location.second[1] << " "
                      << known_location.second[2] << std::endl;

  descriptor_variant dv = get_descriptor_from_location(known_location, tm);

  bool is_query_location_in_face = false;

  if(dv.which() == 0)
  {
    vertex_descriptor vd = boost::get<vertex_descriptor>(dv);
    std::cout << "known_location on a vertex: " << vd << std::endl;
    halfedge_descriptor hd = halfedge(vd, tm);

    BOOST_FOREACH(face_descriptor fd, CGAL::faces_around_target(hd, tm))
    {
      // check if query can be found in that face
      query_location = locate(fd, query, tm);

      // @tmp keep or not ?
      internal::snap_location_to_border<TriangleMesh>(query_location, 1e-15);

      is_query_location_in_face = is_in_face(query_location, tm);

      std::cout << "trying face: " << fd << std::endl;
      std::cout << "location: " << query_location.first
                << " bc: " << query_location.second[0] << " "
                           << query_location.second[1] << " "
                           << query_location.second[2] << std::endl;

      if(is_query_location_in_face)
        break;
    }
  }
  else if(dv.which() == 1)
  {
    halfedge_descriptor hd = boost::get<halfedge_descriptor>(dv);
    std::cout << "known_location on a halfedge: " << hd << std::endl;
    face_descriptor fd = face(hd, tm);

    query_location = locate(fd, query, tm);
    internal::snap_location_to_border<TriangleMesh>(query_location); // @tmp keep or not ?
    is_query_location_in_face = is_in_face(query_location, tm);

    if(!is_query_location_in_face)
    {
      fd = face(opposite(hd, tm), tm);
      query_location = locate(fd, query, tm);
      is_query_location_in_face = is_in_face(query_location, tm);
    }
  }
  else // dv.which() == 2
  {
    CGAL_assertion(dv.which() == 2);
    face_descriptor fd = boost::get<face_descriptor>(dv);
    std::cout << "known_location on a face: " << fd << std::endl;

    query_location = locate(fd, query, tm);
    // @tmp keep or not ?
    internal::snap_location_to_border<TriangleMesh>(query_location);
    is_query_location_in_face = is_in_face(query_location, tm);
  }

  // if this is not the same face as for 'known_query', change 'known_location'
  if(is_query_location_in_face &&
     query_location.first != known_location.first)
  {
    known_location = locate(known_location, query_location.first, tm);
  }

  return is_query_location_in_face;
}

// Finding a common face to two locations
// - both locations must be known
template <typename TriangleMesh>
bool
locate_in_common_face(typename internal::Locate_types<TriangleMesh>::Face_location& first_location,
                      typename internal::Locate_types<TriangleMesh>::Face_location& second_location,
                      const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor      face_descriptor;

  std::cout << "locate in common face: " << std::endl;
  std::cout << "first location: " << first_location.first
            << " bc: " << first_location.second[0] << " "
                       << first_location.second[1] << " "
                       << first_location.second[2] << std::endl;

  std::cout << "second location: " << second_location.first
            << " bc: " << second_location.second[0] << " "
                       << second_location.second[1] << " "
                       << second_location.second[2] << std::endl;

  // Check that we actually have something to do
  if(first_location.first == second_location.first)
    return true;

  bool is_first_location_on_border = is_on_face_border(first_location, tm);
  bool is_second_location_on_border = is_on_face_border(second_location, tm);

  // We have already checked that they have different faces, if neither are on
  // a border, then it's hopeless
  if(!is_first_location_on_border && !is_second_location_on_border)
    return false;

  // Find a common face in the sets of incident faces of each location
  std::set<face_descriptor> first_incident_faces;
  std::set<face_descriptor> second_incident_faces;

  internal::incident_faces(first_location, first_incident_faces, tm);
  internal::incident_faces(second_location, second_incident_faces, tm);

  typename std::set<face_descriptor>::const_iterator fit = first_incident_faces.begin();
  typename std::set<face_descriptor>::const_iterator fend = first_incident_faces.end();
  typename std::set<face_descriptor>::const_iterator sit = second_incident_faces.begin();
  typename std::set<face_descriptor>::const_iterator send = second_incident_faces.end();

  while(fit!=fend && sit!=send)
  {
    if(*fit == *sit)
      break;
    else if(*fit < *sit)
      ++fit;
    else
      ++sit;
  }

  if(fit == fend || sit == send) // no common face...
    return false;

  CGAL_assertion(*fit == *sit);
  face_descriptor common_fd = *fit;

  if(first_location.first != common_fd)
    first_location = locate(first_location, common_fd, tm);

  if(second_location.first != common_fd)
    second_location = locate(second_location, common_fd, tm);

  return true;
}

/// @}

/// \name Nearest Face Location Queries
/// @{

/// \brief Creates an `AABB_tree` suitable for use with `locate`.
///
/// \tparam AABBTraits A model of `AABBTraits` used to define a \cgal `AABB_tree`.
///
/// \param outTree Output parameter to store the computed `AABB_tree`
///
template <class AABBTraits, typename TriangleMesh, typename NamedParameters>
void build_aabb_tree(const TriangleMesh& tm,
                     AABB_tree<AABBTraits>& outTree,
                     const NamedParameters& np)
{
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::type VertexPointMap;

  using boost::choose_param;
  using boost::get_param;

  VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                                      get_const_property_map(vertex_point, tm));

  typename boost::graph_traits<TriangleMesh>::face_iterator facesStart, facesEnd;
  boost::tie(facesStart, facesEnd) = faces(tm);
  outTree.rebuild(facesStart, facesEnd, tm, vpmap);
  outTree.build();
}

template <class AABBTraits, typename TriangleMesh>
void build_aabb_tree(const TriangleMesh& tm,
                     AABB_tree<AABBTraits>& outTree,
                     const typename AABBTraits::Primitive* /*dummy for sfinae*/ = 0)
{
  return build_aabb_tree(tm, outTree, parameters::all_default());
}

/// \brief Returns the face location nearest to the given point.
///
/// \tparam AABBTraits A model of `AABBTraits` used to define a \cgal `AABB_tree`.
///
/// \param p the point to locate on the input face graph
/// \param tree An `AABB_tree` containing the triangular faces of the input surface mesh to perform the point location with
///
template <class AABBTraits, typename TriangleMesh, typename NamedParameters>
typename internal::Locate_types<TriangleMesh>::Face_location
locate(const typename AABBTraits::Point_3& p,
       const AABB_tree<AABBTraits>& tree,
       const TriangleMesh& tm,
       const NamedParameters& np)
{
  typename AABB_tree<AABBTraits>::Point_and_primitive_id result =
    tree.closest_point_and_primitive(p);

  return locate(result.second, result.first, tm, np);
}

template <class AABBTraits, typename TriangleMesh>
typename internal::Locate_types<TriangleMesh>::Face_location
locate(const typename property_map_value<TriangleMesh, CGAL::vertex_point_t>::type& p,
       const AABB_tree<AABBTraits>& tree,
       const TriangleMesh& tm,
       const typename AABBTraits::Primitive* /*dummy for sfinae*/ = 0)
{
  return locate(p, tree, tm, parameters::all_default());
}

/// \brief Returns the nearest face location to the given point.
///  Note that this will (re-)build an `AABB_tree` on each call. If you need
///  to  call this function more than once, use `build_aabb_tree()` to cache a
///  copy of the `AABB_tree`, and use the overloads of this function
///  that accept a reference to an `AABB_tree` as input.
///
/// \tparam TriangleMesh A model of `FaceGraph`.
/// \tparam AABBTraits A model of `AABBTraits` used to define a \cgal `AABB_tree`.
///
/// \param p the point to locate on the input face graph
/// \param tm the input face graph
///
template <typename AABBTraits, typename TriangleMesh, typename NamedParameters>
typename internal::Locate_types<TriangleMesh>::Face_location
locate(const typename boost::property_traits<
         typename GetVertexPointMap<TriangleMesh, NamedParameters>::type>::value_type& p,
       const TriangleMesh& tm,
       const NamedParameters& np,
       const typename AABBTraits::Primitive* /* dummy for sfinae */ = 0)
{
  AABB_tree<AABBTraits> tree;
  build_aabb_tree(tm, tree, np);
  return locate(p, tree, tm, np);
}

template <typename TriangleMesh>
typename internal::Locate_types<TriangleMesh>::Face_location
locate(const typename property_map_value<TriangleMesh, CGAL::vertex_point_t>::type& p,
       const TriangleMesh& tm)
{
  // @fixme doesn't work for 2D mesh yet
  AABB_tree<AABB_face_graph_triangle_primitive<TriangleMesh> > tree;
  build_aabb_tree(tm, tree);
  return locate(p, tm, tree, parameters::all_default());
}

#if 0
/// \brief Returns the face location along `ray` nearest to its source point.
///
/// \tparam AABBTraits A model of `AABBTraits` used to define a \cgal `AABB_tree`.
///
/// \param ray Ray to intersect with the input face graph
/// \param tree A `AABB_tree` containing the triangular faces of the input surface mesh to perform the point location with
///
template <class AABBTraits>
Face_location locate(const Ray_3& ray,
                     const AABB_tree<AABBTraits>& tree,
                     const Triangle_mesh& tm,
                     Vertex_point_map vertexPointMap,
                     const Traits& traits = Traits())
{
  typedef AABB_tree<AABBTraits> AABB_face_graph_tree;
  typename Traits::Construct_barycentric_coordinates_in_triangle_3 cbcit3(traits.construct_barycentric_coordinates_in_triangle_3_object());
  typename Traits::Construct_barycentric_coordinates cbc(traits.construct_barycentric_coordinates_object());
  typename Traits::Compute_squared_distance_3 csd3(traits.compute_squared_distance_3_object());
  typedef typename AABB_face_graph_tree::template Intersection_and_primitive_id<Ray_3>::Type Intersection_type;
  typedef boost::optional<Intersection_type> Ray_intersection;

  std::vector<Ray_intersection> intersections;

  tree.all_intersections(ray, std::back_inserter(intersections));

  bool foundOne = false;
  FT nearestDistance = 0;
  Point_3 nearestPoint = CGAL::ORIGIN;
  face_descriptor nearestFace;

  for (std::size_t i = 0; i < intersections.size(); ++i)
  {
    if (intersections[i])
    {
      Point_3* intersectionPoint = boost::get<Point_3>(&(intersections[i]->first));

      if (intersectionPoint)
      {
        FT distance = csd3(*intersectionPoint, ray.source());

        if (!foundOne || distance < nearestDistance)
        {
          foundOne = true;
          nearestPoint = *intersectionPoint;
          nearestDistance = distance;
          nearestFace = intersections[i]->second;
        }
      }
    }
  }

  if (foundOne)
  {
    Barycentric_coordinates b = cbcit3(triangle_from_face(nearestFace, tm, vertexPointMap), nearestPoint);
    return std::make_pair(nearestFace, b);
  }
  else
  {
    return std::make_pair(boost::graph_traits<TriangleMesh>::null_face(),
                          CGAL::make_array(FT(0.0), FT(0.0), FT(0.0)));
  }
}

///
/// \brief Returns the face location along `ray` nearest to its source point.
///        Note that this will (re-)build an `AABB_tree` on each call. If you need
///        to call this function more than once, use `build_aabb_tree()` to cache a
///        copy of the `AABB_tree`, and use the overloads of this function
///        that accept a reference to an `AABB_tree` as input.
///
/// \tparam AABBTraits A model of `AABBTraits` used to define an `AABB_tree`.
///
/// \param ray Ray to intersect with the input face graph
///
template <class AABBTraits>
Face_location locate(const Ray_3& ray,
                     const Triangle_mesh& tm,
                     Vertex_point_map vertexPointMap,
                     const Traits& traits = Traits())
{
  AABB_tree<AABBTraits> tree;
  build_aabb_tree(tm, tree);
  return locate(ray, tree, tm, vertexPointMap, traits);
}
#endif

/// @}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_LOCATE_H
