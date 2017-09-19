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
#include <CGAL/Kernel_traits.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/variant.hpp>

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

// @todo check that all compile (sfinae probably to use a few times)

namespace CGAL{
namespace Polygon_mesh_processing {

namespace internal {

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

// return the number of 'next' one has to apply to get hd, starting from
// halfedge(face(hd, tm), tm)
template <typename TriangleMesh>
int edge_index_in_face(typename boost::graph_traits<TriangleMesh>::halfedge_descriptor he,
                       const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor     face_descriptor;

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

template <typename TriangleMesh>
boost::variant<typename boost::graph_traits<TriangleMesh>::vertex_descriptor,
               typename boost::graph_traits<TriangleMesh>::halfedge_descriptor,
               typename boost::graph_traits<TriangleMesh>::face_descriptor>
get_descriptor_from_location(const std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
                                             CGAL::cpp11::array<
                                               typename CGAL::Kernel_traits<
                                                 typename property_map_value<
                                                   TriangleMesh, CGAL::vertex_point_t>::type>::Kernel::FT, 3> >& loc,
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

template <typename TriangleMesh, typename NamedParameters>
typename boost::property_traits<
  typename GetVertexPointMap<TriangleMesh, NamedParameters>::type>::value_type
loc_to_point(const std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
                             CGAL::cpp11::array<
                               typename CGAL::Kernel_traits<
                                 typename boost::property_traits<
                                   typename GetVertexPointMap<
                                     TriangleMesh, NamedParameters>::type>::value_type>::Kernel::FT, 3> >& loc,
             const TriangleMesh& tm,
             const NamedParameters& np)
{
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::type VertexPointMap;
  typedef typename VertexPointMap::value_type                             Point;

  using boost::choose_param;
  using boost::get_param;

  VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                                      get_const_property_map(vertex_point, tm));

  CGAL_assertion(false);
  return Point(); // @todo
}

template <typename TriangleMesh>
typename property_map_value<TriangleMesh, CGAL::vertex_point_t>::type
loc_to_point(const std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
                             CGAL::cpp11::array<
                               typename CGAL::Kernel_traits<
                                 typename property_map_value<TriangleMesh,
                                                             CGAL::vertex_point_t>::type>::Kernel::FT, 3> >& loc,
             const TriangleMesh& tm)
{
  return internal::loc_to_point(loc, tm, parameters::all_default());
}

} // namespace internal

// Given a location, returns whether a point is on the border of the mesh or not.
template <typename TriangleMesh>
bool
is_on_mesh_border(const std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
                                  CGAL::cpp11::array<
                                    typename CGAL::Kernel_traits<
                                      typename property_map_value<TriangleMesh,
                                                                  CGAL::vertex_point_t>::type>::Kernel::FT, 3> >& loc,
                  const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor     face_descriptor;

  typedef typename CGAL::Kernel_traits<
            typename property_map_value<TriangleMesh,
              CGAL::vertex_point_t>::type>::Kernel                        K;
  typedef typename K::FT                                                  FT;

  typedef std::pair<face_descriptor, CGAL::cpp11::array<FT, 3> >          Face_location;
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
std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
          CGAL::cpp11::array<
            typename CGAL::Kernel_traits<
              typename property_map_value<TriangleMesh,
                                          CGAL::vertex_point_t>::type>::Kernel::FT, 3> >
locate(typename boost::graph_traits<TriangleMesh>::vertex_descriptor v,
       const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor     face_descriptor;

  typedef typename CGAL::Kernel_traits<
            typename property_map_value<TriangleMesh,
              CGAL::vertex_point_t>::type>::Kernel                       K;
  typedef typename K::FT                                                 FT;

  halfedge_descriptor he = next(halfedge(v, tm), tm);
  face_descriptor fd = face(he, tm);

  FT coords[3] = { FT(0.0), FT(0.0), FT(0.0) };
  std::size_t edge_local_index = internal::edge_index_in_face(he, tm);
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
std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
          CGAL::cpp11::array<
            typename CGAL::Kernel_traits<
              typename property_map_value<TriangleMesh,
                                          CGAL::vertex_point_t>::type>::Kernel::FT, 3> >
locate(typename boost::graph_traits<TriangleMesh>::halfedge_descriptor he,
         typename CGAL::Kernel_traits<
           typename property_map_value<TriangleMesh,
                                       CGAL::vertex_point_t>::type>::Kernel::FT t,
       const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor     face_descriptor;

  typedef typename CGAL::Kernel_traits<
            typename property_map_value<TriangleMesh,
              CGAL::vertex_point_t>::type>::Kernel                       K;
  typedef typename K::FT                                                 FT;

  face_descriptor fd = face(he, tm);
  std::size_t edge_local_index = internal::edge_index_in_face(he, tm);

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
std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
          CGAL::cpp11::array<
            typename CGAL::Kernel_traits<
              typename boost::property_traits<
                typename GetVertexPointMap<
                  TriangleMesh, NamedParameters>::type>::value_type>::Kernel::FT, 3> >
locate(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
       const typename boost::property_traits<typename GetVertexPointMap<
         TriangleMesh, NamedParameters>::type>::value_type& query,
       const TriangleMesh& tm,
       const NamedParameters& np)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor   vertex_descriptor;

  typedef typename CGAL::Kernel_traits<
            typename property_map_value<TriangleMesh,
              CGAL::vertex_point_t>::type>::Kernel                        K;
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

  CGAL::cpp11::array<FT, 3> coords = internal::barycentric_coordinates<K>(p0, p1, p2, query);

  if(coords[0] < 0. || coords[0] > 1. ||
     coords[1] < 0. || coords[1] > 1. ||
     coords[2] < 0. || coords[2] > 1.)
  {
    std::cerr << "Warning: point " << query << " is not in the face " << f << std::endl;
    std::cerr << "Coordinates: " << coords[0] << " " << coords[1] << " " << coords[2] << std::endl;
  }

  return std::make_pair(f, coords);
}

template <typename TriangleMesh>
std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
          CGAL::cpp11::array<
            typename CGAL::Kernel_traits<
              typename property_map_value<TriangleMesh,
                                          CGAL::vertex_point_t>::type>::Kernel::FT, 3> >
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
template <typename TriangleMesh>
std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
          CGAL::cpp11::array<
            typename CGAL::Kernel_traits<
              typename property_map_value<TriangleMesh,
                                          CGAL::vertex_point_t>::type>::Kernel::FT, 3> >
locate(const std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
                         CGAL::cpp11::array<
                           typename CGAL::Kernel_traits<
                             typename property_map_value<TriangleMesh,
                                                         CGAL::vertex_point_t>::type>::Kernel::FT, 3> >& loc,
       const typename boost::graph_traits<TriangleMesh>::face_descriptor f,
       const TriangleMesh& tm)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor   vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor     face_descriptor;

  typedef boost::variant<vertex_descriptor,
                         halfedge_descriptor,
                         face_descriptor>                                descriptor_variant;

  typedef typename CGAL::Kernel_traits<
            typename property_map_value<TriangleMesh,
              CGAL::vertex_point_t>::type>::Kernel                       K;
  typedef typename K::FT                                                 FT;

  typedef std::pair<face_descriptor, CGAL::cpp11::array<FT, 3> >         Face_location;

  Face_location loc_in_f = std::make_pair(f, CGAL::make_array(0.,0.,0.));

  descriptor_variant dv = internal::get_descriptor_from_location(loc, tm);

  // This function is meaningless if the point is on the border of the face
  CGAL_precondition(dv.which() != 2);

  if(dv.which() == 0) // we're on a vertex
  {
    // vertex_descriptor vd = boost::get<vertex_descriptor>(dv);
    CGAL_assertion(false); // @todo
  }
  else // if(dv.which() == 1) // we're on a halfedge (and importantly, not on a vertex!)
  {
    halfedge_descriptor hd = boost::get<halfedge_descriptor>(dv);
    halfedge_descriptor opp_hd = opposite(hd, tm);
    CGAL_assertion(face(hd, tm) == loc.first);
    CGAL_assertion(face(opp_hd, tm) == f);

    int index_of_hd = internal::edge_index_in_face(hd, tm);
    int index_of_opp_hd = internal::edge_index_in_face(opp_hd, tm);

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

  std::cout << "Ini loc: " << loc.first << " b: " << loc.second[0] << " " << loc.second[1] << " " << loc.second[2] << std::endl;
  std::cout << "Out loc: " << loc_in_f.first << " b: " << loc_in_f.second[0] << " " << loc_in_f.second[1] << " " << loc_in_f.second[2] << std::endl;

  CGAL_postcondition(loc_in_f.first == f);
  return loc_in_f;
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
void build_aabb_tree(const TriangleMesh& tm, AABB_tree<AABBTraits>& outTree)
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
std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
          CGAL::cpp11::array<
            typename CGAL::Kernel_traits<
              typename boost::property_traits<
                typename GetVertexPointMap<
                  TriangleMesh, NamedParameters>::type>::value_type>::Kernel::FT, 3> >
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
std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
          CGAL::cpp11::array<
            typename CGAL::Kernel_traits<
              typename property_map_value<TriangleMesh,
                                          CGAL::vertex_point_t>::type>::Kernel::FT, 3> >
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
std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
          CGAL::cpp11::array<
            typename CGAL::Kernel_traits<
              typename boost::property_traits<
                typename GetVertexPointMap<
                  TriangleMesh, NamedParameters>::type>::value_type>::Kernel::FT, 3> >
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
std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
          CGAL::cpp11::array<
            typename CGAL::Kernel_traits<
              typename property_map_value<TriangleMesh,
                                          CGAL::vertex_point_t>::type>::Kernel::FT, 3> >
locate(const typename property_map_value<TriangleMesh, CGAL::vertex_point_t>::type& p,
       const TriangleMesh& tm)
{
  AABB_tree<AABB_face_graph_triangle_primitive<TriangleMesh> > tree;
  build_aabb_tree(tm, tree);
  return locate(p, tm, tree, parameters::all_default());
}

#if 0
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
Face_location locate(const Ray_3& ray, const Triangle_mesh& tm, Vertex_point_map vertexPointMap, const Traits& traits = Traits())
{
  AABB_tree<AABBTraits> tree;
  build_aabb_tree(tm, tree);
  return locate(ray, tree, tm, vertexPointMap, traits);
}

/// \brief Returns the face location along `ray` nearest to its source point.
///
/// \tparam AABBTraits A model of `AABBTraits` used to define a \cgal `AABB_tree`.
///
/// \param ray Ray to intersect with the input face graph
/// \param tree A `AABB_tree` containing the triangular faces of the input surface mesh to perform the point location with
///
template <class AABBTraits>
Face_location locate(const Ray_3& ray, const AABB_tree<AABBTraits>& tree, const Triangle_mesh& tm, Vertex_point_map vertexPointMap, const Traits& traits = Traits())
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
#endif

/// @}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_LOCATE_H
