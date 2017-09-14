// Copyright (c) 2014, 2017 GeometryFactory (France).
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

#include <CGAL/array.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/variant.hpp>

#include <utility>

namespace CGAL{
namespace Polygon_mesh_processing {

namespace internal {

// return the number of 'next' one has to apply to get hd, starting from
// halfedge(face(hd, mesh), mesh)
template <class TriangleMesh>
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

template<typename Face_location, typename PolygonMesh>
boost::variant<typename boost::graph_traits<PolygonMesh>::vertex_descriptor,
               typename boost::graph_traits<PolygonMesh>::halfedge_descriptor,
               typename boost::graph_traits<PolygonMesh>::face_descriptor>
get_descriptor_from_location(const Face_location& loc, const PolygonMesh& mesh)
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor     face_descriptor;

  typedef typename Face_location::second_type                            Barycentric_coordinates;

  const face_descriptor fd = loc.first;
  const Barycentric_coordinates& baryc = loc.second;

  // the first barycentric coordinate corresponds to source(halfedge(fd, mesh), mesh)
  halfedge_descriptor hd = prev(halfedge(fd, mesh), mesh);

  // check if the point is a vertex
  for(int i=0; i<3; ++i)
  {
    if(baryc[i] == 1.) // coordinate at target(hd, mesh)
      return target(hd, mesh);
    hd = next(hd, mesh);
  }
  CGAL_assertion(hd == prev(halfedge(fd, mesh), mesh));

  // check if the point is on an edge
  for(int i=0; i<3; ++i)
  {
    if(baryc[i] == 0) // coordinate at target(hd, mesh)
      return prev(hd, mesh);
    hd = next(hd, mesh);
  }

  return fd;
}

#if 0
template<typename Point, typename Face_location, typename PolygonMesh>
Point loc_to_point(const Face_location& loc, const PolygonMesh& mesh)
{
  CGAL_assertion(false);
  return Point(); // @todo
}
#endif

} // namespace internal

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
template<typename TriangleMesh>
std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
          CGAL::cpp11::array<
            typename CGAL::Kernel_traits<
              typename property_map_value<TriangleMesh,
                                          CGAL::vertex_point_t>::type>::Kernel::FT, 3> >
face_location(typename boost::graph_traits<TriangleMesh>::vertex_descriptor v,
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
template<typename TriangleMesh>
std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
          CGAL::cpp11::array<
            typename CGAL::Kernel_traits<
              typename property_map_value<TriangleMesh,
                                          CGAL::vertex_point_t>::type>::Kernel::FT, 3> >
face_location(typename boost::graph_traits<TriangleMesh>::halfedge_descriptor he,
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
template<typename TriangleMesh, typename NamedParameters>
std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
          CGAL::cpp11::array<
            typename CGAL::Kernel_traits<
              typename property_map_value<TriangleMesh,
                                          CGAL::vertex_point_t>::type>::Kernel::FT, 3> >
face_location(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
              typename property_map_value<TriangleMesh, CGAL::vertex_point_t>::type& query,
              const TriangleMesh& tm,
              const NamedParameters& np)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor   vertex_descriptor;

  typedef typename CGAL::Kernel_traits<
            typename property_map_value<TriangleMesh,
              CGAL::vertex_point_t>::type>::Kernel                       K;
  typedef typename K::FT                                                 FT;

  // VertexPointMap
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type       VPMap;
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::type::value_type Point;

  VPMap vpmap = choose_param(get_param(np, vertex_point),
                             get_const_property_map(vertex_point, tm));

  vertex_descriptor vd0 = source(halfedge(f, tm), tm);
  vertex_descriptor vd1 = target(halfedge(f, tm), tm);
  vertex_descriptor vd2 = target(next(halfedge(f, tm), tm), tm);

  const Point& v0 = get(vpmap, vd0);
  const Point& v1 = get(vpmap, vd1);
  const Point& v2 = get(vpmap, vd2);

  FT inverted_area = 1./CGAL::area(v0, v1, v2);
  FT b0 = CGAL::area(v1, v2, query) * inverted_area;
  FT b1 = CGAL::area(v2, v0, query) * inverted_area;
  FT b2 = 1. - b0 - b1;

  return std::make_pair(f, CGAL::make_array(b0, b1, b2));
}

template<typename TriangleMesh>
std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
          CGAL::cpp11::array<
            typename CGAL::Kernel_traits<
              typename property_map_value<TriangleMesh,
                                          CGAL::vertex_point_t>::type>::Kernel::FT, 3> >
face_location(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
              typename property_map_value<TriangleMesh, CGAL::vertex_point_t>::type& query,
              const TriangleMesh& tm)
{
  return face_location(f, query, tm, parameters::all_default());
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
template<typename TriangleMesh>
std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
          CGAL::cpp11::array<
            typename CGAL::Kernel_traits<
              typename property_map_value<TriangleMesh,
                                          CGAL::vertex_point_t>::type>::Kernel::FT, 3> >
face_location(const std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
                              CGAL::cpp11::array<
                                typename CGAL::Kernel_traits<
                                  typename property_map_value<TriangleMesh,
                                                              CGAL::vertex_point_t>::type>::Kernel::FT, 3> >& loc,
              const typename boost::graph_traits<TriangleMesh>::face_descriptor f,
              const TriangleMesh& mesh)
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

  descriptor_variant dv = internal::get_descriptor_from_location(loc, mesh);

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
    halfedge_descriptor opp_hd = opposite(hd, mesh);
    CGAL_assertion(face(hd, mesh) == loc.first);
    CGAL_assertion(face(opp_hd, mesh) == f);

    int index_of_hd = internal::edge_index_in_face(hd, mesh);
    int index_of_opp_hd = internal::edge_index_in_face(opp_hd, mesh);

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

// Given a location, returns whether a point is on the border of the mesh or not.
template<typename TriangleMesh>
bool is_on_mesh_border(const std::pair<typename boost::graph_traits<TriangleMesh>::face_descriptor,
                                       CGAL::cpp11::array<
                                         typename CGAL::Kernel_traits<
                                           typename property_map_value<TriangleMesh,
                                                                       CGAL::vertex_point_t>::type>::Kernel::FT, 3> >& loc,
                       const TriangleMesh& mesh)
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

  // the first barycentric coordinate corresponds to source(halfedge(fd, mesh), mesh)
  halfedge_descriptor hd = prev(halfedge(fd, mesh), mesh);

  // check if the point is a vertex
  for(int i=0; i<3; ++i)
  {
    if(baryc[i] == 1.) // coordinate at target(hd, mesh)
      return CGAL::is_border(target(hd, mesh), mesh);
    hd = next(hd, mesh);
  }
  CGAL_assertion(hd == prev(halfedge(fd, mesh), mesh));

  // check if the point is on an edge
  for(int i=0; i<3; ++i)
  {
    if(baryc[i] == 0) // coordinate at target(hd, mesh)
      return CGAL::is_border(prev(hd, mesh), mesh);
    hd = next(hd, mesh);
  }

  // point is strictly within the face, so it's not on the border
  return false;
}

#if 0 // not adapted yet
/// \name Nearest Face Location Queries
/// @{

/// \brief Returns the nearest face location to the given point.
///  Note that this will (re-)build an `AABB_tree` on each call. If you need
///  to  call this function more than once, use `build_aabb_tree()` to cache a
///  copy of the `AABB_tree`, and use the overloads of this function
///  that accept a reference to an `AABB_tree` as input.
///
/// \tparam AABBTraits A model of `AABBTraits` used to define a \cgal `AABB_tree`.
///
/// \param p Point to locate on the input face graph
///
template <class AABBTraits>
Face_location
locate(const Point_3& location,
       const Triangle_mesh& tm,
       Vertex_point_map vertexPointMap)
{
  AABB_tree<AABBTraits> tree;
  build_aabb_tree(tm, tree);
  return locate(location, tree, tm, vertexPointMap, traits);
}


/// \brief Returns the face location nearest to the given point.
///
///
/// \tparam AABBTraits A model of `AABBTraits` used to define a \cgal `AABB_tree`.
///
/// \param p Point to locate on the input face graph
/// \param tree A `AABB_tree` containing the triangular faces of the input surface mesh to perform the point location with
///
template <class AABBTraits>
Face_location
locate(const Point_3& location,
       const AABB_tree<AABBTraits>& tree,
       const Triangle_mesh& tm,
       Vertex_point_map vertexPointMap)
{
  typename Traits::Construct_barycentric_coordinates_in_triangle_3 cbcit3(traits.construct_barycentric_coordinates_in_triangle_3_object());
  typename AABB_tree<AABBTraits>::Point_and_primitive_id result = tree.closest_point_and_primitive(location);

  face_descriptor f = result.second;
  Barycentric_coordinates b = cbcit3(triangle_from_face(f, tm, vertexPointMap), result.first);
  return std::make_pair(f, b);
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

/// @}

/// \brief Creates an `AABB_tree` suitable for use with `locate`.
///
/// \tparam AABBTraits A model of `AABBTraits` used to define a \cgal `AABB_tree`.
///
/// \param outTree Output parameter to store the computed `AABB_tree`
///
template <class AABBTraits>
void build_aabb_tree(const Triangle_mesh& tm, AABB_tree<AABBTraits>& outTree)
{
  face_iterator facesStart, facesEnd;
  boost::tie(facesStart, facesEnd) = faces(tm);
  outTree.rebuild(facesStart, facesEnd, tm);
  outTree.build();
}

#endif

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_LOCATE_H
