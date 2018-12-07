// Copyright (c) 2017 GeometryFactory
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
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_POLYLINE_TRACING_MOTORCYCLE_GRAPH_TRAITS_2_H
#define CGAL_POLYLINE_TRACING_MOTORCYCLE_GRAPH_TRAITS_2_H

#include <CGAL/Polyline_tracing/internal/robust_collinear.h>

#include <CGAL/array.h>
#include <CGAL/Bbox_2.h>

#include <boost/graph/graph_traits.hpp>

#include <utility>

namespace CGAL {

namespace Polyline_tracing {

/*!
\ingroup PkgPolylineTracingTraitsClasses

\brief A model of the concept `MotorcycleGraphTraits`
as required by the `Motorcycle_graph` class.

\tparam K A \cgal Kernel
\tparam InputMesh A model of `FaceListGraph`
\tparam OutputGraph A model of `MutableHalfedgeGraph`

\cgalModels `MotorcycleGraphTraits`
*/
template <typename K,
          typename InputMesh,
          typename OutputGraph = InputMesh>
class Motorcycle_graph_traits_2
  : public K
{
  typedef K                                                   Base;

public:
  /// Kernel type
  typedef K                                                   Kernel;

  /// Triangle mesh type
  typedef InputMesh                                           Triangle_mesh;

  /// Motorcycle Graph output type
  typedef OutputGraph                                         Face_graph;

  typedef typename Kernel::FT                                 FT;
  typedef typename Kernel::Point_2                            Point_d;
  typedef typename Kernel::Segment_2                          Segment_d;
  typedef typename Kernel::Vector_2                           Vector_d;
  typedef typename Kernel::Ray_2                              Ray_d;
  typedef typename Kernel::Triangle_2                         Triangle_d;

  typedef CGAL::Bbox_2                                        Bbox_d;

private:
  // Graph traits
  typedef typename boost::graph_traits<Triangle_mesh>::face_descriptor      face_descriptor;

public:
  /// Barycentric coordinates type
  ///
  // Points are not described through a Point_d, but through an ordered pair
  // specifying a location on the surface of the `Triangle_mesh`.
  //
  // If `tm` is the input graph and given the pair (`f`, `bc`)
  // such that `bc` is `(w0, w1, w2)`, the correspondance with the weights in `bc`
  // and the vertices of the face `f` is the following:
  // - `w0 = source(halfedge(f, tm), tm)`
  // - `w1 = target(halfedge(f, tm), tm)`
  // - `w2 = target(next(halfedge(f, tm), tm), tm)`
  typedef typename CGAL::cpp11::array<FT, 3>                  Barycentric_coordinates;
  typedef std::pair<face_descriptor, Barycentric_coordinates> Face_location;

  // 2D robust predicates
  typedef internal::Robust_collinear_are_strictly_ordered_along_line_2<K>
                                                              Collinear_are_strictly_ordered_along_line_2;

  // Predicates types
  typedef typename K::Are_ordered_along_line_2                Are_ordered_along_line_d;
  typedef typename K::Collinear_2                             Collinear_d;
  typedef typename K::Collinear_are_ordered_along_line_2      Collinear_are_ordered_along_line_d;
  typedef typename K::Do_intersect_2                          Do_intersect_d;

  // Constructions types
  typedef typename K::Angle_2                                 Angle_d;
  typedef typename K::Compute_squared_distance_2              Compute_squared_distance_d;
  typedef typename K::Construct_midpoint_2                    Construct_midpoint_d;
  typedef typename K::Construct_segment_2                     Construct_segment_d;
  typedef typename K::Intersect_2                             Intersect_d;

  // Constructor
  explicit Motorcycle_graph_traits_2(const K& k = K()) : Base(k) { }

  static int dimension() { return dim; }

  // 2D robust predicate
  Collinear_are_strictly_ordered_along_line_2
  collinear_are_strictly_ordered_along_line_2_object() const
  { return Collinear_are_strictly_ordered_along_line_2(static_cast<const Base&>(*this)); }

  // Predicates
  Are_ordered_along_line_d
  are_ordered_along_line_d_object() const
  { return Base::are_ordered_along_line_2_object(); }

  Collinear_d
  collinear_d_object() const
  { return Base::collinear_2_object(); }

  Collinear_are_ordered_along_line_d
  collinear_are_ordered_along_line_d_object() const
  { return Base::collinear_are_ordered_along_line_2_object(); }

  Do_intersect_d
  do_intersect_d_object() const
  { return Base::do_intersect_2_object(); }

  // Constructions
  Angle_d
  angle_d_object() const
  { return Base::angle_2_object(); }

  Compute_squared_distance_d
  compute_squared_distance_d_object() const
  { return Base::compute_squared_distance_2_object(); }

  Construct_midpoint_d
  construct_midpoint_d_object() const
  { return Base::construct_midpoint_2_object(); }

  Construct_segment_d
  construct_segment_d_object() const
  { return Base::construct_segment_2_object(); }

  Intersect_d
  intersect_d_object() const
  { return Base::intersect_2_object(); }

public:
  static const int dim = 2;
};

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_MOTORCYCLE_GRAPH_TRAITS_2_H
