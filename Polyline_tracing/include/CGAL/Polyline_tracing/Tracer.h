// Copyright (c) 2017 GeometryFactory (France).
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
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_POLYLINE_TRACING_TRACER_H
#define CGAL_POLYLINE_TRACING_TRACER_H

#include <CGAL/Polyline_tracing/Dictionary.h>

#include <boost/tuple/tuple.hpp>
#include <boost/variant.hpp>

namespace CGAL {

namespace Polyline_tracing {

template<typename K, typename PolygonMesh, typename Visitor>
class Tracer
{
public:
  typedef typename K::FT                                  FT;

  typedef Dictionary<K, PolygonMesh>                      Dictionary;
  typedef typename Dictionary::DEC_it                     DEC_it;
  typedef Dictionary_entry<K, PolygonMesh>                Dictionary_entry;
  typedef typename Dictionary_entry::Face_location        Face_location;

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor  halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor      face_descriptor;

  typedef boost::variant<vertex_descriptor,
                         halfedge_descriptor,
                         face_descriptor>                 descriptor_variant;

  template<typename Motorcycle>
  boost::tuple<DEC_it, DEC_it, FT> trace(const Motorcycle& mc, Dictionary& points, const PolygonMesh& mesh);
};

// -----------------------------------------------------------------------------

template<typename K, typename PolygonMesh, typename Visitor>
template<typename Motorcycle>
boost::tuple<typename Tracer<K, PolygonMesh, Visitor>::DEC_it,
             typename Tracer<K, PolygonMesh, Visitor>::DEC_it,
             typename Tracer<K, PolygonMesh, Visitor>::FT>
Tracer<K, PolygonMesh, Visitor>::
trace(const Motorcycle& mc, Dictionary& points, const PolygonMesh& mesh)
{
  const Face_location& loc = mc.current_location();
  descriptor_variant dv = internal::get_descriptor_from_location(loc, mesh);
  return boost::apply_visitor(Visitor(&mc, points, mesh), dv);
}

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_TRACER_H
