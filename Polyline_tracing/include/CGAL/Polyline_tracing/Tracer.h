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

#include <CGAL/Polygon_mesh_processing/locate.h>

#include <boost/tuple/tuple.hpp>
#include <boost/variant.hpp>

namespace CGAL {

namespace Polyline_tracing {

template<typename MotorcycleGraphTraits, typename Visitor>
class Tracer
{
public:
  typedef MotorcycleGraphTraits                               Geom_traits;
  typedef typename Geom_traits::Triangle_mesh                 Triangle_mesh;

  typedef typename Geom_traits::FT                            FT;

  typedef Dictionary<Geom_traits>                             Dictionary;
  typedef typename Dictionary::DEC_it                         DEC_it;
  typedef Dictionary_entry<Geom_traits>                       Dictionary_entry;

  typedef typename Geom_traits::Face_location                 Face_location;

  typedef typename Geom_traits::vertex_descriptor             vertex_descriptor;
  typedef typename Geom_traits::halfedge_descriptor           halfedge_descriptor;
  typedef typename Geom_traits::face_descriptor               face_descriptor;

  typedef boost::variant<vertex_descriptor,
                         halfedge_descriptor,
                         face_descriptor>                 descriptor_variant;

  // templated by 'Motorcycle' to avoid a circular dependency
  template<typename Motorcycle>
  boost::tuple<bool, DEC_it, DEC_it, FT, bool> trace(const Motorcycle& mc,
                                                     Dictionary& points,
                                                     const Triangle_mesh& mesh);
};

// -----------------------------------------------------------------------------

template<typename MotorcycleGraphTraits, typename Visitor>
template<typename Motorcycle>
boost::tuple<bool, // successfuly computed a next path or not
             typename Tracer<MotorcycleGraphTraits, Visitor>::DEC_it, // next_source
             typename Tracer<MotorcycleGraphTraits, Visitor>::DEC_it, // next_destination
             typename Tracer<MotorcycleGraphTraits, Visitor>::FT, // time_at_next_destination
             bool> // whether the destination is final or not
Tracer<MotorcycleGraphTraits, Visitor>::
trace(const Motorcycle& mc, Dictionary& points, const Triangle_mesh& mesh)
{
#ifdef CGAL_MOTORCYCLE_GRAPH_VERBOSE
  std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*" << std::endl;
  std::cout << "Computing the next path for motorcycle: " << mc.id() << std::endl;
  std::cout << "Current position: " << mc.position()->point() << std::endl
            << "Location: " << mc.current_location().first << " b: "
            << mc.current_location().second [0] << " "
            << mc.current_location().second [1] << " "
            << mc.current_location().second [2] << std::endl;
#endif

  // just to get rid of a degenerate case
  CGAL_precondition(mc.direction()); // direction must be known
  if(*(mc.direction()) == CGAL::NULL_VECTOR)
  {
    std::cerr << "Warning: the motorcycle direction is null and "
              << "the next destination is thus the current position" << std::endl;

    return boost::make_tuple(true, mc.position(), mc.position(),
                             mc.current_time(), true /*final destination*/);
  }

  const Face_location& loc = mc.current_location();
  descriptor_variant dv =
    CGAL::Polygon_mesh_processing::internal::get_descriptor_from_location(loc, mesh);
  return boost::apply_visitor(Visitor(&mc, points, mesh), dv);
}

} // namespace Polyline_tracing

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_TRACER_H
