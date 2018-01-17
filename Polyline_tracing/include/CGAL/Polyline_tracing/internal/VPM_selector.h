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

#ifndef CGAL_POLYLINE_TRACING_VPM_SELECTOR_H
#define CGAL_POLYLINE_TRACING_VPM_SELECTOR_H

#include <CGAL/Kernel_traits.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

namespace CGAL {

namespace internal {

template<typename PolygonMesh>
struct P2_or_P3_to_P3
{
  typedef typename CGAL::Kernel_traits<
            typename property_map_value<PolygonMesh,
              CGAL::vertex_point_t>::type>::Kernel                          K;
  typedef typename K::Point_2                                               Point_2;
  typedef typename K::Point_3                                               Point_3;

  P2_or_P3_to_P3() { }

  Point_3 operator()(const Point_2& p) const { return Point_3(p.x(), p.y(), 0.); }
  const Point_3& operator()(const Point_3& p) const { return p; }
};

} // namespace internal

template<typename PolygonMesh/*, typename NamedParameters*/>
struct P2_to_P3_VPM
{
private:
  typedef P2_to_P3_VPM<PolygonMesh/*, NamedParameters*/>    Self;

public:
//  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type      VertexPointMap;
  typedef typename property_map_selector<PolygonMesh, boost::vertex_point_t>::const_type VertexPointMap;

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor      vertex_descriptor;

  typedef internal::P2_or_P3_to_P3<PolygonMesh>                             P2_or_P3_to_P3;

  typedef typename CGAL::Kernel_traits<
            typename property_map_value<PolygonMesh,
              CGAL::vertex_point_t>::type>::Kernel                          K;
  typedef typename K::Point_2                                               Point_2;
  typedef typename K::Point_3                                               Point_3;

  // classical typedefs
  typedef vertex_descriptor                   key_type;
  typedef Point_3                             value_type;
  typedef value_type                          reference;
  typedef boost::readable_property_map_tag    category;

  // Constructor
  P2_to_P3_VPM() { }
  P2_to_P3_VPM(const PolygonMesh& mesh)
    : conv(), vpmap(get_const_property_map(vertex_point, mesh))
  { }

//  P2_to_P3_VPM(const PolygonMesh& mesh, const NamedParameters& np)
//    :
//      vpmap(choose_param(get_param(np, internal_np::vertex_point),
//                         get_property_map(vertex_point, mesh)))
//  { }

  // Access
  const P2_or_P3_to_P3& converter() const { return conv; }
  const VertexPointMap& vpm() const { return vpmap; }

  // get function for property map
  inline friend
  reference
  get(const Self& pmap, key_type v)
  {
    return pmap.converter()(get(pmap.vpm(), v));
  }

private:
  P2_or_P3_to_P3 conv;
  VertexPointMap vpmap;
};

template<typename PolygonMesh, typename NamedParameters>
struct PolygonMesh_dimension_vpmap_type_selector
{
  typedef typename Polygon_mesh_processing::GetVertexPointMap<PolygonMesh,
                                                              NamedParameters>::type  VPMNative;
  typedef P2_to_P3_VPM<PolygonMesh/*, NamedParameters*/>                              VPM2D;

  // if the dimension is 2, the vpmap transforms from Point_2 to Point_3
  // by assigning 0 to axis Oz

  // otherwise, returns the default pmap

};

} // namespace CGAL

#endif // CGAL_POLYLINE_TRACING_VPM_SELECTOR_H
