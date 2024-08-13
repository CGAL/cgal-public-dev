// Copyright (c) 2018  Carnegie Mellon University (USA), GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Christina Vaz, Keenan Crane, Andreas Fabri

#ifndef CGAL_INTRINSIC_MOLLIFICATION_3_H
#define CGAL_INTRINSIC_MOLLIFICATION_3_H

#include <CGAL/license/Heat_method_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Dynamic_property_map.h>
#include <CGAL/Heat_method_3/internal/V2V.h>
#include <CGAL/Iterator_range.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/double.h>
#include <CGAL/number_utils.h>
#include <CGAL/property_map.h>
#include <CGAL/squared_distance_3.h>

#include <algorithm>
#include <boost/iterator/transform_iterator.hpp>

#include <cmath>

namespace CGAL
{
namespace Heat_method_3
{
/**
 * A tag class used to specify that the heat method applies the intrinsic mollification to the input
 * mesh.
 */
namespace internal
{
struct Mollification
{
};

template <typename T>
static constexpr bool is_mollification_scheme_v = std::is_base_of<Mollification, T>::value;
}  // namespace internal

struct Intrinsic_mollification_constant : internal::Mollification
{
};

namespace internal
{
template <typename T>
static constexpr bool always_false = false;

template <typename Scheme>
struct mollification_scheme
{
  template <
    typename TriangleMesh,
    typename VertexPointMap,
    typename HalfEdgeLengthMap,
    typename NamedParameters = CGAL::parameters::Default_named_parameters,
    typename Traits = typename Kernel_traits<typename boost::property_traits<
      typename boost::property_map<TriangleMesh, vertex_point_t>::const_type>::value_type>::Kernel>
  static void apply(const TriangleMesh &input_tm,
                    const VertexPointMap &vpm,
                    HalfEdgeLengthMap &halfedge_length_map,
                    const NamedParameters &np = CGAL::parameters::default_values())
  {
    // mode doesn't need mollification
    // static_assert(internal::always_false<Scheme>, "mollification scheme not implemented");
  }
};

/*!
\ingroup PkgHeatMethodRef
@brief mollify sliver triangles.

 @tparam TriangleMesh a triangulated surface mesh, model of `FaceListGraph` and `HalfedgeListGraph`
 @tparam Traits a model of `HeatMethodTraits_3`

  @tparam PolygonMesh must be model of `MutableFaceGraph`
  @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  @param pmesh polygon mesh which has the hole
  @param border_halfedge a border halfedge incident to the hole
  @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones
listed below

  \cgalNamedParamsBegin

    \cgalParamNBegin{face_output_iterator}
      \cgalParamDescription{iterator over patch faces}
      \cgalParamType{a model of `OutputIterator`
    holding `boost::graph_traits<PolygonMesh>::%face_descriptor` for patch faces}
      \cgalParamDefault{`Emptyset_iterator`}
    \cgalParamNEnd

  \cgalNamedParamsEnd

  @return void
  */
template <>
struct mollification_scheme<Intrinsic_mollification_constant>
{
  template <
    typename TriangleMesh,
    typename VertexPointMap,
    typename HalfEdgeLengthMap,
    typename NamedParameters = CGAL::parameters::Default_named_parameters,
    typename Traits = typename Kernel_traits<typename boost::property_traits<
      typename boost::property_map<TriangleMesh, vertex_point_t>::const_type>::value_type>::Kernel>

  static void apply(const TriangleMesh &input_tm,
                    const VertexPointMap &vpm,
                    HalfEdgeLengthMap &halfedge_length_map,
                    const NamedParameters &np = CGAL::parameters::default_values())
  {
    typedef boost::graph_traits<TriangleMesh> graph_traits;
    typedef typename graph_traits::halfedge_descriptor halfedge_descriptor;
    /// Geometric typedefs
    typedef typename Traits::Point_3 Point_3;
    typedef typename Traits::FT FT;
    typedef typename Traits::Vector_3 Vector_3;

    typedef int Index;

    using parameters::choose_parameter;
    using parameters::get_parameter;

    FT delta = choose_parameter(
      get_parameter(np, internal_np::delta),
      FT{1e-4
         * (*std::min_element(halfedges(input_tm).begin(),
                              halfedges(input_tm).end(),
                              [&](halfedge_descriptor hd1, halfedge_descriptor hd2) {
                                return get(halfedge_length_map, hd1)
                                       < get(halfedge_length_map, hd2);
                              }))});

    // compute smallest length epsilon we can add to
    // all edges to ensure that the strict triangle
    // inequality holds with a tolerance of delta
    FT epsilon = 0;
    for (halfedge_descriptor hd : halfedges(input_tm))
      {
        halfedge_descriptor hd2 = next(hd, input_tm);
        halfedge_descriptor hd3 = next(hd2, input_tm);
        FT ineq = get(halfedge_length_map, hd2) + get(halfedge_length_map, hd3)
                  - get(halfedge_length_map, hd);
        epsilon = (std::max)(epsilon, (std::max)(0., delta - ineq));
      }
    std::cout << "im epsilon = " << epsilon << "\n";
    // update edge lengths
    for (halfedge_descriptor hd : halfedges(input_tm))
      {
        get(halfedge_length_map, hd) += epsilon;
      }
  }
};
}  // namespace internal
}  // namespace Heat_method_3
}  // namespace CGAL

#include <CGAL/enable_warnings.h>
#endif  // CGAL_INTRINSIC_MOLLIFICATION_3_H
