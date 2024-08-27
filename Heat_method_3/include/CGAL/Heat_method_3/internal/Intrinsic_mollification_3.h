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
#include <CGAL/Weights/utils.h>
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

namespace CGAL {
namespace Heat_method_3 {
/*!
 * \ingroup PkgHeatMethodRef
 *
 * \brief A concept for a mollification scheme.
 *
 * This concept defines the requirements for a mollification scheme used in mesh processing. A mollification
 * scheme is a technique applied to triangle meshes to improve their quality, particularly to address
 * issues such as degenerate triangles. Types that model this concept must provide implementations of 
 * two functions: `apply` and `cotangent`.
 */
struct Mollification_scheme_identity {
  /*!
   * \brief Applies the mollification scheme to the triangle mesh.
   *
   * This function modifies or initializes the provided triangle mesh, vertex point map, and halfedge
   * length map according to the mollification scheme. The specific behavior depends on the implementation
   * of the mollification scheme. Typically, this function adjusts vertex positions, edge lengths,
   * or other properties to enhance the quality of the mesh.
   *
   * \param tm The input triangle mesh to be processed.
   * \param vpm The property map that provides the points of the vertices in the mesh.
   * \param np Optional named parameters for customizing the function's behavior. Defaults to
   * CGAL::parameters::default_values().
   *
   * \tparam TriangleMesh The type of the triangle mesh.
   * \tparam VertexPointMap The type of the vertex point property map.
   * \tparam NamedParameters The type of named parameters (defaults to CGAL::parameters::Default_named_parameters).
   * \tparam Traits The traits class used for geometric operations (defaults to a traits class inferred from the
   * vertex point map).
   *
   * \return A property map that provides the length of the halfedges in the mesh.
   */
  template <typename TriangleMesh,
            typename VertexPointMap,
            typename NamedParameters = CGAL::parameters::Default_named_parameters,
            typename Traits = typename Kernel_traits<typename boost::property_traits<typename boost::property_map<TriangleMesh, vertex_point_t>::const_type>::value_type>::Kernel>
  static auto apply(const TriangleMesh &tm,
                    const VertexPointMap &vpm,
                    const NamedParameters &np = CGAL::parameters::default_values())
  {
    typedef typename Traits::FT FT;
    typedef CGAL::dynamic_halfedge_property_t<FT> Halfedge_length_tag;
    typedef typename boost::property_map<TriangleMesh, Halfedge_length_tag>::const_type
        HalfedgeLengthMap;
    return HalfedgeLengthMap{};
  }

  /*!
   * \brief Computes the cotangent weight for a given triangle.
   *
   * This function calculates the cotangent weight for the angle opposite to the edge between `v1` and `v2` 
   * in the triangle defined by vertices `v1`, `v2`, and `v3`. The cotangent weight is used in various 
   * geometric algorithms, including those involving the discrete Laplacian.
   *
   * \param tm The input triangle mesh.
   * \param vpm The property map for vertex points.
   * \param he_length_map The property map for the lengths of halfedges in the mesh.
   * \param v1 The first vertex of the triangle.
   * \param v2 The second vertex of the triangle.
   * \param v3 The third vertex of the triangle.
   * \param traits The traits class providing geometric operations.
   *
   * \tparam TriangleMesh The type of the triangle mesh.
   * \tparam VertexPointMap The type of the vertex point property map.
   * \tparam HalfedgeLengthMap The type of the halfedge length property map.
   * \tparam Traits The traits class used for geometric operations.
   *
   * \return The cotangent weight of the angle opposite to the edge between `v1` and `v2`.
   */
  template <typename TriangleMesh,
            typename VertexPointMap,
            typename HalfedgeLengthMap,
            typename Traits>
  static auto cotangent(const TriangleMesh &tm,
                        const VertexPointMap &vpm,
                        const HalfedgeLengthMap &he_length_map,
                        typename boost::graph_traits<TriangleMesh>::vertex_descriptor v1,
                        typename boost::graph_traits<TriangleMesh>::vertex_descriptor v2,
                        typename boost::graph_traits<TriangleMesh>::vertex_descriptor v3,
                        Traits traits)
  {
    typedef typename Traits::FT FT;
    typedef typename boost::property_traits<VertexPointMap>::reference VertexPointMap_reference;

    VertexPointMap_reference p_i = get(vpm, v1);
    VertexPointMap_reference p_j = get(vpm, v2);
    VertexPointMap_reference p_k = get(vpm, v3);
    return CGAL::Weights::cotangent(p_i, p_j, p_k, traits);
  }
};

/*!
 * \brief Applies the constant mollification scheme to the triangle mesh.
 *
 * This function modifies the provided triangle mesh by adding a constant epsilon to each edge
 * length to ensure the strict triangle inequality holds with a tolerance of `delta`. The epsilon
 * value is computed based on the smallest edge length in the mesh and the specified `delta`.
 *
 */
struct mollification_scheme_constant {
  template <typename TriangleMesh,
            typename VertexPointMap,
            typename NamedParameters = CGAL::parameters::Default_named_parameters,
            typename Traits = typename Kernel_traits<typename boost::property_traits<typename boost::
              property_map<TriangleMesh, vertex_point_t>::const_type>::value_type>::Kernel>
  static auto apply(const TriangleMesh &tm,
      const VertexPointMap &vpm,
      const NamedParameters &np = CGAL::parameters::default_values())
  {
    typedef boost::graph_traits<TriangleMesh> graph_traits;
    typedef typename graph_traits::halfedge_descriptor halfedge_descriptor;
    typedef typename graph_traits::vertex_descriptor vertex_descriptor;
    typedef typename graph_traits::edge_descriptor edge_descriptor;
    /// Geometric typedefs
    typedef typename Traits::Point_3 Point_3;
    typedef typename Traits::FT FT;
    typedef typename Traits::Vector_3 Vector_3;
    typename Traits::Compute_squared_distance_3 squared_distance =
        Traits().compute_squared_distance_3_object();
    typedef typename boost::property_traits<VertexPointMap>::reference VertexPointMap_reference;

    typedef int Index;

    using parameters::choose_parameter;
    using parameters::get_parameter;

    typedef CGAL::dynamic_halfedge_property_t<FT> Halfedge_length_tag;
    typedef typename boost::property_map<TriangleMesh, Halfedge_length_tag>::const_type
        HalfedgeLengthMap;

    HalfedgeLengthMap he_length_map(get(Halfedge_length_tag(), tm));

    FT min_length = std::numeric_limits<double>::max();
    for (auto e : edges(tm)) {
      halfedge_descriptor hd1 = halfedge(e, tm);
      halfedge_descriptor hd2 = opposite(hd1, tm);
      vertex_descriptor v1 = target(e, tm);
      vertex_descriptor v2 = source(e, tm);
      VertexPointMap_reference p1 = get(vpm, v1);
      VertexPointMap_reference p2 = get(vpm, v2);
      FT e_length = CGAL::approximate_sqrt(squared_distance(p1, p2));
      put(he_length_map, hd1, e_length);
      put(he_length_map, hd2, e_length);
      min_length = e_length > FT(0) ? CGAL::min(min_length, e_length) : min_length;
    }

    FT delta = choose_parameter(get_parameter(np, internal_np::delta), FT{1e-4 * min_length});
    std::cout << "delta = " << delta << "\n";
    // compute smallest length epsilon we can add to
    // all edges to ensure that the strict triangle
    // inequality holds with a tolerance of delta
    FT epsilon = 0;
    for (halfedge_descriptor hd : halfedges(tm)) {
      halfedge_descriptor hd2 = next(hd, tm);
      halfedge_descriptor hd3 = next(hd2, tm);
      FT ineq = get(he_length_map, hd2) + get(he_length_map, hd3) - get(he_length_map, hd);
      epsilon = (std::max)(epsilon, (std::max)(0., delta - ineq));
    }
    std::cout << "im epsilon = " << epsilon << "\n";

    // update edge lengths
    for (halfedge_descriptor hd : halfedges(tm)) {
      put(he_length_map, hd, epsilon + get(he_length_map, hd));
    }

    return he_length_map;
  }

  template <typename TriangleMesh,
            typename VertexPointMap,
            typename HalfedgeLengthMap,
            typename Traits>
  static auto cotangent(const TriangleMesh &tm,
      const VertexPointMap &vpm,
      const HalfedgeLengthMap &he_length_map,
      typename boost::graph_traits<TriangleMesh>::vertex_descriptor v1,
      typename boost::graph_traits<TriangleMesh>::vertex_descriptor v2,
      typename boost::graph_traits<TriangleMesh>::vertex_descriptor v3,
      Traits traits)
  {
    typedef typename Traits::FT FT;

    // Retrieve edge lengths from the halfedge length map
    const FT a = get(he_length_map, halfedge(v1, v2, tm).first);
    const FT b = get(he_length_map, halfedge(v2, v3, tm).first);
    const FT c = get(he_length_map, halfedge(v3, v1, tm).first);

    // Compute the semi-perimeter of the triangle
    const FT S = (a + b + c) / 2;

    // Compute the area of the triangle using Heron's formula
    const FT area = CGAL::approximate_sqrt(S * (S - a) * (S - b) * (S - c));

    // If the area is zero, return an undefined value
    if (is_zero(area)) {
      return FT(0);  // Undefined cotangent
    }

    assert(std::isfinite(area));
    // Compute the cotangent of the angle opposite to edge length c
    return (CGAL::square(a) + CGAL::square(b) - CGAL::square(c)) / (4 * area);
  }
};

struct mollification_scheme_local_one_by_one {
  template <typename TriangleMesh,
            typename VertexPointMap,
            typename NamedParameters = CGAL::parameters::Default_named_parameters,
            typename Traits = typename Kernel_traits<typename boost::property_traits<typename boost::
              property_map<TriangleMesh, vertex_point_t>::const_type>::value_type>::Kernel>
  static auto apply(const TriangleMesh &tm,
      const VertexPointMap &vpm,
      const NamedParameters &np = CGAL::parameters::default_values())
  {
    typedef boost::graph_traits<TriangleMesh> graph_traits;
    typedef typename graph_traits::halfedge_descriptor halfedge_descriptor;
    typedef typename graph_traits::vertex_descriptor vertex_descriptor;
    typedef typename graph_traits::edge_descriptor edge_descriptor;
    /// Geometric typedefs
    typedef typename Traits::Point_3 Point_3;
    typedef typename Traits::FT FT;
    typedef typename Traits::Vector_3 Vector_3;
    typename Traits::Compute_squared_distance_3 squared_distance =
        Traits().compute_squared_distance_3_object();
    typedef typename boost::property_traits<VertexPointMap>::reference VertexPointMap_reference;

    typedef int Index;

    using parameters::choose_parameter;
    using parameters::get_parameter;

    typedef CGAL::dynamic_halfedge_property_t<FT> Halfedge_length_tag;
    typedef typename boost::property_map<TriangleMesh, Halfedge_length_tag>::const_type
        HalfedgeLengthMap;

    HalfedgeLengthMap he_length_map(get(Halfedge_length_tag(), tm));

    FT min_length = std::numeric_limits<double>::max();
    for (auto e : edges(tm)) {
      halfedge_descriptor hd1 = halfedge(e, tm);
      halfedge_descriptor hd2 = opposite(hd1, tm);
      vertex_descriptor v1 = target(e, tm);
      vertex_descriptor v2 = source(e, tm);
      VertexPointMap_reference p1 = get(vpm, v1);
      VertexPointMap_reference p2 = get(vpm, v2);
      FT e_length = CGAL::approximate_sqrt(squared_distance(p1, p2));
      put(he_length_map, hd1, e_length);
      put(he_length_map, hd2, e_length);
      min_length = e_length > FT(0) ? CGAL::min(min_length, e_length) : min_length;
    }

    // TODO: add threshold parameter instead of constant 1e-4
    FT delta = choose_parameter(get_parameter(np, internal_np::delta), FT{1e-4 * min_length});
    std::cout << "delta = " << delta << "\n";
    // compute smallest length epsilon we can add to
    // all edges to ensure that the strict triangle
    // inequality holds with a tolerance of delta
    FT epsilon = 0;
    for (halfedge_descriptor hd : halfedges(tm)) {
      halfedge_descriptor hd2 = next(hd, tm);
      halfedge_descriptor hd3 = next(hd2, tm);
      std::array<halfedge_descriptor, 3> he = {hd, hd2, hd3};
      std::array<double, 3> L = {
          get(he_length_map, hd), get(he_length_map, hd2), get(he_length_map, hd3)};
      // Create indices vector and sort based on values
      std::vector<int> indices(3);
      std::iota(indices.begin(), indices.end(), 0);  // Fill indices with 0, 1, ..., size-1

      // Sort indices based on values in `row`
      std::sort(indices.begin(), indices.end(), [&](int i1, int i2) { return L[i1] < L[i2]; });
      // Reorder a, b, c so that c <= b <= a
      auto a = L[indices[2]];
      auto b = L[indices[1]];
      auto c = L[indices[0]];

      if (a - b - c > 1e-4)
        continue;




      // Mollify
      c = std::max(c, delta + a - b);
      c = std::max(c, delta + b - a);
      b = std::max(b, c);
      a = std::max(a, b);

      // Reorder back to original order
      put(he_length_map, he[indices[0]], c);
      put(he_length_map, he[indices[1]], b);
      put(he_length_map, he[indices[2]], a);
    }

    return he_length_map;
  }

  template <typename TriangleMesh,
            typename VertexPointMap,
            typename HalfedgeLengthMap,
            typename Traits>
  static auto cotangent(const TriangleMesh &tm,
      const VertexPointMap &vpm,
      const HalfedgeLengthMap &he_length_map,
      typename boost::graph_traits<TriangleMesh>::vertex_descriptor v1,
      typename boost::graph_traits<TriangleMesh>::vertex_descriptor v2,
      typename boost::graph_traits<TriangleMesh>::vertex_descriptor v3,
      Traits traits)
  {
    typedef typename Traits::FT FT;

    // Retrieve edge lengths from the halfedge length map
    const FT a = get(he_length_map, halfedge(v1, v2, tm).first);
    const FT b = get(he_length_map, halfedge(v2, v3, tm).first);
    const FT c = get(he_length_map, halfedge(v3, v1, tm).first);

    // Compute the semi-perimeter of the triangle
    const FT S = (a + b + c) / 2;

    // Compute the area of the triangle using Heron's formula
    const FT area = CGAL::approximate_sqrt(S * (S - a) * (S - b) * (S - c));

    // If the area is zero, return an undefined value
    if (is_zero(area)) {
      return FT(0);  // Undefined cotangent
    }

    assert(std::isfinite(area));
    // Compute the cotangent of the angle opposite to edge length c
    return (CGAL::square(a) + CGAL::square(b) - CGAL::square(c)) / (4 * area);
  }
};
}  // namespace Heat_method_3
}  // namespace CGAL

#include <CGAL/enable_warnings.h>
#endif  // CGAL_INTRINSIC_MOLLIFICATION_3_H
