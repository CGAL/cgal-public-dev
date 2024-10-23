// Copyright (c) 2024 GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s): Islam Omar, Hossam Saeed, Mael Rouxel-Labbé

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
#include <limits>

#include <CGAL/MP_Float.h>
#include <CGAL/QP_functions.h>
#include <CGAL/QP_models.h>

#include <osqp++.h>
#include <Eigen/SparseCore>

namespace CGAL {
namespace Heat_method_3 {
// Default delta for the mollification scheme
constexpr auto Kdelta = 1e-2;

/*!
 * \ingroup PkgHeatMethodRef
 *
 * \brief A concept for a mollification scheme.
 *
 * This concept defines the requirements for a mollification scheme used in mesh processing.
 * A mollification scheme is a technique applied to triangle meshes to improve their quality,
 * particularly in addressing issues such as degenerate triangles. Types that model this concept
 * must provide implementations of two functions: `apply` and `cotangent`.
 */
struct Mollification_scheme_identity {
  /*!
   * \brief Applies the mollification scheme to the triangle mesh.
   *
   * This function modifies or initializes the provided triangle mesh, vertex point map, and
   * halfedge length map according to the specified mollification scheme. The specific behavior 
   * depends on the implementation of the mollification scheme. Typically, this function adjusts 
   * vertex positions, edge lengths, or other properties to enhance the mesh's quality.
   *
   * \param tm The input triangle mesh to be processed.
   * \param vpm The property map that provides the points of the vertices in the mesh.
   * \param np An optional named parameter for specifying the mollification behavior, with details below:
   * \cgalNamedParamsBegin
   * \cgalParamNBegin{delta}
   *   \cgalParamDescription{A parameter used to set the relative length change (relative to 
   *   average edge lengths) in various mollification schemes.}
   *   \cgalParamType{float}
   *   \cgalParamDefault{`1e-2`}
   * \cgalParamNEnd
   * \cgalNamedParamsEnd
   *
   * \tparam TriangleMesh The type of the triangle mesh.
   * \tparam VertexPointMap The type of the vertex point property map.
   * \tparam NamedParameters The type of named parameters (defaults to 
   *   CGAL::parameters::Default_named_parameters).
   * \tparam Traits The traits class used for geometric operations (defaults to a traits class 
   *   inferred from the vertex point map).
   *
   * \return A property map that provides the lengths of the halfedges in the mesh.
   */
  template <typename TriangleMesh,
            typename VertexPointMap,
            typename NamedParameters = CGAL::parameters::Default_named_parameters,
            typename Traits = typename Kernel_traits<typename boost::property_traits<typename boost::
                property_map<TriangleMesh, vertex_point_t>::const_type>::value_type>::Kernel>
  static auto apply(const TriangleMesh &tm,
                    const VertexPointMap &vpm,
                    const NamedParameters &np = CGAL::parameters::default_values())
  {
    using FT = typename Traits::FT;
    using Halfedge_length_tag = CGAL::dynamic_halfedge_property_t<FT>;
    using HalfedgeLengthMap = typename boost::property_map<TriangleMesh, Halfedge_length_tag>::const_type;
    return HalfedgeLengthMap{};
  }

  /*!
   * \brief Computes the cotangent weight for a given triangle.
   *
   * This function calculates the cotangent weight for the angle opposite to the edge between 
   * `v1` and `v2` in the triangle defined by vertices `v1`, `v2`, and `v3`. The cotangent weight 
   * is used in various geometric algorithms, including those involving the discrete Laplacian.
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
    using FT = typename Traits::FT;
    using VertexPointMapReference = typename boost::property_traits<VertexPointMap>::reference;

    VertexPointMapReference p_i = get(vpm, v1);
    VertexPointMapReference p_j = get(vpm, v2);
    VertexPointMapReference p_k = get(vpm, v3);
    return CGAL::Weights::cotangent(p_i, p_j, p_k, traits);
  }

  /*!
   * \brief Computes the area of a triangle defined by three vertices.
   *
   * This function calculates the area of the triangle formed by vertices `v1`, `v2`, and `v3`
   * in the triangle mesh. The area is computed using the cross product of the vectors defined 
   * by the vertices.
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
   * \return The area of the triangle.
   */
  template <typename TriangleMesh,
            typename VertexPointMap,
            typename HalfedgeLengthMap,
            typename Traits>
  static auto face_area(const TriangleMesh &tm,
                        const VertexPointMap &vpm,
                        const HalfedgeLengthMap &he_length_map,
                        typename boost::graph_traits<TriangleMesh>::vertex_descriptor v1,
                        typename boost::graph_traits<TriangleMesh>::vertex_descriptor v2,
                        typename boost::graph_traits<TriangleMesh>::vertex_descriptor v3,
                        Traits traits)
  {
    typename Traits::Construct_vector_3 construct_vector = traits.construct_vector_3_object();

    auto p_i = get(vpm, v1);
    auto p_j = get(vpm, v2);
    auto p_k = get(vpm, v3);
    
    // Compute edge vectors
    Vector_3 v_ij = construct_vector(p_i, p_j);
    Vector_3 v_ik = construct_vector(p_i, p_k);
    
    // Calculate the cross product to get the area
    Vector_3 cross = cross_product(v_ij, v_ik);
    double area = CGAL::sqrt(to_double(scalar_product(cross, cross))) * 0.5;
    return area;
  }

  /*!
   * \brief Computes the vectors of the edges of a triangle defined by three vertices.
   *
   * This function returns the edge vectors of the triangle formed by vertices `v1`, `v2`, and `v3`.
   * The vectors are computed based on the positions of the vertices in the triangle mesh.
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
   * \return An array containing the three edge vectors of the triangle.
   */
  template <typename TriangleMesh,
            typename VertexPointMap,
            typename HalfedgeLengthMap,
            typename Traits>
  static std::array<typename Traits::Vector_3, 3> face_vectors(const TriangleMesh &tm,
                                                                const VertexPointMap &vpm,
                                                                const HalfedgeLengthMap &he_length_map,
                                                                typename boost::graph_traits<TriangleMesh>::vertex_descriptor v1,
                                                                typename boost::graph_traits<TriangleMesh>::vertex_descriptor v2,
                                                                typename boost::graph_traits<TriangleMesh>::vertex_descriptor v3,
                                                                Traits traits)
  {
    using Vector_3 = typename Traits::Vector_3;
    typename Traits::Construct_vector_3 construct_vector = traits.construct_vector_3_object();

    auto p_i = get(vpm, v1);
    auto p_j = get(vpm, v2);
    auto p_k = get(vpm, v3);
    
    // Compute edge vectors
    Vector_3 v_ij = construct_vector(p_i, p_j);
    Vector_3 v_jk = construct_vector(p_j, p_k);
    Vector_3 v_ki = construct_vector(p_k, p_i);
    
    return {v_ij, v_jk, v_ki};
  }

  /*!
   * \brief Computes the length of the edge defined by two vertices.
   *
   * This function calculates the length of the edge between vertices `v1` and `v2`.
   *
   * \param tm The input triangle mesh.
   * \param vpm The property map for vertex points.
   * \param he_length_map The property map for the lengths of halfedges in the mesh.
   * \param v1 The first vertex of the edge.
   * \param v2 The second vertex of the edge.
   * \param traits The traits class providing geometric operations.
   *
   * \tparam TriangleMesh The type of the triangle mesh.
   * \tparam VertexPointMap The type of the vertex point property map.
   * \tparam HalfedgeLengthMap The type of the halfedge length property map.
   * \tparam Traits The traits class used for geometric operations.
   *
   * \return The length of the edge between `v1` and `v2`.
   */
  template <typename TriangleMesh,
            typename VertexPointMap,
            typename HalfedgeLengthMap,
            typename Traits>
  static auto length(const TriangleMesh &tm,
                     const VertexPointMap &vpm,
                     const HalfedgeLengthMap &he_length_map,
                     typename boost::graph_traits<TriangleMesh>::vertex_descriptor v1,
                     typename boost::graph_traits<TriangleMesh>::vertex_descriptor v2,
                     Traits traits)
  {
    typename Traits::Compute_squared_distance_3 squared_distance =
        traits.compute_squared_distance_3_object();
    auto p1 = get(vpm, v1);
    auto p2 = get(vpm, v2);
    return std::sqrt(to_double(squared_distance(p1, p2)));
  }
};

// Common functions calculated from halfedge lengths
struct Mollification_scheme_common {
  /*!
   * \brief Computes the cotangent of the angle opposite to edge c in a triangle.
   *
   * This function retrieves the lengths of the edges of the triangle formed by vertices
   * v1, v2, and v3, and calculates the cotangent of the angle opposite to the edge 
   * length defined by vertex v3.
   *
   * \param tm The input triangle mesh.
   * \param vpm The property map for vertex points.
   * \param he_length_map The property map for halfedge lengths.
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
   * \return The cotangent of the angle opposite to edge length c.
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

    // Retrieve edge lengths from the halfedge length map
    const FT a = get(he_length_map, halfedge(v1, v2, tm).first);
    const FT b = get(he_length_map, halfedge(v2, v3, tm).first);
    const FT c = get(he_length_map, halfedge(v3, v1, tm).first);

    // Compute the semi-perimeter of the triangle
    const FT S = (a + b + c) / 2;

    // Compute the area of the triangle using Heron's formula
    const FT area = CGAL::approximate_sqrt(S * (S - a) * (S - b) * (S - c));

    // If the area is zero, return an undefined value
    if (is_zero(area) || !std::isfinite(area)) {
      return FT(0);  // Undefined cotangent
    }

    assert(std::isfinite(area));
    // Compute the cotangent of the angle opposite to edge length c
    return (CGAL::square(a) + CGAL::square(b) - CGAL::square(c)) / (4 * area);
  }

  /*!
   * \brief Computes the area of a triangle defined by three vertices.
   *
   * This function calculates the area of the triangle formed by vertices v1, v2, and v3
   * using Heron's formula based on the lengths of its edges.
   *
   * \param tm The input triangle mesh.
   * \param vpm The property map for vertex points.
   * \param he_length_map The property map for halfedge lengths.
   * \param v1 The first vertex of the triangle.
   * \param v2 The second vertex of the triangle.
   * \param v3 The third vertex of the triangle.
   *
   * \tparam TriangleMesh The type of the triangle mesh.
   * \tparam VertexPointMap The type of the vertex point property map.
   * \tparam HalfedgeLengthMap The type of the halfedge length property map.
   * \tparam Traits The traits class used for geometric operations.
   *
   * \return The area of the triangle.
   */
  template <typename TriangleMesh,
      typename VertexPointMap,
      typename HalfedgeLengthMap,
      typename Traits>
  static auto face_area(const TriangleMesh &tm,
      const VertexPointMap &vpm,
      const HalfedgeLengthMap &he_length_map,
      typename boost::graph_traits<TriangleMesh>::vertex_descriptor v1,
      typename boost::graph_traits<TriangleMesh>::vertex_descriptor v2,
      typename boost::graph_traits<TriangleMesh>::vertex_descriptor v3,
      Traits)
  {
    auto a = get(he_length_map, halfedge(v1, v2, tm).first);
    auto b = get(he_length_map, halfedge(v2, v3, tm).first);
    auto c = get(he_length_map, halfedge(v3, v1, tm).first);
    auto S = (a + b + c) / 2;
    return CGAL::sqrt(S * (S - a) * (S - b) * (S - c));
  }

  /*!
   * \brief Computes the edge vectors of a triangle defined by three vertices.
   *
   * This function retrieves the vertices of the triangle formed by vertices v1, v2, and v3
   * and calculates the vectors corresponding to its edges.
   *
   * \param tm The input triangle mesh.
   * \param vpm The property map for vertex points.
   * \param he_length_map The property map for halfedge lengths.
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
   * \return An array containing the three edge vectors of the triangle.
   */
  template <typename TriangleMesh,
      typename VertexPointMap,
      typename HalfedgeLengthMap,
      typename Traits>
  static std::array<typename Traits::Vector_3, 3> face_vectors(const TriangleMesh &tm,
      const VertexPointMap &vpm,
      const HalfedgeLengthMap &he_length_map,
      typename boost::graph_traits<TriangleMesh>::vertex_descriptor v1,
      typename boost::graph_traits<TriangleMesh>::vertex_descriptor v2,
      typename boost::graph_traits<TriangleMesh>::vertex_descriptor v3,
      Traits traits)
  {
    typedef typename Traits::Vector_3 Vector_3;
    typedef typename Traits::Point_3 Point_3;
    typename Traits::Construct_vector_3 construct_vector = traits.construct_vector_3_object();

    auto a = get(he_length_map, halfedge(v1, v2, tm).first);
    auto b = get(he_length_map, halfedge(v2, v3, tm).first);
    auto c = get(he_length_map, halfedge(v3, v1, tm).first);
    auto S = (a + b + c) / 2;
    double angle_a = -(b * b) + c * c + a * a;
    angle_a = acos(std::clamp(angle_a / (2 * c * a), -1.0, 1.0));
    Point_3 p0(0, 0, 0);
    Point_3 p2(c * std::cos(angle_a), c * std::sin(angle_a), 0);
    Point_3 p1(a, 0, 0);
    Vector_3 v_ij = construct_vector(p0, p1);
    Vector_3 v_jk = construct_vector(p1, p2);
    Vector_3 v_ki = construct_vector(p2, p0);
    return {v_ij, v_jk, v_ki};
  }

  /*!
   * \brief Retrieves the length of an edge defined by two vertices.
   *
   * This function fetches the length of the edge between the vertices v1 and v2 
   * from the halfedge length map.
   *
   * \param tm The input triangle mesh.
   * \param vpm The property map for vertex points.
   * \param he_length_map The property map for halfedge lengths.
   * \param v1 The first vertex of the edge.
   * \param v2 The second vertex of the edge.
   * \param traits The traits class providing geometric operations.
   *
   * \tparam TriangleMesh The type of the triangle mesh.
   * \tparam VertexPointMap The type of the vertex point property map.
   * \tparam HalfedgeLengthMap The type of the halfedge length property map.
   * \tparam Traits The traits class used for geometric operations.
   *
   * \return The length of the edge between `v1` and `v2`.
   */
  template <typename TriangleMesh,
      typename VertexPointMap,
      typename HalfedgeLengthMap,
      typename Traits>
  static auto length(const TriangleMesh &tm,
      const VertexPointMap &vpm,
      const HalfedgeLengthMap &he_length_map,
      typename boost::graph_traits<TriangleMesh>::vertex_descriptor v1,
      typename boost::graph_traits<TriangleMesh>::vertex_descriptor v2,
      Traits traits)
  {
    return get(he_length_map, halfedge(v1, v2, tm).first);
  }
};

// Global mollification scheme with a constant adjustment to edge lengths
struct Mollification_scheme_global_constant : public Mollification_scheme_common {
    /*!
     * \brief Applies the constant mollification scheme to the triangle mesh.
     *
     * This function modifies the provided triangle mesh by adding a constant epsilon to each edge
     * length to ensure the strict triangle inequality holds with a tolerance of `delta`. The epsilon
     * value is computed based on the smallest edge length in the mesh and the specified `delta`.
     *
     * \param tm The input triangle mesh.
     * \param vpm The property map for vertex points.
     * \param np Named parameters, allowing for customization of parameters.
     *
     * \tparam TriangleMesh The type of the triangle mesh.
     * \tparam VertexPointMap The type of the vertex point property map.
     * \tparam NamedParameters Optional named parameters type (default is CGAL::parameters::Default_named_parameters).
     * \tparam Traits The traits class used for geometric operations.
     *
     * \return A property map containing updated halfedge lengths.
     */
    template <typename TriangleMesh,
              typename VertexPointMap,
              typename NamedParameters = CGAL::parameters::Default_named_parameters,
              typename Traits = typename Kernel_traits<typename boost::property_traits<typename boost::property_map<TriangleMesh, vertex_point_t>::const_type>::value_type>::Kernel>
    static auto apply(const TriangleMesh &tm,
                      const VertexPointMap &vpm,
                      const NamedParameters &np = CGAL::parameters::default_values())
    {
        typedef boost::graph_traits<TriangleMesh> graph_traits;
        typedef typename graph_traits::halfedge_descriptor halfedge_descriptor;
        typedef typename graph_traits::vertex_descriptor vertex_descriptor;
        typedef typename graph_traits::edge_descriptor edge_descriptor;

        // Geometric type definitions
        typedef typename Traits::Point_3 Point_3;
        typedef typename Traits::FT FT;
        typedef typename Traits::Vector_3 Vector_3;
        typename Traits::Compute_squared_distance_3 squared_distance = 
            Traits().compute_squared_distance_3_object();
        typedef typename boost::property_traits<VertexPointMap>::reference VertexPointMap_reference;

        using parameters::choose_parameter;
        using parameters::get_parameter;

        // Dynamic property to store halfedge lengths
        typedef CGAL::dynamic_halfedge_property_t<FT> Halfedge_length_tag;
        typedef typename boost::property_map<TriangleMesh, Halfedge_length_tag>::const_type HalfedgeLengthMap;

        HalfedgeLengthMap he_length_map(get(Halfedge_length_tag(), tm));

        FT min_length = std::numeric_limits<double>::max();
        FT avg_length = 0;

        // Calculate edge lengths and update min/average lengths
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
            avg_length += e_length;
        }
        avg_length /= tm.number_of_edges();

        // Determine the delta value based on the parameters
        FT delta = choose_parameter(get_parameter(np, internal_np::delta), FT{avg_length * Kdelta});

        // Compute the smallest length epsilon that can be added to all edges
        FT epsilon = 0;
        for (halfedge_descriptor hd : halfedges(tm)) {
            halfedge_descriptor hd2 = next(hd, tm);
            halfedge_descriptor hd3 = next(hd2, tm);
            FT ineq = get(he_length_map, hd2) + get(he_length_map, hd3) - get(he_length_map, hd);
            epsilon = (std::max)(epsilon, (std::max)(0., delta - ineq));
        }

        // Update edge lengths by adding epsilon
        for (halfedge_descriptor hd : halfedges(tm)) {
            put(he_length_map, hd, epsilon + get(he_length_map, hd));
        }

        return he_length_map;
    }
};

// Local mollification scheme with one-by-one adjustment of edge lengths
struct Mollification_scheme_local_one_by_one : public Mollification_scheme_common {
    /*!
     * \brief Applies a local one-by-one mollification scheme to the triangle mesh.
     *
     * This function modifies the triangle mesh by adjusting the lengths of the edges
     * in a face one at a time, ensuring that the strict triangle inequality holds with
     * a specified tolerance `delta`. The adjustments are performed sequentially for each 
     * face in the mesh, based on the current lengths of the edges and the computed epsilon.
     *
     * \param tm The input triangle mesh.
     * \param vpm The property map for vertex points.
     * \param np Named parameters, allowing for customization of parameters.
     *
     * \tparam TriangleMesh The type of the triangle mesh.
     * \tparam VertexPointMap The type of the vertex point property map.
     * \tparam NamedParameters Optional named parameters type (default is CGAL::parameters::Default_named_parameters).
     * \tparam Traits The traits class used for geometric operations.
     *
     * \return A property map containing updated halfedge lengths.
     */
    template <typename TriangleMesh,
              typename VertexPointMap,
              typename NamedParameters = CGAL::parameters::Default_named_parameters,
              typename Traits = typename Kernel_traits<typename boost::property_traits<typename boost::property_map<TriangleMesh, vertex_point_t>::const_type>::value_type>::Kernel>
    static auto apply(const TriangleMesh &tm,
                      const VertexPointMap &vpm,
                      const NamedParameters &np = CGAL::parameters::default_values())
    {
        typedef boost::graph_traits<TriangleMesh> graph_traits;
        typedef typename graph_traits::face_descriptor face_descriptor;
        typedef typename graph_traits::halfedge_descriptor halfedge_descriptor;
        typedef typename graph_traits::vertex_descriptor vertex_descriptor;
        typedef typename graph_traits::edge_descriptor edge_descriptor;

        // Geometric type definitions
        typedef typename Traits::Point_3 Point_3;
        typedef typename Traits::FT FT;
        typedef typename Traits::Vector_3 Vector_3;
        typename Traits::Compute_squared_distance_3 squared_distance = 
            Traits().compute_squared_distance_3_object();
        typedef typename boost::property_traits<VertexPointMap>::reference VertexPointMap_reference;

        using parameters::choose_parameter;
        using parameters::get_parameter;

        // Dynamic property to store halfedge lengths
        typedef CGAL::dynamic_halfedge_property_t<FT> Halfedge_length_tag;
        typedef typename boost::property_map<TriangleMesh, Halfedge_length_tag>::const_type HalfedgeLengthMap;

        HalfedgeLengthMap he_length_map(get(Halfedge_length_tag(), tm));

        FT min_length = std::numeric_limits<double>::max();
        FT avg_length = 0;

        // Calculate lengths for each edge and update min/average lengths
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
            avg_length += e_length;
        }
        avg_length /= tm.number_of_edges();

        // Determine the delta value based on the parameters
        FT delta = choose_parameter(get_parameter(np, internal_np::delta), FT{avg_length * Kdelta});

        // Compute the smallest length epsilon that can be added to all edges
        FT epsilon = 0;
        CGAL::Vertex_around_face_iterator<TriangleMesh> vbegin, vend;
        for (face_descriptor f : faces(tm)) {
            boost::tie(vbegin, vend) = vertices_around_face(halfedge(f, tm), tm);
            vertex_descriptor v0 = *(vbegin);
            vertex_descriptor v1 = *(++vbegin);
            vertex_descriptor v2 = *(++vbegin);

            halfedge_descriptor hd = halfedge(v0, v1, tm).first;
            halfedge_descriptor hd2 = halfedge(v1, v2, tm).first;
            halfedge_descriptor hd3 = halfedge(v2, v0, tm).first;

            std::array<halfedge_descriptor, 3> he = {hd, hd2, hd3};
            std::array<double, 3> L = {
                get(he_length_map, hd), get(he_length_map, hd2), get(he_length_map, hd3)};
            
            // Create indices vector and sort based on values
            std::vector<int> indices(3);
            std::iota(indices.begin(), indices.end(), 0);  // Fill indices with 0, 1, ..., size-1

            // Sort indices based on lengths
            std::sort(indices.begin(), indices.end(), [&](int i1, int i2) { return L[i1] < L[i2]; });

            // Reorder a, b, c so that c <= b <= a
            auto a = L[indices[2]];
            auto b = L[indices[1]];
            auto c = L[indices[0]];

            // Check triangle inequality and adjust if necessary
            if ((a + b - c >= delta) && (a + c - b >= delta) && (b + c - a >= delta))
                continue;

            // Mollify edges to satisfy triangle inequality
            c = std::max(c, (delta + a - b));
            b = std::max(b, c);
            a = std::max(a, b);

            // Update edge lengths
            put(he_length_map, he[indices[0]], c);
            put(he_length_map, he[indices[1]], b);
            put(he_length_map, he[indices[2]], a);
        }

        return he_length_map;
    }
};

// Local mollification scheme with one-by-one edge length interpolation
struct Mollification_scheme_local_one_by_one_interpolation : public Mollification_scheme_common {
  /*!
   * \brief Applies the one-by-one interpolation mollification scheme to a triangle mesh.
   *
   * This function computes and updates the lengths of the halfedges in the given triangle mesh
   * using a one-by-one interpolation strategy to ensure that the triangle inequality is satisfied
   * with a specified tolerance. It adjusts edge lengths based on the average edge length and a delta 
   * parameter.
   *
   * \param tm The input triangle mesh.
   * \param vpm The property map for vertex points.
   * \param np Named parameters, allowing for customization of parameters.
   *
   * \tparam TriangleMesh The type of the triangle mesh.
   * \tparam VertexPointMap The type of the vertex point property map.
   * \tparam NamedParameters Optional named parameters type (default is CGAL::parameters::Default_named_parameters).
   * \tparam Traits The traits class used for geometric operations.
   *
   * \return A property map containing updated halfedge lengths.
   */
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
    typedef typename graph_traits::face_descriptor face_descriptor;
    typedef typename graph_traits::edge_descriptor edge_descriptor;
    
    // Geometric type definitions
    typedef typename Traits::Point_3 Point_3;
    typedef typename Traits::FT FT;
    typedef typename Traits::Vector_3 Vector_3;
    typename Traits::Compute_squared_distance_3 squared_distance =
        Traits().compute_squared_distance_3_object();
    typedef typename boost::property_traits<VertexPointMap>::reference VertexPointMap_reference;

    typedef int Index;

    using parameters::choose_parameter;
    using parameters::get_parameter;

    // Define a halfedge length property map
    typedef CGAL::dynamic_halfedge_property_t<FT> Halfedge_length_tag;
    typedef typename boost::property_map<TriangleMesh, Halfedge_length_tag>::const_type
        HalfedgeLengthMap;

    HalfedgeLengthMap he_length_map(get(Halfedge_length_tag(), tm));

    FT min_length = std::numeric_limits<double>::max();
    FT avg_length = 0;

    // Compute lengths for each edge and store them
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
      avg_length += e_length;
    }
    avg_length /= tm.number_of_edges();
    
    // Calculate delta using the average length
    FT delta = avg_length * choose_parameter(get_parameter(np, internal_np::delta), FT{Kdelta});

    // Compute the smallest length epsilon to maintain triangle inequality
    FT epsilon = 0;
    auto x = 0;
    CGAL::Vertex_around_face_iterator<TriangleMesh> vbegin, vend;
    
    // Iterate through each face to ensure triangle inequality
    for (face_descriptor f : faces(tm)) {
      boost::tie(vbegin, vend) = vertices_around_face(halfedge(f, tm), tm);
      vertex_descriptor v0 = *(vbegin);
      vertex_descriptor v1 = *(++vbegin);
      vertex_descriptor v2 = *(++vbegin);

      halfedge_descriptor hd = halfedge(v0, v1, tm).first;
      halfedge_descriptor hd2 = halfedge(v1, v2, tm).first;
      halfedge_descriptor hd3 = halfedge(v2, v0, tm).first;

      std::array<halfedge_descriptor, 3> he = {hd, hd2, hd3};
      std::array<double, 3> L = {
          get(he_length_map, hd), get(he_length_map, hd2), get(he_length_map, hd3)};
      
      // Create indices vector and sort based on values
      std::vector<int> indices(3);
      std::iota(indices.begin(), indices.end(), 0);  // Fill indices with 0, 1, ..., size-1

      // Sort indices based on values in `L`
      std::sort(indices.begin(), indices.end(), [&](int i1, int i2) { return L[i1] < L[i2]; });
      
      // Reorder a, b, c so that c <= b <= a
      auto a = L[indices[2]];
      auto b = L[indices[1]];
      auto c = L[indices[0]];

      if ((a + b - c >= delta) && (a + c - b >= delta) && (b + c - a >= delta))
        continue;

      // Adjust edge lengths to satisfy triangle inequality
      auto c_prev = c;
      c = std::max(c, delta + a - b);
      c = std::max(c, delta + b - a);
      if (a - c_prev > 1e-6 * delta)
        b = std::max(c, b + (a - b) * (c - c_prev) / (a - c_prev));
      else
        b = std::max(b, c);
      a = std::max(a, b);

      // Reorder back to original order
      put(he_length_map, he[indices[0]], c);
      put(he_length_map, he[indices[1]], b);
      put(he_length_map, he[indices[2]], a);
    }

    return he_length_map;
  }
};

// Local mollification scheme with constant edge length adjustment
struct Mollification_scheme_local_constant : public Mollification_scheme_common {
  /*!
   * \brief Applies a constant mollification scheme to a triangle mesh.
   *
   * This function adjusts the lengths of the halfedges in the given triangle mesh
   * to ensure that the triangle inequality is satisfied with a uniform constant adjustment 
   * across all edges. It modifies the lengths based on the average edge length and a delta 
   * parameter.
   *
   * \param tm The input triangle mesh.
   * \param vpm The property map for vertex points.
   * \param np Named parameters, allowing for customization of parameters.
   *
   * \tparam TriangleMesh The type of the triangle mesh.
   * \tparam VertexPointMap The type of the vertex point property map.
   * \tparam NamedParameters Optional named parameters type (default is CGAL::parameters::Default_named_parameters).
   * \tparam Traits The traits class used for geometric operations.
   *
   * \return A property map containing updated halfedge lengths.
   */
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
    typedef typename graph_traits::face_descriptor face_descriptor;
    typedef typename graph_traits::edge_descriptor edge_descriptor;
    
    // Geometric type definitions
    typedef typename Traits::Point_3 Point_3;
    typedef typename Traits::FT FT;
    typedef typename Traits::Vector_3 Vector_3;
    typename Traits::Compute_squared_distance_3 squared_distance =
        Traits().compute_squared_distance_3_object();
    typedef typename boost::property_traits<VertexPointMap>::reference VertexPointMap_reference;

    typedef int Index;

    using parameters::choose_parameter;
    using parameters::get_parameter;

    // Define a halfedge length property map
    typedef CGAL::dynamic_halfedge_property_t<FT> Halfedge_length_tag;
    typedef typename boost::property_map<TriangleMesh, Halfedge_length_tag>::const_type
        HalfedgeLengthMap;

    HalfedgeLengthMap he_length_map(get(Halfedge_length_tag(), tm));

    FT min_length = std::numeric_limits<double>::max();
    FT avg_length = 0;

    // Compute lengths for each edge and store them
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
      avg_length += e_length;
    }
    avg_length /= tm.number_of_edges();
    
    // Calculate delta using the average length
    FT delta = avg_length * choose_parameter(get_parameter(np, internal_np::delta), FT{Kdelta});

    // Adjust edge lengths to satisfy triangle inequality
    CGAL::Vertex_around_face_iterator<TriangleMesh> vbegin, vend;
    for (face_descriptor f : faces(tm)) {
      boost::tie(vbegin, vend) = vertices_around_face(halfedge(f, tm), tm);
      vertex_descriptor v0 = *(vbegin);
      vertex_descriptor v1 = *(++vbegin);
      vertex_descriptor v2 = *(++vbegin);

      halfedge_descriptor hd = halfedge(v0, v1, tm).first;
      halfedge_descriptor hd2 = halfedge(v1, v2, tm).first;
      halfedge_descriptor hd3 = halfedge(v2, v0, tm).first;

      std::array<halfedge_descriptor, 3> he = {hd, hd2, hd3};
      std::array<double, 3> L = {
          get(he_length_map, hd), get(he_length_map, hd2), get(he_length_map, hd3)};
      
      // Reorder a, b, c so that c <= b <= a
      auto a = L[0];
      auto b = L[1];
      auto c = L[2];

      if ((a + b - c >= delta) && (a + c - b >= delta) && (b + c - a >= delta))
        continue;
      
      // Calculate necessary epsilon adjustment for triangle inequality
      FT ineq = a + b - c;
      FT epsilon = (std::max)(0., delta - ineq);
      ineq = c + b - a;
      epsilon = (std::max)(epsilon, (std::max)(0., delta - ineq));
      ineq = a + c - b;
      epsilon = (std::max)(epsilon, (std::max)(0., delta - ineq));
      
      // Adjust edge lengths with the calculated epsilon
      a = L[0] + epsilon;
      b = L[1] + epsilon;
      c = L[2] + epsilon;
      put(he_length_map, he[0], a);
      put(he_length_map, he[1], b);
      put(he_length_map, he[2], c);
    }

    return he_length_map;
  }
};

// Local mollification scheme with linear adjustment of edge lengths
struct Mollification_scheme_local_minimal_distance_linear : public Mollification_scheme_common {
  /*!
   * \brief Applies a linear mollification scheme to a triangle mesh.
   *
   * This function adjusts the lengths of the halfedges in the given triangle mesh
   * to ensure that the triangle inequality is satisfied with a minimal distance 
   * adjustment. It computes the edge lengths based on the vertex positions and 
   * applies a linear programming approach to find optimal edge lengths that 
   * satisfy the triangle inequality with a tolerance defined by the delta parameter.
   *
   * \param tm The input triangle mesh.
   * \param vpm The property map for vertex points.
   * \param np Named parameters, allowing for customization of parameters.
   *
   * \tparam TriangleMesh The type of the triangle mesh.
   * \tparam VertexPointMap The type of the vertex point property map.
   * \tparam NamedParameters Optional named parameters type (default is CGAL::parameters::Default_named_parameters).
   * \tparam Traits The traits class used for geometric operations.
   *
   * \return A property map containing updated halfedge lengths.
   */
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
    typedef typename graph_traits::face_descriptor face_descriptor;
    typedef typename graph_traits::edge_descriptor edge_descriptor;

    // Geometric type definitions
    typedef typename Traits::Point_3 Point_3;
    typedef typename Traits::FT FT;
    typedef typename Traits::Vector_3 Vector_3;
    typename Traits::Compute_squared_distance_3 squared_distance =
        Traits().compute_squared_distance_3_object();
    typedef typename boost::property_traits<VertexPointMap>::reference VertexPointMap_reference;

    using parameters::choose_parameter;
    using parameters::get_parameter;

    typedef CGAL::MP_Float ET;

    // Program and solution types
    typedef CGAL::Quadratic_program<double> Program;
    typedef CGAL::Quadratic_program_solution<ET> Solution;
    typedef typename Solution::Variable_value_iterator Variable_value_iterator;
    typedef CGAL::Real_embeddable_traits<typename Variable_value_iterator::value_type> RE_traits;
    typename RE_traits::To_double to_double;

    // Define a halfedge length property map
    typedef CGAL::dynamic_halfedge_property_t<FT> Halfedge_length_tag;
    typedef typename boost::property_map<TriangleMesh, Halfedge_length_tag>::const_type
        HalfedgeLengthMap;

    HalfedgeLengthMap he_length_map(get(Halfedge_length_tag(), tm));

    double min_length = std::numeric_limits<double>::max();
    double max_length = -1;
    FT avg_length = 0;

    // Compute lengths for each edge and store them
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
      max_length = CGAL::max(max_length, e_length);
      avg_length += e_length;
    }
    avg_length /= tm.number_of_edges();
    
    // Calculate delta using the average length
    FT delta = avg_length * choose_parameter(get_parameter(np, internal_np::delta), FT{Kdelta});

    // Adjust edge lengths to satisfy triangle inequality
    CGAL::Vertex_around_face_iterator<TriangleMesh> vbegin, vend;
    for (face_descriptor f : faces(tm)) {
      boost::tie(vbegin, vend) = vertices_around_face(halfedge(f, tm), tm);
      vertex_descriptor v0 = *(vbegin);
      vertex_descriptor v1 = *(++vbegin);
      vertex_descriptor v2 = *(++vbegin);

      halfedge_descriptor hd = halfedge(v0, v1, tm).first;
      halfedge_descriptor hd2 = halfedge(v1, v2, tm).first;
      halfedge_descriptor hd3 = halfedge(v2, v0, tm).first;

      std::array<halfedge_descriptor, 3> he = {hd, hd2, hd3};
      std::array<double, 3> L = {
          get(he_length_map, hd), get(he_length_map, hd2), get(he_length_map, hd3)};

      // Reorder a, b, c so that c <= b <= a
      auto a = L[0];
      auto b = L[1];
      auto c = L[2];

      if ((a + b - c >= delta) && (a + c - b >= delta) && (b + c - a >= delta))
        continue;

      // Set up the linear program to adjust lengths
      Program lp(CGAL::LARGER, true, 0, false, max_length);
      const int L0 = 0;
      const int L1 = 1;
      const int L2 = 2;
      lp.set_a(L0, 0, 1.0);
      lp.set_a(L1, 0, 1.0);
      lp.set_a(L2, 0, -1.0);
      lp.set_b(0, delta);  // L0 + L1 - L2 >= delta
      lp.set_a(L0, 1, 1.0);
      lp.set_a(L1, 1, -1.0);
      lp.set_a(L2, 1, 1.0);
      lp.set_b(1, delta);  // L0 - L1 + L2 >= delta
      lp.set_a(L0, 2, -1.0);
      lp.set_a(L1, 2, 1.0);
      lp.set_a(L2, 2, 1.0);
      lp.set_b(2, delta);  // -L0 + L1 + L2 >= delta
      lp.set_a(L0, 3, 1.0);
      lp.set_b(3, L[0]);  // L0 >= L0_old
      lp.set_a(L1, 4, 1.0);
      lp.set_b(4, L[1]);  // L1 >= L1_old
      lp.set_a(L2, 5, 1.0);
      lp.set_b(5, L[2]);  // L2 >= L2_old
      lp.set_c(L0, 1.0);
      lp.set_c(L1, 1.0);
      lp.set_c(L2, 1.0);  // L0 + L1 + L2

      // Solve the program, using ET as the exact type
      Solution solution = CGAL::solve_linear_program(lp, ET());
      // Get variables
      auto X = solution.variable_values_begin();

      L = {to_double(X[0]), to_double(X[1]), to_double(X[2])};
      put(he_length_map, he[0], L[0]);
      put(he_length_map, he[1], L[1]);
      put(he_length_map, he[2], L[2]);
      a = L[0];
      b = L[1];
      c = L[2];
      assert((a + b - c >= delta - 1e-5) && (a + c - b >= delta - 1e-5) &&
             (b + c - a >= delta - 1e-5));
    }

    return he_length_map;
  }
};

// Local mollification scheme with quadratic adjustment of edge lengths
struct Mollification_scheme_local_minimal_distance_quadratic : public Mollification_scheme_common {
  /*!
   * \brief Applies a quadratic mollification scheme to a triangle mesh.
   *
   * This function adjusts the lengths of the halfedges in the given triangle mesh
   * to ensure that the triangle inequality is satisfied with a quadratic adjustment 
   * based on the minimal distance between edges. It modifies the lengths using a 
   * quadratic programming approach that minimizes the squared differences from the 
   * original edge lengths while satisfying the triangle inequality with a defined 
   * tolerance.
   *
   * \param tm The input triangle mesh.
   * \param vpm The property map for vertex points.
   * \param np Named parameters, allowing for customization of parameters.
   *
   * \tparam TriangleMesh The type of the triangle mesh.
   * \tparam VertexPointMap The type of the vertex point property map.
   * \tparam NamedParameters Optional named parameters type (default is CGAL::parameters::Default_named_parameters).
   * \tparam Traits The traits class used for geometric operations.
   *
   * \return A property map containing updated halfedge lengths.
   */
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
    typedef typename graph_traits::face_descriptor face_descriptor;
    typedef typename graph_traits::edge_descriptor edge_descriptor;

    // Geometric type definitions
    typedef typename Traits::Point_3 Point_3;
    typedef typename Traits::FT FT;
    typedef typename Traits::Vector_3 Vector_3;
    typename Traits::Compute_squared_distance_3 squared_distance =
        Traits().compute_squared_distance_3_object();
    typedef typename boost::property_traits<VertexPointMap>::reference VertexPointMap_reference;

    using parameters::choose_parameter;
    using parameters::get_parameter;

    typedef CGAL::MP_Float ET;

    // Program and solution types
    typedef CGAL::Quadratic_program<double> Program;
    typedef CGAL::Quadratic_program_solution<ET> Solution;
    typedef typename Solution::Variable_value_iterator Variable_value_iterator;
    typedef CGAL::Real_embeddable_traits<typename Variable_value_iterator::value_type> RE_traits;
    typename RE_traits::To_double to_double;

    // Define a halfedge length property map
    typedef CGAL::dynamic_halfedge_property_t<FT> Halfedge_length_tag;
    typedef typename boost::property_map<TriangleMesh, Halfedge_length_tag>::const_type
        HalfedgeLengthMap;

    HalfedgeLengthMap he_length_map(get(Halfedge_length_tag(), tm));

    double min_length = std::numeric_limits<double>::max();
    double max_length = -1;
    FT avg_length = 0;

    // Compute lengths for each edge and store them
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
      max_length = CGAL::max(max_length, e_length);
      avg_length += e_length;
    }
    avg_length /= tm.number_of_edges();
    
    // Calculate delta using the average length
    FT delta = avg_length * choose_parameter(get_parameter(np, internal_np::delta), FT{Kdelta});

    // Adjust edge lengths to satisfy triangle inequality
    CGAL::Vertex_around_face_iterator<TriangleMesh> vbegin, vend;
    for (face_descriptor f : faces(tm)) {
      boost::tie(vbegin, vend) = vertices_around_face(halfedge(f, tm), tm);
      vertex_descriptor v0 = *(vbegin);
      vertex_descriptor v1 = *(++vbegin);
      vertex_descriptor v2 = *(++vbegin);

      halfedge_descriptor hd = halfedge(v0, v1, tm).first;
      halfedge_descriptor hd2 = halfedge(v1, v2, tm).first;
      halfedge_descriptor hd3 = halfedge(v2, v0, tm).first;

      std::array<halfedge_descriptor, 3> he = {hd, hd2, hd3};
      std::array<double, 3> L = {
          get(he_length_map, hd), get(he_length_map, hd2), get(he_length_map, hd3)};

      // Reorder a, b, c so that c <= b <= a
      auto a = L[0];
      auto b = L[1];
      auto c = L[2];

      if ((a + b - c >= delta) && (a + c - b >= delta) && (b + c - a >= delta))
        continue;

      // Set up the quadratic program to adjust lengths
      Program lp(CGAL::LARGER, true, 0, false, max_length);
      const int L0 = 0;
      const int L1 = 1;
      const int L2 = 2;
      lp.set_a(L0, 0, 1.0);
      lp.set_a(L1, 0, 1.0);
      lp.set_a(L2, 0, -1.0);
      lp.set_b(0, delta);  // L0 + L1 - L2 >= delta
      lp.set_a(L0, 1, 1.0);
      lp.set_a(L1, 1, -1.0);
      lp.set_a(L2, 1, 1.0);
      lp.set_b(1, delta);  // L0 - L1 + L2 >= delta
      lp.set_a(L0, 2, -1.0);
      lp.set_a(L1, 2, 1.0);
      lp.set_a(L2, 2, 1.0);
      lp.set_b(2, delta);  // -L0 + L1 + L2 >= delta
      lp.set_a(L0, 3, 1.0);
      lp.set_b(3, L[0]);  // L0 >= L0_old
      lp.set_a(L1, 4, 1.0);
      lp.set_b(4, L[1]);  // L1 >= L1_old
      lp.set_a(L2, 5, 1.0);
      lp.set_b(5, L[2]);  // L2 >= L2_old
      lp.set_c(L0, 1.0);
      lp.set_c(L1, 1.0);
      lp.set_c(L2, 1.0);  // L0 + L1 + L2

      // Solve the program, using ET as the exact type
      Solution solution = CGAL::solve_linear_program(lp, ET());
      // Get variables
      auto X = solution.variable_values_begin();

      L = {to_double(X[0]), to_double(X[1]), to_double(X[2])};
      put(he_length_map, he[0], L[0]);
      put(he_length_map, he[1], L[1]);
      put(he_length_map, he[2], L[2]);
      a = L[0];
      b = L[1];
      c = L[2];
      assert((a + b - c >= delta - 1e-5) && (a + c - b >= delta - 1e-5) &&
             (b + c - a >= delta - 1e-5));
    }

    return he_length_map;
  }
};

struct Mollification_scheme_global_minimal_distance_linear : public Mollification_scheme_common {
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
    typedef typename graph_traits::face_descriptor face_descriptor;
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

    double min_length = std::numeric_limits<double>::max();
    double max_length = -1;
    FT avg_length = 0;
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
      max_length = CGAL::max(max_length, e_length);
      avg_length += e_length;
    }
    avg_length /= tm.number_of_edges();
    // TODO: add threshold parameter instead of constant 1e-4
    FT delta = avg_length * choose_parameter(get_parameter(np, internal_np::delta), FT{Kdelta});

    // compute smallest length epsilon we can add to
    // all edges to ensure that the strict triangle
    // inequality holds with a tolerance of delta
    FT epsilon = 0;
    CGAL::Vertex_around_face_iterator<TriangleMesh> vbegin, vend;
    typedef CGAL::dynamic_edge_property_t<Index> Edge_property_tag;
    typedef typename boost::property_map<TriangleMesh, Edge_property_tag>::type Edge_id_map;
    auto edge_id_map = get(Edge_property_tag(), tm);
    Index edge_i = 0;

    auto edges_cnt = num_edges(tm);
    auto faces_cnt = num_faces(tm);
    Eigen::SparseMatrix<double> constraint_matrix(3 * faces_cnt + edges_cnt, edges_cnt);
    Eigen::SparseMatrix<double> objective_matrix(edges_cnt, edges_cnt);
    constraint_matrix.setZero();
    objective_matrix.setZero();
    osqp::OsqpInstance instance;
    instance.objective_matrix = objective_matrix;
    instance.objective_vector.resize(edges_cnt);
    instance.constraint_matrix = constraint_matrix;
    instance.lower_bounds.resize(edges_cnt + 3 * faces_cnt);
    instance.upper_bounds.resize(edges_cnt + 3 * faces_cnt);
    for (edge_descriptor ed : edges(tm)) {
      auto L = get(he_length_map, halfedge(ed, tm));
      // L0 + L1 + ...
      // instance.objective_matrix.coeffRef(edge_i, edge_i) = 2.0;
      instance.objective_vector(edge_i) = 1.0;
      // L_0 <= L <= infy
      instance.constraint_matrix.coeffRef(edge_i, edge_i) = 1.0;
      instance.lower_bounds(edge_i) = L;
      instance.upper_bounds(edge_i) = std::numeric_limits<double>::infinity();
      put(edge_id_map, ed, edge_i++);
    }

    for (face_descriptor f : faces(tm)) {
      boost::tie(vbegin, vend) = vertices_around_face(halfedge(f, tm), tm);
      vertex_descriptor v0 = *(vbegin);
      vertex_descriptor v1 = *(++vbegin);
      vertex_descriptor v2 = *(++vbegin);

      halfedge_descriptor hd = halfedge(v0, v1, tm).first;
      halfedge_descriptor hd2 = halfedge(v1, v2, tm).first;
      halfedge_descriptor hd3 = halfedge(v2, v0, tm).first;

      Index L0 = get(edge_id_map, edge(hd, tm));
      Index L1 = get(edge_id_map, edge(hd2, tm));
      Index L2 = get(edge_id_map, edge(hd3, tm));
      assert(L0 != L1 && L1 != L2);
      std::array<double, 3> L = {
          get(he_length_map, hd), get(he_length_map, hd2), get(he_length_map, hd3)};
      // delta <= L0 + L1 - L2
      instance.constraint_matrix.coeffRef(edge_i, L0) = 1.0;
      instance.constraint_matrix.coeffRef(edge_i, L1) = 1.0;
      instance.constraint_matrix.coeffRef(edge_i, L2) = -1.0;
      instance.lower_bounds(edge_i) = delta;
      instance.upper_bounds(edge_i) = std::numeric_limits<double>::infinity();
      ++edge_i;
      // delta <= L0 - L1 + L2
      instance.constraint_matrix.coeffRef(edge_i, L0) = 1.0;
      instance.constraint_matrix.coeffRef(edge_i, L1) = -1.0;
      instance.constraint_matrix.coeffRef(edge_i, L2) = 1.0;
      instance.lower_bounds(edge_i) = delta;
      instance.upper_bounds(edge_i) = std::numeric_limits<double>::infinity();
      ++edge_i;
      // delta <= -L0 + L1 + L2
      instance.constraint_matrix.coeffRef(edge_i, L0) = -1.0;
      instance.constraint_matrix.coeffRef(edge_i, L1) = 1.0;
      instance.constraint_matrix.coeffRef(edge_i, L2) = 1.0;
      instance.lower_bounds(edge_i) = delta;
      instance.upper_bounds(edge_i) = std::numeric_limits<double>::infinity();
      ++edge_i;
    }
    instance.constraint_matrix.makeCompressed();
    instance.objective_matrix.makeCompressed();

    osqp::OsqpSolver solver;
    osqp::OsqpSettings settings;
    settings.verbose = false;
    settings.warm_start = true;
    // Edit settings if appropriate.
    auto status = solver.Init(instance, settings);
    assert(status.ok());
    // Assuming status.ok().
    osqp::OsqpExitCode exit_code = solver.Solve();
    // Assuming exit_code == OsqpExitCode::kOptimal.
    double optimal_objective = solver.objective_value();
    Eigen::VectorXd optimal_solution = solver.primal_solution();
    assert(exit_code == osqp::OsqpExitCode::kOptimal);

    for (edge_descriptor ed : edges(tm)) {
      halfedge_descriptor hd1 = halfedge(ed, tm);
      halfedge_descriptor hd2 = opposite(hd1, tm);
      Index idx = get(edge_id_map, ed);
      put(he_length_map, hd1, optimal_solution(idx));
      put(he_length_map, hd2, optimal_solution(idx));
    }

    return he_length_map;
  }
};

struct Mollification_scheme_global_minimal_distance_quadratic : public Mollification_scheme_common {
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
    typedef typename graph_traits::face_descriptor face_descriptor;
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

    double min_length = std::numeric_limits<double>::max();
    double max_length = -1;
    FT avg_length = 0;
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
      max_length = CGAL::max(max_length, e_length);
      avg_length += e_length;
    }
    avg_length /= tm.number_of_edges();
    // TODO: add threshold parameter instead of constant 1e-4
    FT delta = avg_length * choose_parameter(get_parameter(np, internal_np::delta), FT{Kdelta});

    // compute smallest length epsilon we can add to
    // all edges to ensure that the strict triangle
    // inequality holds with a tolerance of delta
    FT epsilon = 0;
    CGAL::Vertex_around_face_iterator<TriangleMesh> vbegin, vend;
    typedef CGAL::dynamic_edge_property_t<Index> Edge_property_tag;
    typedef typename boost::property_map<TriangleMesh, Edge_property_tag>::type Edge_id_map;
    auto edge_id_map = get(Edge_property_tag(), tm);
    Index edge_i = 0;

    auto edges_cnt = num_edges(tm);
    auto faces_cnt = num_faces(tm);
    Eigen::SparseMatrix<double> constraint_matrix(3 * faces_cnt + edges_cnt, edges_cnt);
    Eigen::SparseMatrix<double> objective_matrix(edges_cnt, edges_cnt);
    constraint_matrix.setZero();
    objective_matrix.setZero();
    osqp::OsqpInstance instance;
    instance.objective_matrix = objective_matrix;
    instance.objective_vector.resize(edges_cnt);
    instance.constraint_matrix = constraint_matrix;
    instance.lower_bounds.resize(edges_cnt + 3 * faces_cnt);
    instance.upper_bounds.resize(edges_cnt + 3 * faces_cnt);
    for (edge_descriptor ed : edges(tm)) {
      auto L = get(he_length_map, halfedge(ed, tm));
      // (L0 - L0_old)^2 = L0^2 - 2 L0 L0_old + L0_old^2
      // NOTE: constant L0_old^2 is ignored.
      instance.objective_matrix.coeffRef(edge_i, edge_i) = 2.0;
      instance.objective_vector(edge_i) = -2.0 * L;
      // L_0 <= L <= infy
      instance.constraint_matrix.coeffRef(edge_i, edge_i) = 1.0;
      instance.lower_bounds(edge_i) = L;
      instance.upper_bounds(edge_i) = std::numeric_limits<double>::infinity();
      put(edge_id_map, ed, edge_i++);
    }

    for (face_descriptor f : faces(tm)) {
      boost::tie(vbegin, vend) = vertices_around_face(halfedge(f, tm), tm);
      vertex_descriptor v0 = *(vbegin);
      vertex_descriptor v1 = *(++vbegin);
      vertex_descriptor v2 = *(++vbegin);

      halfedge_descriptor hd = halfedge(v0, v1, tm).first;
      halfedge_descriptor hd2 = halfedge(v1, v2, tm).first;
      halfedge_descriptor hd3 = halfedge(v2, v0, tm).first;

      Index L0 = get(edge_id_map, edge(hd, tm));
      Index L1 = get(edge_id_map, edge(hd2, tm));
      Index L2 = get(edge_id_map, edge(hd3, tm));
      assert(L0 != L1 && L1 != L2);
      std::array<double, 3> L = {
          get(he_length_map, hd), get(he_length_map, hd2), get(he_length_map, hd3)};
      // delta <= L0 + L1 - L2
      instance.constraint_matrix.coeffRef(edge_i, L0) = 1.0;
      instance.constraint_matrix.coeffRef(edge_i, L1) = 1.0;
      instance.constraint_matrix.coeffRef(edge_i, L2) = -1.0;
      instance.lower_bounds(edge_i) = delta;
      instance.upper_bounds(edge_i) = std::numeric_limits<double>::infinity();
      ++edge_i;
      // delta <= L0 - L1 + L2
      instance.constraint_matrix.coeffRef(edge_i, L0) = 1.0;
      instance.constraint_matrix.coeffRef(edge_i, L1) = -1.0;
      instance.constraint_matrix.coeffRef(edge_i, L2) = 1.0;
      instance.lower_bounds(edge_i) = delta;
      instance.upper_bounds(edge_i) = std::numeric_limits<double>::infinity();
      ++edge_i;
      // delta <= -L0 + L1 + L2
      instance.constraint_matrix.coeffRef(edge_i, L0) = -1.0;
      instance.constraint_matrix.coeffRef(edge_i, L1) = 1.0;
      instance.constraint_matrix.coeffRef(edge_i, L2) = 1.0;
      instance.lower_bounds(edge_i) = delta;
      instance.upper_bounds(edge_i) = std::numeric_limits<double>::infinity();
      ++edge_i;
    }
    instance.constraint_matrix.makeCompressed();
    instance.objective_matrix.makeCompressed();

    osqp::OsqpSolver solver;
    osqp::OsqpSettings settings;
    settings.verbose = false;
    // Edit settings if appropriate.
    auto status = solver.Init(instance, settings);

    // Assuming status.ok().
    osqp::OsqpExitCode exit_code = solver.Solve();
    // Assuming exit_code == OsqpExitCode::kOptimal.
    double optimal_objective = solver.objective_value();
    Eigen::VectorXd optimal_solution = solver.primal_solution();

    for (edge_descriptor ed : edges(tm)) {
      halfedge_descriptor hd1 = halfedge(ed, tm);
      halfedge_descriptor hd2 = opposite(hd1, tm);
      Index idx = get(edge_id_map, ed);
      put(he_length_map, hd1, optimal_solution(idx));
      put(he_length_map, hd2, optimal_solution(idx));
    }

    return he_length_map;
  }
};
}  // namespace Heat_method_3
}  // namespace CGAL

#include <CGAL/enable_warnings.h>
#endif  // CGAL_INTRINSIC_MOLLIFICATION_3_H
