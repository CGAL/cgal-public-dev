// Copyright (c) 2019 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is a part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY, AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Dmitry Anisimov, Simon Giraudot, Pierre Alliez, Florent Lafarge, and Andreas Fabri
//

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_2_DEPR_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_2_DEPR_H

// Example:

  // void compute_visibility_2(
  //   const FT visibility_scale_2) {

  //   if (!m_boundaries_detected) return;
  //   if (m_partition_2.empty()) return;

  //   Points points;
  //   points.reserve(m_interior_points.size() + m_boundary_points.size());
  //   for (const auto& idx : m_interior_points)
  //     points.push_back(idx);
  //   for (const auto& idx : m_boundary_points)
  //     points.push_back(idx);

  //   Visibility_query visibility_query(
  //     points, visibility_scale_2, m_data.point_map_2);
  //   const Visibility_2 visibility(
  //     points,
  //     visibility_query,
  //     m_data.point_map_2,
  //     m_data.visibility_map_d);
  //   visibility.compute(m_partition_2);
  // }

// STL includes.
#include <vector>
#include <utility>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/barycenter.h>

// LOD includes.
#include <CGAL/Levels_of_detail/enum.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange,
  typename NeighborQuery,
  typename PointMap2,
  typename VisibilityMap>
  class Visibility_2_depr {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Neighbor_query = NeighborQuery;
    using Point_map_2 = PointMap2;
    using Visibility_map = VisibilityMap;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Vector_2 = typename Traits::Vector_2;
    using Triangle_2 = typename Traits::Triangle_2;

    using Partition_2 = internal::Partition_2<Traits>;
    using Face = typename Partition_2::Face;

    Visibility_2_depr(
      const Input_range& input_range,
      const Neighbor_query& neighbor_query,
      const Point_map_2& point_map_2,
      const Visibility_map& visibility_map) :
    m_input_range(input_range),
    m_neighbor_query(neighbor_query),
    m_point_map_2(point_map_2),
    m_visibility_map(visibility_map),
    m_num_probes(100)
    { }

    void compute(Partition_2& partition) const {

      if (partition.empty()) return;
      label_exterior_faces(partition.faces);
      for (auto& face : partition.faces)
        if (face.exterior) {
          face.visibility = Visibility_label::OUTSIDE;
          face.inside = FT(0);
          face.outside = FT(1);
        }
        else compute_monte_carlo_label(face);
    }

  private:
    const Input_range& m_input_range;
    const Neighbor_query& m_neighbor_query;
    const Point_map_2& m_point_map_2;
    const Visibility_map& m_visibility_map;
    const std::size_t m_num_probes;

    void label_exterior_faces(
      std::vector<Face>& faces) const {

      // std::vector<Point_2> bbox;
      // internal::bounding_box_2(m_input_range, m_point_map_2, bbox);
      for (auto& face : faces) {

        face.exterior = false;
        const auto& neighbors = face.neighbors;
        for (const int idx : neighbors) {
          if (idx < 0) {
            face.exterior = true;
            break;
          }
        }
        if (face.exterior)
          continue;
        // if (is_outside_bbox(face, bbox))
          // face.exterior = true;
      }
    }

    // THIS CODE IS WORKING BUT I FOUND A BETTER AND FASTER WAY TO DO THE SAME.
    // STILL, IT CAN BE INCLUDED FOR SLIGHTLY BETTER QUALITY.
    // JUST OUTCOMMENT ALL LINES BELOW AND ABOVE.

    // bool is_outside_bbox(
    //   const Face& face, const std::vector<Point_2>& bbox) const {

    //   const Point_2& bottom = bbox[0];
    //   const Point_2& top = bbox[2];
    //   const auto& tri = face.base.delaunay;

    //   for (auto fh = tri.finite_faces_begin();
    //   fh != tri.finite_faces_end(); ++fh) {

    //     const Point_2 b = CGAL::barycenter(
    //       fh->vertex(0)->point(), FT(1),
    //       fh->vertex(1)->point(), FT(1),
    //       fh->vertex(2)->point(), FT(1));

    //     if (b.x() < bottom.x() || b.y() < bottom.y() ||
    //         b.x() > top.x() || b.y() > top.y())
    //       return true;
    //   }
    //   return false;
    // }

    void compute_monte_carlo_label(Face& face) const {

      std::vector< std::pair<Triangle_2, FT> > probability;
      create_probability(face, probability);

      const FT mean_value = compute_mean_value(probability);
      CGAL_assertion(mean_value >= FT(0) && mean_value <= FT(1));
      if (mean_value > FT(1) / FT(2))
        face.visibility = Visibility_label::INSIDE;
      else
        face.visibility = Visibility_label::OUTSIDE;

      face.inside = mean_value;
      face.outside = FT(1) - mean_value;
    }

    void create_probability(
      const Face& face,
      std::vector< std::pair<Triangle_2, FT> >& probability) const {

      const auto& tri = face.base.delaunay;
      probability.clear();
      probability.reserve(tri.number_of_faces());

      FT area = FT(0); Triangle_2 triangle;
      for (auto fh = tri.finite_faces_begin();
      fh != tri.finite_faces_end(); ++fh) {

        triangle = Triangle_2(
          fh->vertex(0)->point(),
          fh->vertex(1)->point(),
          fh->vertex(2)->point());
        probability.push_back(std::make_pair(triangle, area));
        area += CGAL::abs(triangle.area());
      }
      probability.push_back(std::make_pair(Triangle_2(), area));
    }

    FT compute_mean_value(
      const std::vector< std::pair<Triangle_2, FT> >& probability) const {

      Point_2 point;
      FT mean_value = FT(0);
      std::vector<std::size_t> neighbors;
      for (std::size_t i = 0; i < m_num_probes; ++i) {

        compute_random_point_in_triangles(probability, point);
        m_neighbor_query(point, neighbors);
        FT value = FT(0);
        if (neighbors.size() > 0)
          value += get_function_value(point, neighbors);
        mean_value += value;
      }

      mean_value /= static_cast<FT>(m_num_probes);
      return mean_value;
    }

    void compute_random_point_in_triangles(
      const std::vector< std::pair<Triangle_2, FT> >& probability,
      Point_2& point) const {

      const FT key = static_cast<FT>(
        CGAL::to_double(probability.back().second) *
        (rand() / static_cast<double>(RAND_MAX)));

      for (std::size_t i = 0; i < probability.size() - 1; ++i) {
        if (probability[i].second < key && key <= probability[i + 1].second) {
          compute_random_point_in_triangle(probability[i].first, point);
          return;
        }
      }
      std::cerr <<
        "Error (compute_random_point_in_triangles()): probability is out of range!"
      << std::endl;
      point = Point_2(FT(0), FT(0));
    }

    void compute_random_point_in_triangle(
      const Triangle_2& triangle,
      Point_2& point) const {

      const Vector_2 v(triangle[0], triangle[1]);
      const Vector_2 w(triangle[0], triangle[2]);
      point = triangle[0];
      double rv = rand() / static_cast<double>(RAND_MAX);
      double rw = rand() / static_cast<double>(RAND_MAX);

      if (rv + rw > 1.0) {
        rv = 1.0 - rv;
        rw = 1.0 - rw;
      }

      const FT bv = static_cast<FT>(rv);
      const FT bw = static_cast<FT>(rw);
      point += bv * v;
      point += bw * w;
    }

    FT get_function_value(
      const Point_2& p,
      const std::vector<std::size_t>& neighbors) const {

      std::size_t final_idx = 0; FT min_dist = internal::max_value<FT>();
      for (std::size_t i = 0; i < neighbors.size(); ++i) {
        const Point_2& q =
        get(m_point_map_2, *(m_input_range.begin() + neighbors[i]));
        const FT sq_dist = CGAL::squared_distance(p, q);
        if (sq_dist < min_dist) {
          min_dist = sq_dist;
          final_idx = i;
        }
      }

      const double value =
      get(m_visibility_map, *(m_input_range.begin() + neighbors[final_idx]));
      return static_cast<FT>(value);
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_2_DEPR_H
