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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_3_EXP_3_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_3_EXP_3_H

// Example:

// void compute_visibility_3_exp_3(
//   const FT visibility_scale_3) {

//   using Visibility_3_exp_3 = internal::Visibility_3_exp_3<
//   Traits, Points_3, Point_map_3>;

//   if (m_partition_3.empty()) return;
//   Visibility_3_exp_3 visibility(
//     m_cluster,
//     m_data.point_map_3,
//     m_building,
//     m_roof_points_3);
//   visibility.compute(m_partition_3);

//   apply_graphcut_3(
//     m_data.parameters.buildings.graphcut_beta_3);
//   compute_roofs_and_corresponding_walls();

//   std::cout << "visibility finished" << std::endl;
// }

// STL includes.
#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <memory>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/point_generators_3.h>

// Internal includes.
#include <CGAL/Levels_of_detail/enum.h>
#include <CGAL/Levels_of_detail/internal/struct.h>
#include <CGAL/Levels_of_detail/internal/utils.h>

// Spatial search.
#include <CGAL/Levels_of_detail/internal/Spatial_search/K_neighbor_query.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits,
  typename InputRange,
  typename PointMap>
  class Visibility_3_exp_3 {

  public:
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Vector_3 = typename Traits::Vector_3;
    using Line_3 = typename Traits::Line_3;
    using Plane_3 = typename Traits::Plane_3;
    using Triangle_2 = typename Traits::Triangle_2;

    using Indices = std::vector<std::size_t>;
    using Partition_3 = internal::Partition_3<Traits>;
    using Stats = std::pair<FT, FT>;
    using Face = typename Partition_3::Face;
    using Building = internal::Building<Traits>;

    using Local_traits = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Local_point_3 = typename Local_traits::Point_3;
    using Delaunay_3 = CGAL::Delaunay_triangulation_3<Local_traits>;
    using Polyhedron = CGAL::Polyhedron_3<Local_traits>;
    using Generator = CGAL::Random_points_in_triangle_mesh_3<Polyhedron>;

    using Saver = Saver<Traits>;
    using Color = CGAL::Color;

    using Pair = std::pair<Point_2, FT>;
    using Point_map_2 = CGAL::First_of_pair_property_map<Pair>;
    using K_neighbor_query =
    internal::K_neighbor_query<Traits, std::vector<Pair>, Point_map_2>;

    Visibility_3_exp_3(
      const Input_range& input_range,
      const Point_map& point_map,
      const Building& building,
      const std::vector<Indices>& roof_points_3) :
    m_input_range(input_range),
    m_point_map(point_map),
    m_building(building),
    m_roof_points_3(roof_points_3),
    m_num_samples(100), // num samples per tetrahedron
    m_k(3)
    { }

    void compute(Partition_3& partition) {

      m_samples.clear();
      create_tree();

      if (partition.empty()) return;
      label_exterior_faces(partition.faces);
      for (auto& face : partition.faces) {
        if (face.exterior) {
          face.visibility = Visibility_label::OUTSIDE;
          face.inside = FT(0);
          face.outside = FT(1);
        } else compute_face_label(face);
      }

      // Saver saver;

      // const std::string sname =
      // "/Users/monet/Documents/lod/logs/buildings/samples";
      // const Color color(0, 0, 0);
      // saver.export_points(m_samples, color, sname);

      // std::vector<Point_3> queries;
      // queries.reserve(m_queries.size());
      // for (const auto& p : m_queries)
      //   queries.push_back(Point_3(p.x(), p.y(), FT(0)));
      // const std::string qname =
      // "/Users/monet/Documents/lod/logs/buildings/queries";
      // const Color color(0, 0, 0);
      // saver.export_points(queries, color, qname);
    }

  private:
    const Input_range& m_input_range;
    const Point_map& m_point_map;
    const Building& m_building;
    const std::vector<Indices>& m_roof_points_3;
    const std::size_t m_num_samples;
    std::vector<Point_3> m_samples;
    const std::size_t m_k;

    std::vector< std::pair<Point_2, FT> > m_queries;
    Point_map_2 m_point_map_2;
    std::shared_ptr<K_neighbor_query> m_neighbor_query_ptr;

    void create_tree() {

      std::size_t num_points = 0;
      for (const auto& region : m_roof_points_3)
        num_points += region.size();

      m_queries.clear();
      m_queries.reserve(num_points);
      for (const auto& region : m_roof_points_3) {
        for (const std::size_t idx : region) {
          const Point_3& p = get(m_point_map, *(m_input_range.begin() + idx));
          const Point_2 q = internal::point_2_from_point_3(p);
          m_queries.push_back(std::make_pair(q, p.z()));
        }
      }

      m_neighbor_query_ptr = std::make_shared<K_neighbor_query>(
        m_queries, FT(m_k), m_point_map_2);
    }

    bool is_horizontal(const std::vector<std::size_t>& points) const {

	    Plane_3 plane;
      const Vector_3 orth = Vector_3(FT(0), FT(0), FT(1));
      internal::plane_from_points_3(points, m_point_map, plane);
      Vector_3 norm = plane.orthogonal_vector();
      internal::normalize(norm);

	    const FT angle = internal::angle_3d(orth, norm);
      return (angle >= FT(0) && angle <= FT(1));
    }

    void label_exterior_faces(
      std::vector<Face>& faces) const {

      for (auto& face : faces) {
        face.exterior = false;
        const auto& neighbors = face.neighbors;
        for (const int idx : neighbors) {
          if (idx < 0) {
            face.exterior = true;
            break;
          }
        }
      }
    }

    void compute_face_label(Face& face) {

      const Stats stats = estimate_in_out_values(face);
      CGAL_assertion(stats.first  >= FT(0) && stats.first  <= FT(1));
      CGAL_assertion(stats.second >= FT(0) && stats.second <= FT(1));
      CGAL_assertion(
        CGAL::abs(stats.first + stats.second - FT(1)) < internal::tolerance<FT>());

      if (stats.first > FT(1) / FT(2))
        face.visibility = Visibility_label::INSIDE;
      else
        face.visibility = Visibility_label::OUTSIDE;
      face.inside = stats.first;
      face.outside = stats.second;
    }

    Stats estimate_in_out_values(const Face& polyhedron) {

      Point_3 b;
      internal::compute_barycenter_3(polyhedron.vertices, b);
      if (is_above_building(b))
        return std::make_pair(FT(0), FT(1));
      if (is_below_building(b))
        return std::make_pair(FT(0), FT(1));
      if (is_outside_boundary(b))
        return std::make_pair(FT(0), FT(1));
      return estimate_in_out_values_statistically(polyhedron);
    }

    bool is_above_building(const Point_3& query) const {
      return query.z() > m_building.top_z;
    }

    bool is_below_building(const Point_3& query) const {
      return query.z() < m_building.bottom_z;
    }

    bool is_outside_boundary(const Point_3& query) const {

      const Point_2 p = internal::point_2_from_point_3(query);
      const auto& tri = m_building.base1.triangulation.delaunay;
      const auto fh = tri.locate(p);
      return !fh->info().tagged;
    }

    Stats estimate_in_out_values_statistically(
      const Face& polyhedron) {

      std::size_t in = 0, out = 0;
      m_samples.clear();
      create_samples(polyhedron, m_samples);
      CGAL_assertion(m_samples.size() >= m_num_samples);

      compute_stats(m_samples, in, out);
      if (in == 0 && out == 0) in = 1;

      const FT tmp_in = static_cast<FT>(in);
      const FT tmp_out = static_cast<FT>(out);
      const FT sum = tmp_in + tmp_out;
      CGAL_assertion(sum > FT(0));

      const FT final_in = tmp_in  / sum;
      const FT final_out = tmp_out / sum;

      return std::make_pair(final_in, final_out);
    }

    void create_samples(
      const Face& polyhedron,
      std::vector<Point_3>& samples) const {

      const auto& vertices = polyhedron.vertices;
      Delaunay_3 delaunay_3;
      for (const auto& p : vertices)
        delaunay_3.insert(Local_point_3(
          CGAL::to_double(p.x()),
          CGAL::to_double(p.y()),
          CGAL::to_double(p.z())));

      std::vector<Local_point_3> points;
      for (auto cit = delaunay_3.finite_cells_begin();
      cit != delaunay_3.finite_cells_end(); ++cit) {
        const auto& tet = delaunay_3.tetrahedron(cit);
        Polyhedron tmp;
        tmp.make_tetrahedron(
          tet.vertex(0), tet.vertex(1), tet.vertex(2), tet.vertex(3));

        Generator generator(tmp);
        std::copy_n(generator, m_num_samples, std::back_inserter(points));
      }

      for (const auto& p : points)
        samples.push_back(Point_3(
          static_cast<FT>(p.x()),
          static_cast<FT>(p.y()),
          static_cast<FT>(p.z())));
    }

    void compute_stats(
      const std::vector<Point_3>& samples,
      std::size_t& in, std::size_t& out) const {

      std::vector<std::size_t> neighbors;
      for (const Point_3& p : samples) {
        const Point_2 q = internal::point_2_from_point_3(p);
        (*m_neighbor_query_ptr)(q, neighbors);
        CGAL_assertion(neighbors.size() == m_k);

        FT z = FT(0);
        for (const std::size_t idx : neighbors) {
          z += m_queries[idx].second;
        }
        z /= static_cast<FT>(neighbors.size());

        if (p.z() < z) in += 1;
        else out += 1;
      }
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_3_EXP_3_H
