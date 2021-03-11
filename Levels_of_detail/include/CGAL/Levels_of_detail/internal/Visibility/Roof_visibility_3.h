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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_ROOF_VISIBILITY_3_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_ROOF_VISIBILITY_3_H

// STL includes.
#include <vector>
#include <utility>
#include <iostream>
#include <fstream>
#include <memory>
#include <algorithm>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>

// Internal includes.
#include <CGAL/Levels_of_detail/enum.h>
#include <CGAL/Levels_of_detail/internal/struct.h>
#include <CGAL/Levels_of_detail/internal/utils.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Roof_visibility_3 {

  public:
    using Traits = GeomTraits;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Plane_3 = typename Traits::Plane_3;

    using Indices = std::vector<std::size_t>;
    using Partition_2 = internal::Partition_2<Traits>;
    using Partition_3 = internal::Partition_3<Traits>;
    using Face_2 = typename Partition_2::Face;
    using Face_3 = typename Partition_3::Face;
    using Stats = std::pair<FT, FT>;

    using Local_traits = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Local_point_3 = typename Local_traits::Point_3;
    using Delaunay_3 = CGAL::Delaunay_triangulation_3<Local_traits>;
    using Generator = CGAL::Random_points_in_tetrahedron_3<Local_point_3>;

    using Delaunay = typename Triangulation<Traits>::Delaunay;
    using Location_type = typename Delaunay::Locate_type;

    Roof_visibility_3(
      const Partition_2& partition_2,
      const FT bottom_z, const FT top_z) :
    m_partition_2(partition_2),
    m_bottom_z(bottom_z),
    m_top_z(top_z),
    m_num_samples(100),
    m_random(0)
    { }

    void compute(Partition_3& partition_3) {

      if (partition_3.empty()) return;
      label_exterior_faces(partition_3.faces);
      for (auto& pface : partition_3.faces) {
        if (pface.exterior) {

          pface.visibility = Visibility_label::OUTSIDE;
          pface.inside  = FT(0);
          pface.outside = FT(1);
        } else compute_face_label(pface);
      }
    }

  private:
    const Partition_2& m_partition_2;
    const FT m_bottom_z;
    const FT m_top_z;
    const std::size_t m_num_samples;
    Random m_random;

    void label_exterior_faces(
      std::vector<Face_3>& faces) const {

      for (auto& face : faces) {
        face.exterior = false;
        const auto& neighbors = face.neighbors;
        for (const int idx : neighbors) {
          if (idx < 0) {
            face.exterior = true; break;
          }
        }
      }
    }

    void compute_face_label(Face_3& face) {

      const Stats stats = estimate_in_out_values(face);
      CGAL_assertion(stats.first  >= FT(0) && stats.first  <= FT(1));
      CGAL_assertion(stats.second >= FT(0) && stats.second <= FT(1));
      CGAL_assertion(
        CGAL::abs(stats.first + stats.second - FT(1)) < internal::tolerance<FT>());

      if (stats.first >= FT(1) / FT(2))
        face.visibility = Visibility_label::INSIDE;
      else
        face.visibility = Visibility_label::OUTSIDE;
      face.inside  = stats.first;
      face.outside = stats.second;
    }

    Stats estimate_in_out_values(const Face_3& polyhedron) {

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
      return query.z() > m_top_z;
    }

    bool is_below_building(const Point_3& query) const {
      return query.z() < m_bottom_z;
    }

    bool is_outside_boundary(const Point_3& query) const {

      const Point_2 p = internal::point_2_from_point_3(query);
      for (const auto& pface : m_partition_2.faces) {
        if (pface.visibility == Visibility_label::INSIDE)
          continue;

        const auto& base = pface.base;
        Location_type type; int stub;
        const auto fh = base.delaunay.locate(p, type, stub);
        if (
          type != Delaunay::OUTSIDE_CONVEX_HULL &&
          type != Delaunay::OUTSIDE_AFFINE_HULL)
          return true;
      }
      return false;
    }

    Stats estimate_in_out_values_statistically(
      const Face_3& polyhedron) {

      std::size_t in = 0, out = 0;
      std::vector<Point_3> samples;
      create_samples(polyhedron, samples);
      compute_stats(samples, in, out);
      if (in == 0 && out == 0) {
        in = 1; out = 1;
      }

      const FT tmp_in = static_cast<FT>(in);
      const FT tmp_out = static_cast<FT>(out);
      const FT sum = tmp_in + tmp_out;
      CGAL_assertion(sum > FT(0));

      const FT final_in  = tmp_in  / sum;
      const FT final_out = tmp_out / sum;

      return std::make_pair(final_in, final_out);
    }

    void create_samples(
      const Face_3& polyhedron,
      std::vector<Point_3>& samples) {

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

        const FT volume = CGAL::abs(tet.volume());
        if (volume < FT(1) / FT(1000))
          continue;

        Generator generator(tet, m_random);
        std::copy_n(generator, m_num_samples, std::back_inserter(points));
      }

      samples.clear();
      samples.reserve(points.size());

      for (const auto& p : points)
        samples.push_back(Point_3(
          static_cast<FT>(p.x()),
          static_cast<FT>(p.y()),
          static_cast<FT>(p.z())));
    }

    void compute_stats(
      const std::vector<Point_3>& samples,
      std::size_t& in, std::size_t& out) const {

      for (const auto& p : samples)
        handle_sample_point(p, in, out);
    }

    void handle_sample_point(
      const Point_3& query,
      std::size_t& in, std::size_t& out) const {

      const Point_2 p = internal::point_2_from_point_3(query);
      for (const auto& pface : m_partition_2.faces) {
        if (pface.visibility == Visibility_label::OUTSIDE)
          continue;

        const auto& base = pface.base;
        Location_type type; int stub;
        const auto fh = base.delaunay.locate(p, type, stub);
        if (
          type != Delaunay::OUTSIDE_CONVEX_HULL &&
          type != Delaunay::OUTSIDE_AFFINE_HULL) {

          const auto& plane = pface.plane;
          const Point_3 ref = internal::position_on_plane_3(p, plane);

          if (query.z() < ref.z()) in += 1;
          else out += 1;
          return;
        }
      }
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_ROOF_VISIBILITY_3_H
