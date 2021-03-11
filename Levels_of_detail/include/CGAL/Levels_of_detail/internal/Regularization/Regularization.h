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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_REGULARIZATION_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_REGULARIZATION_H

#include <CGAL/license/Levels_of_detail.h>

#include <CGAL/Levels_of_detail/internal/Regularization/Shape_regularization.h>
#include <CGAL/Levels_of_detail/internal/Regularization/Angle_regularization_2.h>
#include <CGAL/Levels_of_detail/internal/Regularization/Ordinate_regularization_2.h>
#include <CGAL/Levels_of_detail/internal/Regularization/Delaunay_neighbor_query_2.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Regularization {

  public:
    using Traits = GeomTraits;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Segment_2 = typename Traits::Segment_2;

    using Indices = std::vector<std::size_t>;
    using Segment_range = std::vector<Segment_2>;
    using Segment_map = CGAL::Identity_property_map<Segment_2>;

    using Neighbor_query =
    internal::Delaunay_neighbor_query_2<Traits, Segment_range, Segment_map>;
    using RT_angles =
    internal::Angle_regularization_2<Traits, Segment_range, Segment_map>;
    using RT_ordinates =
    internal::Ordinate_regularization_2<Traits, Segment_range, Segment_map>;

    using SR_angles = internal::Shape_regularization
    <Traits, Segment_range, Neighbor_query, RT_angles>;
    using SR_ordinates = internal::Shape_regularization
    <Traits, Segment_range, Neighbor_query, RT_ordinates>;

    using Saver = Saver<Traits>;
    using Color = CGAL::Color;

    struct Polygon_neighbor_query {

      Polygon_neighbor_query(
        const Segment_range& segments) :
      m_segments(segments) { }

      void operator()(
        const std::size_t query_index,
        std::vector<std::size_t>& neighbors) const {

        CGAL_assertion(query_index >= 0);
        CGAL_assertion(query_index < m_segments.size());

        neighbors.clear();
        const std::size_t n = m_segments.size();

        const std::size_t im = (query_index + n - 1) % n;
        const std::size_t ip = (query_index + 1) % n;
        neighbors.push_back(im);
        neighbors.push_back(ip);
      }

    private:
      const Segment_range& m_segments;
    };

    void regularize_polygon_angles(
      std::vector<Segment_range>& contours,
      const FT angle_bound) {

      using PR_angles = internal::Shape_regularization
      <Traits, Segment_range, Polygon_neighbor_query, RT_angles>;

      m_parallel_contours.clear();
      std::vector<std::size_t> group;
      std::vector<Indices> parallel_groups;
      for (std::size_t i = 0; i < contours.size(); ++i) {

        auto& contour = contours[i];
        if (contour.size() <= 2) continue;
        Polygon_neighbor_query neighbor_query(contour);

        group.clear();
        group.resize(contour.size());
        std::iota(group.begin(), group.end(), 0);

        RT_angles rt_angles(contour, angle_bound);
        rt_angles.add_group(group);
        rt_angles.make_bounds();

        PR_angles pr_angles(contour, neighbor_query, rt_angles);
        pr_angles.regularize();

        parallel_groups.clear();
        rt_angles.parallel_groups(std::back_inserter(parallel_groups));
        m_parallel_contours.push_back(parallel_groups);
      }
    }

    void regularize_polygon_ordinates(
      std::vector<Segment_range>& contours,
      const FT ordinate_bound) {

      for (std::size_t i = 0; i < contours.size(); ++i) {
        auto& contour = contours[i];
        const auto& parallel_groups = m_parallel_contours[i];

        if (contour.size() <= 2 || parallel_groups.size() == 0)
          continue;
        std::size_t num_segments = 0;
        for (const auto& parallel_group : parallel_groups)
          num_segments += parallel_group.size();
        if (num_segments <= 2) return;

        RT_ordinates rt_ordinates(contour, ordinate_bound);
        Neighbor_query neighbor_query(contour);
        for(const auto& parallel_group : parallel_groups) {
          neighbor_query.add_group(parallel_group);
          rt_ordinates.add_group(parallel_group);
        }

        SR_ordinates sr_ordinates(contour, neighbor_query, rt_ordinates);
        sr_ordinates.regularize();
      }
    }

    void regularize_angles(
      Segment_range& segments,
      const FT angle_bound) {
      std::cout << "Angle bound: " << angle_bound << std::endl;

      if (segments.size() <= 2)
        return;

      Neighbor_query neighbor_query(segments);

      std::vector<std::size_t> group;
      group.resize(segments.size());
      std::iota(group.begin(), group.end(), 0);

      neighbor_query.add_group(group);

      RT_angles rt_angles(segments, angle_bound);
      rt_angles.add_group(group);

      rt_angles.make_bounds();
      SR_angles sr_angles(segments, neighbor_query, rt_angles);
      sr_angles.regularize();

      m_parallel_groups.clear();
      rt_angles.parallel_groups(
        std::back_inserter(m_parallel_groups));
    }

    void regularize_ordinates(
      Segment_range& segments,
      const FT ordinate_bound) {
      std::cout << "Ordinate bound: " << ordinate_bound << std::endl;

      if (segments.size() <= 2 || m_parallel_groups.size() == 0)
        return;
      std::size_t num_segments = 0;
      for (const auto& parallel_group : m_parallel_groups)
        num_segments += parallel_group.size();
      if (num_segments <= 2) return;

      CGAL_assertion(m_parallel_groups.size() > 0);
      RT_ordinates rt_ordinates(segments, ordinate_bound);

      Neighbor_query neighbor_query(segments);
      for(const auto& parallel_group : m_parallel_groups) {
        neighbor_query.add_group(parallel_group);
        rt_ordinates.add_group(parallel_group);
      }

      rt_ordinates.make_bounds();
      SR_ordinates sr_ordinates(segments, neighbor_query, rt_ordinates);
      sr_ordinates.regularize();
    }

    void save_polylines(
      const Segment_range& segments,
      const std::string name) {

      CGAL_assertion(segments.size() > 0);
      std::vector< std::vector<Point_3> > polylines(segments.size());
      for (std::size_t i = 0; i < segments.size(); ++i) {
        const Point_2& s = segments[i].source();
        const Point_2& t = segments[i].target();

        polylines[i].push_back(Point_3(s.x(), s.y(), FT(0)));
        polylines[i].push_back(Point_3(t.x(), t.y(), FT(0)));
      }

      Saver saver;
      saver.export_polylines(polylines, name);
    }

  private:
    std::vector<Indices> m_parallel_groups;
    std::vector< std::vector<Indices> > m_parallel_contours;
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_REGULARIZATION_H