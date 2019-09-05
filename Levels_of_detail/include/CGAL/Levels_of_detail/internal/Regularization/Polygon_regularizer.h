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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_POLYGON_REGULARIZER_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_POLYGON_REGULARIZER_H

#include <CGAL/license/Levels_of_detail.h>

// STL includes.
#include <vector>
#include <utility>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename GeomTraits>
  class Polygon_regularizer {

  public:
    using Traits = GeomTraits;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Segment_2 = typename Traits::Segment_2;
    using Vector_2 = typename Traits::Vector_2;
    using Line_2 = typename Traits::Line_2;
    using Intersect_2 = typename Traits::Intersect_2;

    using Saver = Saver<Traits>;
    using Color = CGAL::Color;

    Polygon_regularizer(
      const FT min_length,
      const FT angle_bound) :
    m_min_length(min_length),
    m_angle_bound(angle_bound),
    m_pi(static_cast<FT>(CGAL_PI)),
    m_bound_1(   m_angle_bound ),
    m_bound_2( ( m_pi / FT(2)  ) * FT(180) / m_pi - m_bound_1 ),
    m_bound_3(   FT(10) ) { 

      std::cout << "bound 1: " << m_bound_1 << std::endl;
      std::cout << "bound 2: " << m_bound_2 << std::endl;
      std::cout << "bound 3: " << m_bound_3 << std::endl;
    }

    void compute_longest_direction(
      const std::vector< std::vector<Segment_2> >& contours) {

      const auto pair = find_longest_segment(contours);
      m_ref_con = pair.first;
      m_ref_idx = pair.second;
      m_longest = contours[m_ref_con][m_ref_idx];
    }

    void regularize_contours(
      std::vector< std::vector<Segment_2> >& contours) {
      
      for (std::size_t k = 0; k < contours.size(); ++k) {
        auto& contour = contours[k];
        if (contour.size() < 4) 
          continue;
        
        rotate_contour(k, contour);
        correct_contour(k, contour);
        connect_contour(contour);
      }
    }

    void rotate_contour(
      const std::size_t con_idx,
      std::vector<Segment_2>& contour) {

      const std::size_t n = contour.size();
      for (std::size_t i = 0; i < n; ++i) {
        if (con_idx == m_ref_con && i == m_ref_idx) 
          continue;

        auto& si = contour[i];
        rotate_segment(m_longest, si);
      }
    }

    void correct_contour(
      const std::size_t con_idx,
      std::vector<Segment_2>& contour) {

      const std::size_t n = contour.size();
      for (std::size_t i = 0; i < n; ++i) {
        if (con_idx == m_ref_con && i == m_ref_idx) 
          continue;

        const std::size_t im = (i + n - 1) % n;
        const std::size_t ip = (i + 1) % n;

        auto& si = contour[i];
        const auto& sm = contour[im];
        const auto& sp = contour[ip];

        correct_segment(sm, si, sp);
      }
    }

    void connect_contour(
      std::vector<Segment_2>& contour) {
      
      std::vector<Segment_2> segments;
      check_segments(contour, segments);
      intersect_segments(segments);
      contour = segments;
    }

    void check_segments(
      const std::vector<Segment_2>& contour,
      std::vector<Segment_2>& segments) {
      
      segments.clear();
      const std::size_t n = contour.size();
      for (std::size_t i = 0; i < n; ++i) {
        
        const std::size_t im = (i + n - 1) % n;
        const std::size_t ip = (i + 1) % n;
        
        const auto& si = contour[i];
        const auto& sm = contour[im];
        const auto& sp = contour[ip];
        
        const bool success = check_segment(sm, si, sp);
        if (success)
          segments.push_back(si);
      }
    }

    void intersect_segments(
      std::vector<Segment_2>& segments) {

      const std::size_t n = segments.size();
      for (std::size_t i = 0; i < n; ++i) {
        
        const std::size_t im = (i + n - 1) % n;
        const std::size_t ip = (i + 1) % n;
        
        auto& si = segments[i];
        const auto& sm = segments[im];
        const auto& sp = segments[ip];
        
        intersect_segment(sm, si, sp);
      }
    }

  private:
    const FT m_min_length;
    const FT m_angle_bound;

    const FT m_pi;
    const FT m_bound_1, m_bound_2, m_bound_3;

    std::size_t m_ref_con;
    std::size_t m_ref_idx;

    Segment_2 m_longest;

    std::pair<std::size_t, std::size_t> find_longest_segment(
      const std::vector< std::vector<Segment_2> >& contours) {

      std::size_t con_idx = std::size_t(-1);
      std::size_t idx = std::size_t(-1);

      FT max_length = -FT(1);
      for (std::size_t k = 0; k < contours.size(); ++k) {
        for (std::size_t i = 0; i < contours[k].size(); ++i) {
          
          const FT length = contours[k][i].squared_length();
          if (length > max_length) {

            max_length = length;
            con_idx = k;
            idx = i;
          }
        }
      }
      return std::make_pair(con_idx, idx);
    }

    void rotate_segment(
      const Segment_2& longest,
      Segment_2& si) {

      const FT angle = angle_degree_2(longest, si);
      const FT angle_2 = get_angle_2(angle);

      if (CGAL::abs(angle_2) <= m_bound_1) {
        rotate(angle, FT(180), longest, si); // parallel
        return;
      }
      if (CGAL::abs(angle_2) >= m_bound_2) {
        rotate(angle, FT(90), longest, si); // orthogonal
        return;
      }
    }

    void correct_segment(
      const Segment_2& sm, Segment_2& si, const Segment_2& sp) {

      const FT angle_mp = angle_degree_2(sm, sp);
      const FT angle_mp_2 = get_angle_2(angle_mp);

      if (CGAL::abs(angle_mp_2) <= m_bound_3) {
        const FT angle = angle_degree_2(sm, si);
        rotate(angle, FT(90), sm, si); // orthogonal
        return;
      }
    }

    bool check_segment(
      const Segment_2& sm, const Segment_2& si, const Segment_2& sp) {

      // Source.
      const FT angle_mi = angle_degree_2(sm, si);
      const FT angle_mi_2 = get_angle_2(angle_mi);

      if (CGAL::abs(angle_mi_2) <= m_bound_3) {
        if (si.squared_length() > sm.squared_length())
          return true;
        else
          return false;
      }

      // Target.
      const FT angle_pi = angle_degree_2(si, sp);
      const FT angle_pi_2 = get_angle_2(angle_pi);

      if (CGAL::abs(angle_pi_2) <= m_bound_3) {
        if (si.squared_length() > sp.squared_length())
          return true;
        else
          return false;
      }

      return true;
    }

    void intersect_segment(
      const Segment_2& sm, Segment_2& si, const Segment_2& sp) {

      Point_2 source = si.source();
      Point_2 target = si.target();

      const Line_2 line_1 = Line_2(sm.source(), sm.target());
      const Line_2 line_2 = Line_2(si.source(), si.target());
      const Line_2 line_3 = Line_2(sp.source(), sp.target());

      // const bool success1 = intersect_2(line_1, line_2, source);
      // const bool success2 = intersect_2(line_2, line_3, target);

      // if (!success1) source = si.source();
      // if (!success2) target = si.target();

      source = line_1.projection(si.source());
      target = line_3.projection(si.target());

      si = Segment_2(source, target);
    }

    bool intersect_2(
      const Line_2& line_1, const Line_2& line_2,
      Point_2& in_point) {
      
      typename std::result_of<Intersect_2(Line_2, Line_2)>::type result 
      = CGAL::intersection(line_1, line_2);
      if (result) {
        if (const Line_2* line = boost::get<Line_2>(&*result)) 
          return false;
        else {
          const Point_2* point = boost::get<Point_2>(&*result);
          in_point = *point; return true;
        }
      }
      return false;
    }

    FT angle_degree_2(
      const Segment_2& longest, const Segment_2& si) {

      const Vector_2 v1 = si.to_vector();
      const Vector_2 v2 = -longest.to_vector();

		  const FT det = CGAL::determinant(v1, v2);
		  const FT dot = CGAL::scalar_product(v1, v2);
      const FT angle_rad = static_cast<FT>(
        std::atan2(CGAL::to_double(det), CGAL::to_double(dot)));
      const FT angle_deg = angle_rad * FT(180) / m_pi; 
      return angle_deg;
    }

    FT get_angle_2(const FT angle) {
      
      FT angle_2 = angle;
      if (angle_2 > FT(90)) angle_2 = FT(180) - angle_2;
      else if (angle_2 < -FT(90)) angle_2 = FT(180) + angle_2;
      return angle_2;
    }

    void rotate(
      const FT angle_2, 
      const FT ref_angle_2,
      const Segment_2& longest, 
      Segment_2& si) {

      FT angle = angle_2;
      if (angle < FT(0)) angle = angle + ref_angle_2;
      else if (angle > FT(0)) angle = angle - ref_angle_2;

      Point_2 source_i = si.source();
      Point_2 target_i = si.target();
      const Point_2 b = internal::middle_point_2(source_i, target_i);
      const FT angle_rad = angle * m_pi / FT(180); 
      internal::rotate_point_2(angle_rad, b, source_i);
      internal::rotate_point_2(angle_rad, b, target_i);
      si = Segment_2(source_i, target_i);
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_POLYGON_REGULARIZER_H
