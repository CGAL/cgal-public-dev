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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_FACE_SIMPLIFIER_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_FACE_SIMPLIFIER_H

// STL includes.
#include <map>
#include <set>
#include <list>
#include <queue>
#include <vector>
#include <utility>
#include <memory>

// CGAL includes.
#include <CGAL/enum.h>
#include <CGAL/property_map.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

// Other includes.
#include <CGAL/Levels_of_detail/internal/Reconstruction/Image.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"
#include "../../../../../test/Levels_of_detail/include/Utilities.h"

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

template<
typename GeomTraits,
typename Vertex,
typename Edge,
typename Halfedge,
typename Face>
struct Image_face_simplifier {

public:
  using Traits = GeomTraits;

  using FT = typename Traits::FT;
  using Point_2 = typename Traits::Point_2;
  using Point_3 = typename Traits::Point_3;
  using Segment_2 = typename Traits::Segment_2;
  using Segment_3 = typename Traits::Segment_3;
  using Line_2 = typename Traits::Line_2;
  using Line_3 = typename Traits::Line_3;
  using Vector_2 = typename Traits::Vector_2;
  using Plane_3 = typename Traits::Plane_3;
  using Triangle_2 = typename Traits::Triangle_2;
  using Triangle_3 = typename Traits::Triangle_3;
  using Intersect_2 = typename Traits::Intersect_2;
  using Intersect_3 = typename Traits::Intersect_3;

  using Image = internal::Image<Traits>;
  using Point_type = typename Image::Point_type;

  using Indices = std::vector<std::size_t>;
  using Size_pair = std::pair<std::size_t, std::size_t>;

  using Saver = Saver<Traits>;
  using Inserter = Polygon_inserter<Traits>;
  using Indexer = internal::Indexer<Point_3>;
  using Color = CGAL::Color;

  Image_face_simplifier(
    const std::vector<Segment_2>& boundary,
    std::vector<Vertex>& vertices,
    std::vector<Edge>& edges,
    std::vector<Halfedge>& halfedges,
    const Image& image,
    const std::map<std::size_t, Plane_3>& plane_map,
    const FT noise_level_2) :
  m_boundary(boundary),
  m_vertices(vertices),
  m_edges(edges),
  m_halfedges(halfedges),
  m_image(image),
  m_plane_map(plane_map),
  m_noise_level_2(noise_level_2),
  m_pi(static_cast<FT>(CGAL_PI))
  { }

  void simplify_face(const Face& face) {

    Indices free, linear;
    for (const std::size_t he_idx : face.hedges) {
      const auto& he = m_halfedges[he_idx];
      const std::size_t query_idx = he.from_vertex;
      auto& vertex = m_vertices[query_idx];

      if (vertex.type == Point_type::OUTER_CORNER) {
        simplify_both(free, linear); continue;
      }

      if (vertex.type == Point_type::OUTER_BOUNDARY) {
        simplify_both(free, linear);
        project_linear_vertex(vertex); continue;
      }

      if (vertex.type == Point_type::CORNER) {
        simplify_both(free, linear);
        project_linear_vertex(vertex); continue;
      }

      if (vertex.type == Point_type::FREE) {
        if (linear.size() != 0)
          simplify_linear(linear);
        free.push_back(query_idx);
      } else {
        if (free.size() != 0)
          simplify_free(free);
        linear.push_back(query_idx);
      }
    }
    simplify_both(free, linear);
  }

private:
  const std::vector<Segment_2>& m_boundary;
  std::vector<Vertex>& m_vertices;
  std::vector<Edge>& m_edges;
  std::vector<Halfedge>& m_halfedges;
  const Image& m_image;
  const std::map<std::size_t, Plane_3>& m_plane_map;
  const FT m_noise_level_2;
  const FT m_pi;

  void simplify_both(
    Indices& free, Indices& linear) {
    if (free.size() != 0) simplify_free(free);
    if (linear.size() != 0) simplify_linear(linear);
  }

  void simplify_linear(Indices& linear) {

    if (linear.size() == 1) {
      auto& vertex = m_vertices[linear[0]];
      project_linear_vertex(vertex);
      linear.clear();
      return;
    }

    /*
    std::cout << "linear: ";
    for (const std::size_t idx : linear) {
      std::cout << int(m_vertices[idx].type) << " ";
    }
    std::cout << std::endl; */

    /*
    for (const std::size_t idx : linear) {
      auto& vertex = m_vertices[idx];
      project_linear_vertex(vertex);
    } */

    /*
    std::cout.precision(30);
    for (const auto& idx : linear)
      std::cout << m_vertices[idx].point << std::endl;
    std::cout << std::endl; */

    simplify_linear_polyline(linear);
    linear.clear();
  }

  void project_linear_vertex(Vertex& vertex) {

    for (const std::size_t label : vertex.labels) {
      if (label == std::size_t(-1)) {
        /*
        vertex.used = true;
        if (vertex.labels.size() == 1)
          vertex.state = true; */
        return;
      }
    }

    if (vertex.labels.size() == 2) {
      project_onto_line(vertex); return;
    }
    if (vertex.labels.size() >= 3) {
      project_onto_point(vertex); return;
    }
  }

  void project_onto_line(Vertex& vertex) {

    if (vertex.used) return;

    auto it = vertex.labels.begin();
    const std::size_t l1 = *it; ++it;
    const std::size_t l2 = *it;

    const auto& plane1 = m_plane_map.at(l1);
    const auto& plane2 = m_plane_map.at(l2);

    Line_2 line;
    bool success = intersect_planes(plane1, plane2, line);
    if (!success) return;

    if (vertex.type == Point_type::LINEAR) {
      const auto proj = line.projection(vertex.point);
      /* vertex.point = proj; */
      vertex.used = true;
      return;
    }

    if (vertex.type == Point_type::OUTER_BOUNDARY) {
      if (vertex.bd_idx == std::size_t(-1)) return;

      const auto& segment = m_boundary[vertex.bd_idx];
      const auto other = Line_2(segment.source(), segment.target());

      Point_2 res;
      success = intersect_2(line, other, res);
      if (!success) return;

      const FT distance = internal::distance(res, vertex.point);
      if (distance < m_noise_level_2) {
        /* vertex.point = res; */
        vertex.used = true;
        return;
      }
    }
  }

  void project_onto_point(Vertex& vertex) {

    if (vertex.used) return;

    auto it = vertex.labels.begin();
    const std::size_t l1 = *it; ++it;
    const std::size_t l2 = *it; ++it;
    const std::size_t l3 = *it;

    const auto& plane1 = m_plane_map.at(l1);
    const auto& plane2 = m_plane_map.at(l2);
    const auto& plane3 = m_plane_map.at(l3);

    Point_2 res;
    const bool success = intersect_planes(plane1, plane2, plane3, res);
    if (!success) return;

    const FT distance = internal::distance(res, vertex.point);
    if (distance < m_noise_level_2) {

      if (vertex.type == Point_type::CORNER) {
        /* vertex.point = res; */
        vertex.used = true;
        return;
      }

      if (vertex.type == Point_type::OUTER_BOUNDARY) {
        if (vertex.bd_idx == std::size_t(-1)) return;

        const auto& segment = m_boundary[vertex.bd_idx];
        const auto line = Line_2(segment.source(), segment.target());
        const auto proj = line.projection(res);
        /* vertex.point = proj; */
        vertex.used = true;
        return;
      }
    }
  }

  bool intersect_planes(
    const Plane_3& plane1,
    const Plane_3& plane2,
    Line_2& line_2) {

    auto result = CGAL::intersection(plane1, plane2);
    Line_3 line_3; bool found = false;
    if (result) {
      if (const Line_3* l = boost::get<Line_3>(&*result)) {
        found = true; line_3 = *l;
      }
    }
    if (!found) return false;

    const auto p1 = line_3.point(0);
    const auto p2 = line_3.point(1);
    const auto q1 = Point_2(p1.x(), p1.y());
    const auto q2 = Point_2(p2.x(), p2.y());

    line_2 = Line_2(q1, q2);
    return true;
  }

  bool intersect_planes(
    const Plane_3& plane1,
    const Plane_3& plane2,
    const Plane_3& plane3,
    Point_2& point_2) {

    auto result = CGAL::intersection(plane1, plane2, plane3);
    Point_3 point_3; bool found = false;
    if (result) {
      if (const Point_3* p = boost::get<Point_3>(&*result)) {
        found = true; point_3 = *p;
      }
    }
    if (!found) return false;

    point_2 = Point_2(point_3.x(), point_3.y());
    return true;
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

  void simplify_linear_polyline(
    const Indices& polyline) {

    const std::size_t nump = polyline.size();
    auto& p = m_vertices[polyline[0]];
    auto& q = m_vertices[polyline[nump - 1]];

    project_linear_vertex(p);
    project_linear_vertex(q);

    for (std::size_t i = 1; i < nump - 1; ++i) {
      auto& vertex = m_vertices[polyline[i]];
      project_linear_vertex(vertex);
      vertex.state = true;
      vertex.used  = true;
    }
  }

  void simplify_free(Indices& free) {
    if (free.size() <= 2) {
      free.clear(); return;
    }

    /*
    std::cout << "free: ";
    for (const std::size_t idx : free) {
      std::cout << int(m_vertices[idx].type) << " ";
    }
    std::cout << std::endl; */

    simplify_free_polyline(free);
    free.clear();
  }

  void simplify_free_polyline(
    const Indices& polyline) {

    /*
    bool found = false;
    for (const std::size_t idx : polyline) {
      auto& vertex = m_vertices[idx];

      for (const std::size_t label : vertex.labels) {
        if (label == std::size_t(-1)) {
          vertex.used = true;
          if (vertex.labels.size() == 1)
            vertex.state = true;
          found = true; break;
        }
      }
    }
    if (found) return; */

    std::vector<Point_2> points;
    points.reserve(polyline.size());

    for (const std::size_t idx : polyline)
      points.push_back(m_vertices[idx].point);

    using Cost = CGAL::Polyline_simplification_2::Squared_distance_cost;
    using Stop = CGAL::Polyline_simplification_2::Stop_above_cost_threshold;

    const double threshold = m_noise_level_2 / FT(10);

    Cost cost;
    Stop stop(threshold);
    std::vector<Point_2> result;
    CGAL::Polyline_simplification_2::simplify(
      points.begin(), points.end(), cost, stop,
      std::back_inserter(result));

    for (const std::size_t idx : polyline) {
      auto& vertex = m_vertices[idx];
      vertex.used = true;

      if (vertex.type == Point_type::OUTER_CORNER) continue;
      if (vertex.type == Point_type::OUTER_BOUNDARY) continue;
      if (vertex.type == Point_type::BOUNDARY) continue;

      if (is_simplified(vertex.point, result))
        vertex.state = true;
    }
  }

  bool is_simplified(
    const Point_2& query,
    const std::vector<Point_2>& points) {

    for (const auto& point : points) {
      if (internal::are_equal_points_2(query, point)) {
        /* std::cout << "found point" << std::endl; */
        return false;
      }
    }
    return true;
  }
};

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_IMAGE_FACE_SIMPLIFIER_H
