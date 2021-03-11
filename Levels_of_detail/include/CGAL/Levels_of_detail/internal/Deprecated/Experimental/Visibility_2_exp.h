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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_2_EXP_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_2_EXP_H

// Example:

// void compute_visibility_2() {

//   if (!m_boundaries_detected) return;
//   if (m_partition_2.empty()) return;

//   Visibility_2_exp visibility(
//     m_data, m_boundary_points, m_interior_points, m_partition_2);
//   visibility.compute();
//   return;
// }

// STL includes.
#include <vector>
#include <utility>
#include <memory>

// CGAL includes.
#include <CGAL/assertions.h>
#include <CGAL/barycenter.h>
#include <CGAL/property_map.h>

// LOD includes.
#include <CGAL/Levels_of_detail/enum.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename DataStructure>
  class Visibility_2_exp {

  public:
    using Data_structure = DataStructure;
    using Traits = typename Data_structure::Traits;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Vector_2 = typename Traits::Vector_2;
    using Triangle_2 = typename Traits::Triangle_2;

    using Partition_2 = internal::Partition_2<Traits>;
    using Face = typename Partition_2::Face;

    using K_neighbor_query =
    internal::K_neighbor_query<Traits, std::vector< std::pair<Point_2, double> >,
    CGAL::First_of_pair_property_map< std::pair<Point_2, double> > >;
    using Location_type = typename Triangulation<Traits>::Delaunay::Locate_type;

    Visibility_2_exp(
      const Data_structure& data,
      const std::vector<std::size_t>& boundary_points,
      const std::vector<std::size_t>& interior_points,
      Partition_2 &partition_2) :
    m_data(data),
    m_boundary_points(boundary_points),
    m_interior_points(interior_points),
    m_partition_2(partition_2),
    m_num_probes(100)
    { }

    void compute() {

      // Buildings.
      std::vector< std::pair<std::size_t, bool> > items;
      items.reserve(m_interior_points.size() + m_boundary_points.size());
      for (const auto& idx : m_interior_points)
        items.push_back(std::make_pair(idx, false));
      for (const auto& idx : m_boundary_points)
        items.push_back(std::make_pair(idx, false));

      // Ground and vegetation.
      std::vector<std::size_t> gr, veg;
      m_data.points(Semantic_label::GROUND, gr);
      m_data.points(Semantic_label::VEGETATION, veg);

      m_input_range.clear();
      for (const auto& face : m_partition_2.faces) {
        const auto& tri = face.base.delaunay;

        // Ground.
        for (const auto& idx : gr) {
          const auto& p = get(m_data.point_map_2, idx);
          Location_type type; int stub;
          const auto fh = tri.locate(p, type, stub);
          if (type == Triangulation<Traits>::Delaunay::FACE &&
          !tri.is_infinite(fh)) {
            m_input_range.push_back(std::make_pair(
            p, get(m_data.visibility_map_d, idx)));
          }
        }

        // Vegetation.
        for (const auto& idx : veg) {
          const auto& p = get(m_data.point_map_2, idx);
          Location_type type; int stub;
          const auto fh = tri.locate(p, type, stub);
          if (type == Triangulation<Traits>::Delaunay::FACE &&
          !tri.is_infinite(fh)) {
            m_input_range.push_back(std::make_pair(
            p, get(m_data.visibility_map_d, idx)));
          }
        }

        // Buildings.
        for (auto& item : items) {
          if (!item.second) {

            const auto& p = get(m_data.point_map_2, item.first);
            Location_type type; int stub;
            const auto fh = tri.locate(p, type, stub);
            if (type == Triangulation<Traits>::Delaunay::FACE &&
            !tri.is_infinite(fh)) {
                m_input_range.push_back(std::make_pair(
                p, get(m_data.visibility_map_d, item.first)));
              item.second = true;
            }
          }
        }
      }

      // Neighbor query.
      m_neighbor_query = std::make_shared<K_neighbor_query>(
        m_input_range, FT(1), m_point_map_2);
      estimate_visibility();
    }

  private:
    const Data_structure& m_data;
    const std::vector<std::size_t>& m_boundary_points;
    const std::vector<std::size_t>& m_interior_points;
    const std::size_t m_num_probes;

    Partition_2& m_partition_2;
    CGAL::First_of_pair_property_map< std::pair<Point_2, double> > m_point_map_2;
    CGAL::Second_of_pair_property_map< std::pair<Point_2, double> > m_visibility_map;
    std::shared_ptr<K_neighbor_query> m_neighbor_query;
    std::vector< std::pair<Point_2, double> > m_input_range;

    void estimate_visibility() const {

      if (m_partition_2.empty()) return;
      label_exterior_faces(m_partition_2.faces);
      for (auto& face : m_partition_2.faces)
        if (face.exterior) {
          face.visibility = Visibility_label::OUTSIDE;
          face.inside = FT(0);
          face.outside = FT(1);
        }
        else compute_monte_carlo_label(face);
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
        if (face.exterior)
          continue;
      }
    }

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
        (*m_neighbor_query)(point, neighbors);
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
          internal::random_point_in_triangle_2(probability[i].first, point);
          return;
        }
      }
      std::cerr <<
        "Error (compute_random_point_in_triangles()): probability is out of range!"
      << std::endl;
      point = Point_2(FT(0), FT(0));
    }

    FT get_function_value(
      const Point_2& p,
      const std::vector<std::size_t>& neighbors) const {

      const double value =
      get(m_visibility_map, *(m_input_range.begin() + neighbors[0]));
      return static_cast<FT>(value);
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_VISIBILITY_2_EXP_H
