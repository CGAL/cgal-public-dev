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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_BUILDINGS_SITE_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_BUILDINGS_SITE_H

#include <CGAL/license/Levels_of_detail.h>

// STL includes.
#include <memory>
#include <vector>
#include <utility>

// Boost includes.
#include <boost/optional/optional.hpp>

// CGAL includes.
#include <CGAL/property_map.h>

// LOD includes.
#include <CGAL/Levels_of_detail/enum.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utils.h>
#include <CGAL/Levels_of_detail/internal/struct.h>
#include <CGAL/Levels_of_detail/internal/number_utils.h>

// Simplification.
#include <CGAL/Levels_of_detail/internal/Simplification/Thinning_2.h>
#include <CGAL/Levels_of_detail/internal/Simplification/Grid_based_filtering_2.h>
#include <CGAL/Levels_of_detail/internal/Simplification/Alpha_shapes_filtering_2.h>

// Spacial search.
#include <CGAL/Levels_of_detail/internal/Spacial_search/K_neighbor_query.h>
#include <CGAL/Levels_of_detail/internal/Spacial_search/Sphere_neighbor_query.h>

// Shape detection.
#include <CGAL/Levels_of_detail/internal/Shape_detection/Region_growing.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Estimate_normals_2.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Least_squares_line_fit_region.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Least_squares_line_fit_sorting.h>

// Buildings.
#include <CGAL/Levels_of_detail/internal/Buildings/Building_builder.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename DataStructure>
  class Buildings_site {

  public:
    using Data_structure = DataStructure;

    using Traits = typename Data_structure::Traits;
    using Point_map = typename Data_structure::Point_map;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Vector_2 = typename Traits::Vector_2;
    using Vector_3 = typename Traits::Vector_3;
    using Plane_3 = typename Traits::Plane_3;
    using Segment_2 = typename Traits::Segment_2;
    using Segment_3 = typename Traits::Segment_3;
    using Line_2 = typename Traits::Line_2;

    using Points = std::vector<std::size_t>;
    using Point_map_3 = typename Data_structure::Point_map_3;

    using Building = internal::Building<Traits>;
    using Building_ptr = std::shared_ptr<Building>;
    using Building_builder = internal::Building_builder<Traits, Points, Point_map_3>;
    
    using Indexer = internal::Indexer<Point_3>;
    
    using Points_2 = std::vector<Point_2>;
    using Vectors_2 = std::vector<Vector_2>;
    using Pair_item_2 = std::pair<Point_2, Vector_2>;
    using Pair_range_2 = std::vector<Pair_item_2>;
    using Identity_map = CGAL::Identity_property_map<Point_2>;
    using First_of_pair_map = CGAL::First_of_pair_property_map<Pair_item_2>;
    using Second_of_pair_map = CGAL::Second_of_pair_property_map<Pair_item_2>;
    using Boundary_point_map_2 = 
    internal::Item_property_map<Points_2, Identity_map>;
    
    using K_neighbor_query =
    internal::K_neighbor_query<Traits, Points_2, Identity_map>;
    using Sphere_neighbor_query =
    internal::Sphere_neighbor_query<Traits, Points_2, Identity_map>;

    using Grid_based_filtering_2 = internal::Grid_based_filtering_2<Traits>;
    using Alpha_shapes_filtering_2 = internal::Alpha_shapes_filtering_2<Traits>;
    using Thinning_2 = internal::Thinning_2<Traits, Sphere_neighbor_query>;

    using Normal_estimator_2 = 
    internal::Estimate_normals_2<Traits, Points_2, Identity_map, K_neighbor_query>;
    using LSLF_region = 
    internal::Least_squares_line_fit_region<Traits, Pair_range_2, First_of_pair_map, Second_of_pair_map>;
    using LSLF_sorting =
    internal::Least_squares_line_fit_sorting<Traits, Points_2, K_neighbor_query, Identity_map>;
    using Region_growing_2 = 
    internal::Region_growing<Points_2, K_neighbor_query, LSLF_region, typename LSLF_sorting::Seed_map>;

    Buildings_site(
      const Data_structure& data,
      const Points& interior_points,
      const Points& boundary_points,
      const std::size_t site_index) : 
    m_data(data),
    m_interior_points(interior_points),
    m_boundary_points(boundary_points),
    m_site_index(site_index) { 
      CGAL_precondition(m_interior_points.size() > 0);
      CGAL_precondition(m_boundary_points.size() >= 0);
      create_ground_plane();
    }

    void detect_boundaries() {

      extract_boundary_points_2(
        m_data.parameters.buildings.alpha_shape_size_2, 
        m_data.parameters.buildings.grid_cell_width_2);

      extract_wall_points_2(
        m_data.parameters.buildings.region_growing_scale_2,
        m_data.parameters.buildings.region_growing_noise_level_2,
        m_data.parameters.buildings.region_growing_angle_2,
        m_data.parameters.buildings.region_growing_min_length_2);

      compute_approximate_boundaries();
    }

    void compute_footprints() {
      
    }

    void extrude_footprints() {
      
    }

    void detect_roofs() {

    }

    void compute_roofs() {

    }

    void get_buildings(
      std::vector<Building_ptr>& buildings) const {
      
      if (m_buildings.empty()) return;
      buildings.clear();
      for (const auto& building : m_buildings)
        buildings.push_back(std::make_shared<Building>(building));
    }

    template<typename OutputIterator>
    boost::optional<OutputIterator> 
    get_boundary_points(
      OutputIterator output) const {

      if (m_boundary_points_2.empty())
        return boost::none;
      for (const Point_2& p : m_boundary_points_2) {
        const Point_3 q = internal::position_on_plane_3(p, m_ground_plane);
        *(output++) = std::make_pair(q, std::size_t(-1));
      }
      return output;
    }

    template<typename OutputIterator>
    boost::optional<OutputIterator> 
    get_wall_points(
      OutputIterator output) const {

      if (m_wall_points_2.empty())
        return boost::none;
      
      Identity_map identity_map;
      Boundary_point_map_2 point_map_2(m_boundary_points_2, identity_map);

      for (std::size_t i = 0; i < m_wall_points_2.size(); ++i) {
        for (const std::size_t idx : m_wall_points_2[i]) {
          const Point_2& p = get(point_map_2, idx);
          const Point_3 q = internal::position_on_plane_3(p, m_ground_plane);
          *(output++) = std::make_pair(q, i);
        }
      }
      return output;
    }

    template<typename OutputIterator>
    boost::optional<OutputIterator> 
    get_approximate_boundaries(
      OutputIterator output) const {
      
      if (m_approximate_boundaries_2.empty())
        return boost::none;
        
      for (std::size_t i = 0; i < m_approximate_boundaries_2.size(); ++i) {
        const Point_2& s = m_approximate_boundaries_2[i].source();
        const Point_2& t = m_approximate_boundaries_2[i].target();

        const Point_3 p = internal::position_on_plane_3(s, m_ground_plane);
        const Point_3 q = internal::position_on_plane_3(t, m_ground_plane);
        *(output++) = std::make_pair(Segment_3(p, q), i);
      }
      return output;
    }

    /*
    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> > 
    get_building_footprints(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      std::size_t& building_index) const {
      
      for (const auto& building : m_buildings) {
        building.--.output_for_object(
          indexer, num_vertices, vertices, faces, building_index);
        ++building_index;
      }
      return std::make_pair(vertices, faces);
    } 

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> > 
    get_extruded_building_boundaries(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      std::size_t& building_index) const {
      
      for (const auto& building : m_buildings) {
        building.--.output_for_object(
          indexer, num_vertices, vertices, faces, building_index);
        ++building_index;
      }
      return std::make_pair(vertices, faces);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> > 
    get_extruded_building_footprints(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      std::size_t& building_index) const {
      
      for (const auto& building : m_buildings) {
        building.--.output_for_object(
          indexer, num_vertices, vertices, faces, building_index);
        ++building_index;
      }
      return std::make_pair(vertices, faces);
    } */

  private:
    const Data_structure& m_data;
    const Points& m_interior_points;
    const Points& m_boundary_points;
    const std::size_t m_site_index;
    Plane_3 m_ground_plane;

    std::vector<Building> m_buildings;
    Points_2 m_boundary_points_2;
    std::vector< std::vector<std::size_t> > m_wall_points_2;
    std::vector<Segment_2> m_approximate_boundaries_2;

    void create_ground_plane() {

      FT minz = internal::max_value<FT>();
      for (const std::size_t idx : m_interior_points)
        minz = CGAL::min(minz, get(m_data.point_map_3, idx).z());
      for (const std::size_t idx : m_boundary_points)
        minz = CGAL::min(minz, get(m_data.point_map_3, idx).z());
      m_ground_plane = 
      Plane_3(Point_3(FT(0), FT(0), minz), Vector_3(FT(0), FT(0), FT(1)));
    }

    void extract_boundary_points_2(
      const FT alpha_shape_size_2, 
      const FT grid_cell_width_2) {

      m_boundary_points_2.clear();
      const FT sampling_2 = grid_cell_width_2;
      apply_alpha_shapes_filtering_2(alpha_shape_size_2, sampling_2);
      apply_grid_based_filtering_2(grid_cell_width_2);
      apply_thinning_2(m_data.parameters.scale / FT(2));
    }

    void apply_alpha_shapes_filtering_2(
      const FT alpha_shape_size_2, 
      const FT sampling_2) {

      const std::size_t numi = m_interior_points.size();
      const std::size_t numb = m_boundary_points.size();

      if (numi < 3 && numb < 3) return;
      Alpha_shapes_filtering_2 filtering(alpha_shape_size_2);

      if (numi >= 3)
        filtering.add_points(m_interior_points, m_data.point_map_3);
      filtering.get_filtered_points(sampling_2, m_boundary_points_2);
      if (numb >= 3)
        filtering.add_points(m_boundary_points, m_data.point_map_3);
      filtering.get_filtered_points(sampling_2, m_boundary_points_2);
    }

    void apply_grid_based_filtering_2(
      const FT grid_cell_width_2) {

      if (m_boundary_points_2.empty()) return;
      const Grid_based_filtering_2 filtering(grid_cell_width_2);
			filtering.apply(m_boundary_points_2);
    }

    void apply_thinning_2(
      const FT scale) {

      Identity_map identity_map;
      Sphere_neighbor_query neighbor_query(
        m_boundary_points_2, scale, identity_map);
      const Thinning_2 thinning(neighbor_query);
      thinning.apply(m_boundary_points_2);
    }

    void extract_wall_points_2(
      const FT region_growing_scale_2,
      const FT region_growing_noise_level_2,
      const FT region_growing_angle_2,
      const FT region_growing_min_length_2) {

      if (m_boundary_points_2.empty())
        return;

      m_wall_points_2.clear();
      
      Identity_map identity_map;
      K_neighbor_query neighbor_query(
        m_boundary_points_2, region_growing_scale_2, identity_map);

      Vectors_2 normals;
      Normal_estimator_2 estimator(
        m_boundary_points_2, neighbor_query, identity_map);
      estimator.get_normals(normals);

      CGAL_assertion(m_boundary_points_2.size() == normals.size());
      Pair_range_2 range;
      range.reserve(m_boundary_points_2.size());
      for (std::size_t i = 0; i < m_boundary_points_2.size(); ++i)
        range.push_back(std::make_pair(m_boundary_points_2[i], normals[i]));

      First_of_pair_map point_map;
      Second_of_pair_map normal_map;
      LSLF_region region(
        range, 
        region_growing_noise_level_2,
        region_growing_angle_2,
        region_growing_min_length_2,
        point_map,
        normal_map);

      LSLF_sorting sorting(
        m_boundary_points_2, neighbor_query, identity_map);
      sorting.sort();

      Region_growing_2 region_growing(
        m_boundary_points_2,
        neighbor_query,
        region,
        sorting.seed_map());
      region_growing.detect(std::back_inserter(m_wall_points_2));
    }

    void compute_approximate_boundaries() {

      if (m_wall_points_2.empty())
        return;

      m_approximate_boundaries_2.clear();
      m_approximate_boundaries_2.reserve(m_wall_points_2.size());

      Identity_map identity_map;
      Boundary_point_map_2 point_map_2(
        m_boundary_points_2, identity_map);

      Line_2 line; Point_2 p, q; 
      std::vector<std::size_t> indices;
      for (const auto& item_range : m_wall_points_2) {
        indices.clear();
        for (std::size_t i = 0; i < item_range.size(); ++i)
          indices.push_back(i);
        internal::line_from_points_2(
          item_range, point_map_2, line);
        internal::boundary_points_on_line_2(
          item_range, point_map_2, indices, line, p, q);
        m_approximate_boundaries_2.push_back(Segment_2(p, q));
      }
      CGAL_assertion(
        m_approximate_boundaries_2.size() == m_wall_points_2.size());
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_BUILDINGS_SITE_H
