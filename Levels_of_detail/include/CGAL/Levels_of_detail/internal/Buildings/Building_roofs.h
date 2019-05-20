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

#ifndef CGAL_LEVELS_OF_DETAIL_INTERNAL_BUILDING_ROOFS_H
#define CGAL_LEVELS_OF_DETAIL_INTERNAL_BUILDING_ROOFS_H

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

// Spacial search.
#include <CGAL/Levels_of_detail/internal/Spacial_search/Sphere_neighbor_query.h>

// Shape detection.
#include <CGAL/Levels_of_detail/internal/Shape_detection/Region_growing.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Estimate_normals_3.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Least_squares_plane_fit_region.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Least_squares_plane_fit_sorting.h>

// Partitioning.
#include <CGAL/Levels_of_detail/internal/Partitioning/Kinetic_partitioning_3.h>

// Visibility.
#include <CGAL/Levels_of_detail/internal/Visibility/Visibility_3.h>

// Graphcut.
#include <CGAL/Levels_of_detail/internal/Graphcut/Graphcut.h>

// Buildings.
#include <CGAL/Levels_of_detail/internal/Buildings/Building_ground_estimator.h>
#include <CGAL/Levels_of_detail/internal/Buildings/Building_walls_estimator.h>
#include <CGAL/Levels_of_detail/internal/Buildings/Building_roofs_estimator.h>
#include <CGAL/Levels_of_detail/internal/Buildings/Building_builder.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename DataStructure>
  class Building_roofs {

  public:
    using Data_structure = DataStructure;

    using Traits = typename Data_structure::Traits;
    using Point_map = typename Data_structure::Point_map;

    using FT = typename Traits::FT;
    using Point_3 = typename Traits::Point_3;
    using Vector_3 = typename Traits::Vector_3;
    
    using Points_3 = std::vector<std::size_t>;
    using Point_map_3 = typename Data_structure::Point_map_3;

    using Indexer = internal::Indexer<Point_3>;

    using Vectors_3 = std::vector<Vector_3>;
    using Pair_item_3 = std::pair<Point_3, Vector_3>;
    using Pair_range_3 = std::vector<Pair_item_3>;
    using First_of_pair_map = CGAL::First_of_pair_property_map<Pair_item_3>;
    using Second_of_pair_map = CGAL::Second_of_pair_property_map<Pair_item_3>;

    using Sphere_neighbor_query =
    internal::Sphere_neighbor_query<Traits, Points_3, Point_map_3>;
    using Normal_estimator_3 = 
    internal::Estimate_normals_3<Traits, Points_3, Point_map_3, Sphere_neighbor_query>;
    using LSPF_region = 
    internal::Least_squares_plane_fit_region<Traits, Pair_range_3, First_of_pair_map, Second_of_pair_map>;
    using LSPF_sorting =
    internal::Least_squares_plane_fit_sorting<Traits, Points_3, Sphere_neighbor_query, Point_map_3>;
    using Region_growing_3 = 
    internal::Region_growing<Points_3, Sphere_neighbor_query, LSPF_region, typename LSPF_sorting::Seed_map>;

    using Building = internal::Building<Traits>;
    using Triangulation = typename Building::Base::Triangulation::Delaunay;
    using Building_ground_estimator = internal::Building_ground_estimator<Traits, Triangulation>;
    using Approximate_face = internal::Partition_edge_3<Traits>;
    using Building_walls_estimator = internal::Building_walls_estimator<Traits>;
    using Building_roofs_estimator = internal::Building_roofs_estimator<Traits, Points_3, Point_map_3>;

    using Partition_3 = internal::Partition_3<Traits>;
    using Kinetic_partitioning_3 = internal::Kinetic_partitioning_3<Traits>;

    using Visibility_3 = internal::Visibility_3<Traits, Points_3, Point_map_3>;
    using Graphcut_3 = internal::Graphcut<Traits, Partition_3>;

    using Building_builder = internal::Building_builder<Traits, Partition_3, Points_3, Point_map_3>;
    
    Building_roofs(
      const Data_structure& data,
      const Points_3& input,
      Building& building) : 
    m_data(data),
    m_input(input),
    m_building(building),
    m_empty(false) { 
      if (input.empty())
        m_empty = true;
    }

    void detect_roofs() {
      if (empty())
        return;

      create_input_cluster_3(
        m_data.parameters.buildings.region_growing_scale_3,
        m_data.parameters.buildings.region_growing_angle_3);
      extract_roof_regions_3(
        m_data.parameters.buildings.region_growing_scale_3,
        m_data.parameters.buildings.region_growing_noise_level_3,
        m_data.parameters.buildings.region_growing_angle_3,
        m_data.parameters.buildings.region_growing_min_area_3,
        m_data.parameters.buildings.region_growing_distance_to_line_3,
        m_data.parameters.buildings.alpha_shape_size_2);
      make_approximate_bounds();
    }

    void compute_roofs() {
      if (empty())
        return;

      partition_3(
        m_data.parameters.buildings.kinetic_max_intersections_3);
      compute_visibility_3(
        m_data.parameters.buildings.visibility_scale_3);
      apply_graphcut_3(
        m_data.parameters.buildings.graphcut_beta_3);
      // compute_roofs_and_corresponding_walls();
    }

    void set_flat_roofs() {

      m_building.edges2 = m_building.edges1;
      m_building.base2 = m_building.base1;
      m_building.walls2 = m_building.walls1;
      m_building.roofs2 = m_building.roofs1;
    }

    const bool empty() const {
      return m_empty;
    }

    template<typename OutputIterator>
    boost::optional<OutputIterator> 
    get_roof_points(
      OutputIterator output,
      std::size_t& roof_index) const {

      if (m_roof_points_3.empty())
        return boost::none;

      for (std::size_t i = 0; i < m_roof_points_3.size(); ++i) {
        for (const std::size_t idx : m_roof_points_3[i]) {
          const Point_3& p = get(m_data.point_map_3, *(m_cluster.begin() + idx));
          *(output++) = std::make_pair(p, roof_index);
        }
        ++roof_index;
      }
      return output;
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> > 
    get_approximate_bounds(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      const std::size_t building_index) const {
      
      if (m_building_ground.polygon.empty() &&
          m_building_walls.empty() &&
          m_building_roofs.empty())
        return boost::none;

      m_building_ground.output_for_object( 
      indexer, num_vertices, vertices, faces, building_index);
      for (const auto& wall : m_building_walls)
        wall.output_for_object( 
      indexer, num_vertices, vertices, faces, building_index);
      for (const auto& roof : m_building_roofs)
        roof.output_for_object( 
      indexer, num_vertices, vertices, faces, building_index);
      return std::make_pair(vertices, faces);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> > 
    get_partitioning_3(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      const std::size_t building_index) const {
      
      if (m_partition_3.faces.empty())
        return boost::none;

      for (const auto& face : m_partition_3.faces)
        face.output_for_object(
          indexer, num_vertices, vertices, faces, building_index);
      return std::make_pair(vertices, faces);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> > 
    get_walls_corresponding_to_roofs(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      const std::size_t building_index) const {
      
      if (m_building.walls2.empty())
        return boost::none;

      for (const auto& wall : m_building.walls2)
        wall.output_for_object(
          indexer, num_vertices, vertices, faces, building_index);
      return std::make_pair(vertices, faces);
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    boost::optional< std::pair<VerticesOutputIterator, FacesOutputIterator> > 
    get_roofs(
      Indexer& indexer,
      std::size_t& num_vertices,
      VerticesOutputIterator vertices,
      FacesOutputIterator faces,
      const std::size_t building_index) const {
      
      if (m_building.roofs2.empty())
        return boost::none;

      for (const auto& roof : m_building.roofs2)
        roof.output_for_object(
          indexer, num_vertices, vertices, faces, building_index);
      return std::make_pair(vertices, faces);
    }

  private:
    const Data_structure& m_data;
    const Points_3& m_input;
    Building& m_building;
    
    bool m_empty;
    Points_3 m_cluster;
    std::vector< std::vector<std::size_t> > m_roof_points_3;
    Approximate_face m_building_ground;
    std::vector<Approximate_face> m_building_walls;
    std::vector<Approximate_face> m_building_roofs;
    Partition_3 m_partition_3;

    void create_input_cluster_3(
      const FT region_growing_scale_3,
      const FT region_growing_angle_3) {
      
      if (empty()) return;

      // Compute normals.
      Vectors_3 normals;
      Sphere_neighbor_query neighbor_query(
        m_input, region_growing_scale_3, m_data.point_map_3);
      Normal_estimator_3 estimator(
        m_input, neighbor_query, m_data.point_map_3);
      estimator.get_normals(normals);
      CGAL_assertion(normals.size() == m_input.size());

      // Remove vertical points.
      m_cluster.clear();
      const Vector_3 ref = Vector_3(FT(0), FT(0), FT(1));
      for (std::size_t i = 0; i < m_input.size(); ++i) {
        
        const auto& vec = normals[i];
        FT angle = angle_3d(vec, ref);
        if (angle > FT(90)) angle = FT(180) - angle;
        angle = FT(90) - angle;
        if (angle > region_growing_angle_3) 
          m_cluster.push_back(m_input[i]);
      }
    }

    void extract_roof_regions_3(
      const FT region_growing_scale_3,
      const FT region_growing_noise_level_3,
      const FT region_growing_angle_3,
      const FT region_growing_min_area_3,
      const FT region_growing_distance_to_line_3,
      const FT alpha_shape_size_2) {
        
      if (empty()) return;
      m_roof_points_3.clear();

      Sphere_neighbor_query neighbor_query(
        m_cluster, region_growing_scale_3, m_data.point_map_3);

      Vectors_3 normals;
      Normal_estimator_3 estimator(
        m_cluster, neighbor_query, m_data.point_map_3);
      estimator.get_normals(normals);

      CGAL_assertion(m_cluster.size() == normals.size());
      Pair_range_3 range;
      range.reserve(m_cluster.size());
      for (std::size_t i = 0; i < m_cluster.size(); ++i) {
        const Point_3& p = get(m_data.point_map_3, m_cluster[i]);
        range.push_back(std::make_pair(p, normals[i]));
      }

      First_of_pair_map point_map;
      Second_of_pair_map normal_map;
      LSPF_region region(
        range, 
        region_growing_noise_level_3,
        region_growing_angle_3,
        region_growing_min_area_3,
        region_growing_distance_to_line_3,
        alpha_shape_size_2,
        point_map,
        normal_map);

      LSPF_sorting sorting(
        m_cluster, neighbor_query, m_data.point_map_3);
      sorting.sort();

      Region_growing_3 region_growing(
        m_cluster,
        neighbor_query,
        region,
        sorting.seed_map());
      region_growing.detect(std::back_inserter(m_roof_points_3));
    }

    void remove_vertical_roof_regions_3(
      const FT region_growing_angle_3) {
      
      std::vector< std::vector<std::size_t> > regions = m_roof_points_3;
      m_roof_points_3.clear();
      std::vector<Point_3> points;
      for (const auto& region : regions) {
        points.clear();
        for (const std::size_t idx : region) {
          const auto& p = get(m_data.point_map_3, *(m_cluster.begin() + idx));
          points.push_back(p);
        }
        if (!internal::is_vertical_polygon(points, region_growing_angle_3))
          m_roof_points_3.push_back(region);
      }
    }

    void make_approximate_bounds() {
        
      // Roofs.
      const Building_roofs_estimator restimator(
        m_cluster,
        m_data.point_map_3,
        m_roof_points_3);
      restimator.estimate(m_building_roofs);

      if (m_building_roofs.empty()) {
        m_empty = true;
        return;
      }

      // Ground.
      const FT bottom_z = m_building.bottom_z;
      const Building_ground_estimator gestimator(
        m_building.base1.triangulation.delaunay,
        bottom_z);
      gestimator.estimate(m_building_ground);

      // Walls.
      FT top_z = m_building.top_z;
      CGAL_assertion(top_z > bottom_z);
      top_z -= (top_z - bottom_z) / FT(2);

      const Building_walls_estimator westimator(
        m_building.edges1,
        bottom_z, 
        top_z);
      westimator.estimate(m_building_walls);

      // Safety feature to protect buildings with insufficient number of walls.
      if (m_building_walls.size() < 3) {
        m_building_walls.clear();
        m_building_walls.reserve(m_building.edges1.size());

        Approximate_face wall;
        for (const auto& edge : m_building.edges1) {
          westimator.estimate_wall(edge, wall.polygon);
          m_building_walls.push_back(wall);
        }
        CGAL_assertion(m_building_walls.size() == m_building.edges1.size());
      }
    }

    void partition_3(
      const std::size_t kinetic_max_intersections_3) {

      const Kinetic_partitioning_3 kinetic(
        m_building_walls,
        m_building_roofs,
        m_building_ground,
        kinetic_max_intersections_3);
      kinetic.compute(m_partition_3);

      // std::cout << "kinetic finished" << std::endl;
    }

    void compute_visibility_3(
      const FT visibility_scale_3) {

      if (m_partition_3.empty()) return;
      const Visibility_3 visibility(
        m_cluster,
        m_data.point_map_3, 
        m_building,
        m_roof_points_3,
        visibility_scale_3);
      visibility.compute(m_partition_3);

      // std::cout << "visibility finished" << std::endl;
    }

    void apply_graphcut_3(
      const FT graphcut_beta_3) {

      if (m_partition_3.empty()) return;
      const Graphcut_3 graphcut(graphcut_beta_3);
      graphcut.apply(m_partition_3);

      // std::cout << "graphcut finished" << std::endl;
    }

    void compute_roofs_and_corresponding_walls() {

      if (m_partition_3.empty()) return;
      const FT distance_threshold = m_data.parameters.scale / FT(4);
      const Building_builder builder(m_partition_3, distance_threshold);
      builder.add_lod2(m_building);

      if (m_building.roofs2.empty() || m_building.walls2.empty())
        m_empty = true;

      // std::cout << "builder finished" << std::endl;
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_BUILDING_ROOFS_H
