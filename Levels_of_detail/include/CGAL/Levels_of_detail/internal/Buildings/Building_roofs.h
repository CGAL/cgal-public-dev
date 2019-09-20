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

// Partitioning.
#include <CGAL/Levels_of_detail/internal/Partitioning/Kinetic_partitioning_2.h>
#include <CGAL/Levels_of_detail/internal/Partitioning/Kinetic_partitioning_3.h>

// Visibility.
#include <CGAL/Levels_of_detail/internal/Visibility/Visibility_3.h>
#include <CGAL/Levels_of_detail/internal/Visibility/Roof_visibility_2.h>

// Graphcut.
#include <CGAL/Levels_of_detail/internal/Graphcut/Graphcut.h>
#include <CGAL/Levels_of_detail/internal/Graphcut/Roof_graphcut_2.h>

// Buildings.
#include <CGAL/Levels_of_detail/internal/Buildings/Building_ground_estimator.h>
#include <CGAL/Levels_of_detail/internal/Buildings/Building_walls_estimator.h>
#include <CGAL/Levels_of_detail/internal/Buildings/Building_roofs_estimator.h>
#include <CGAL/Levels_of_detail/internal/Buildings/Building_roofs_creator.h>
#include <CGAL/Levels_of_detail/internal/Buildings/Building_builder.h>
#include <CGAL/Levels_of_detail/internal/Buildings/Building_walls_creator.h>

// Experimental.
#include <CGAL/Levels_of_detail/internal/Experimental/Cloud_to_image_converter.h>

// Simplification.
#include <CGAL/Levels_of_detail/internal/Simplification/Generic_simplifier.h>

// Regularization.
#include <CGAL/Levels_of_detail/internal/Regularization/Regularization.h>
#include <CGAL/Levels_of_detail/internal/Regularization/Segment_regularizer.h>

// Testing.
#include "../../../../../test/Levels_of_detail/include/Saver.h"
#include "../../../../../test/Levels_of_detail/include/Utilities.h"

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
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Segment_2 = typename Traits::Segment_2;
    using Plane_3 = typename Traits::Plane_3;
    
    using Points_3 = std::vector<std::size_t>;
    using Point_map_3 = typename Data_structure::Point_map_3;

    using Indices = std::vector<std::size_t>;
    using Indexer = internal::Indexer<Point_3>;

    using Building = internal::Building<Traits>;
    using Triangulation = typename Building::Base::Triangulation::Delaunay;
    using Building_ground_estimator = internal::Building_ground_estimator<Traits, Triangulation>;
    using Approximate_face = internal::Partition_edge_3<Traits>;
    using Boundary = internal::Boundary<Traits>;
    using Building_walls_estimator = internal::Building_walls_estimator<Traits>;
    using Building_roofs_estimator = internal::Building_roofs_estimator<Traits, Points_3, Point_map_3>;
    using Building_roofs_creator = internal::Building_roofs_creator<Traits, Point_map_3>;
    using Building_walls_creator = internal::Building_walls_creator<Traits>;

    using Partition_3 = internal::Partition_3<Traits>;
    using Kinetic_partitioning_3 = internal::Kinetic_partitioning_3<Traits>;
    using Graphcut_3 = internal::Graphcut<Traits, Partition_3>;

    using Building_builder_3 = internal::Building_builder<Traits, Partition_3, Points_3, Point_map_3>;
    
    using Generic_simplifier = internal::Generic_simplifier<Traits, Point_map_3>;
    using Regularization = internal::Regularization<Traits>;

    using Segment_regularizer = internal::Segment_regularizer<Traits>;
    using Partition_2 = internal::Partition_2<Traits>;
    using Kinetic_partitioning_2 = internal::Kinetic_partitioning_2<Traits>;
    using Building_builder_2 = internal::Building_builder<Traits, Partition_2, Points_3, Point_map_3>;

    Building_roofs(
      const Data_structure& data,
      const Points_3& input,
      const std::vector<Point_3>& better_cluster,
      Building& building) : 
    m_data(data),
    m_input(input),
    m_better_cluster(better_cluster),
    m_building(building),
    m_empty(false) { 
      if (input.empty())
        m_empty = true;
    }

    void detect_roofs() {

      detect_roofs_3();
    }

    void detect_roofs_2() {

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

      extract_partitioning_constraints_2(
        m_data.parameters.buildings.grid_cell_width_2,
        m_data.parameters.buildings.alpha_shape_size_2,
        m_data.parameters.buildings.imagecut_beta_2,
        m_data.parameters.buildings.max_height_difference,
        m_data.parameters.buildings.image_noise_2,
        m_data.parameters.buildings.region_growing_min_length_2,
        m_data.parameters.buildings.regularization_angle_bound_2,
        m_data.parameters.buildings.regularization_ordinate_bound_2);
    }

    void detect_roofs_3() {
      
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

      make_approximate_bounds(
        m_data.parameters.buildings.grid_cell_width_2,
        m_data.parameters.buildings.alpha_shape_size_2,
        m_data.parameters.buildings.imagecut_beta_2,
        m_data.parameters.buildings.max_height_difference,
        m_data.parameters.buildings.image_noise_2,
        m_data.parameters.buildings.region_growing_min_length_2,
        m_data.parameters.buildings.regularization_angle_bound_2,
        m_data.parameters.buildings.regularization_ordinate_bound_2);
    }

    void compute_roofs() {

      compute_roofs_3();
    }

    void compute_roofs_2() {
      
      if (empty())
        return;

      partition_2(
        m_data.parameters.buildings.kinetic_min_face_width_2, 
        m_data.parameters.buildings.kinetic_max_intersections_2);

      compute_visibility_2();

      apply_graphcut_2(
        m_data.parameters.buildings.graphcut_beta_2);

      compute_roofs_and_corresponding_walls_2(
        m_data.parameters.buildings.max_height_difference);
    }

    void compute_roofs_3() {
      
      if (empty())
        return;

      partition_3(
        m_data.parameters.buildings.kinetic_max_intersections_3);

      compute_visibility_3();

      apply_graphcut_3(
        m_data.parameters.buildings.graphcut_beta_3);

      compute_roofs_and_corresponding_walls_3(
        m_data.parameters.buildings.max_height_difference);
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
          m_building_outer_walls.empty() &&
          m_building_inner_walls.empty() &&
          m_building_roofs.empty())
        return boost::none;

      m_building_ground.output_for_object( 
      indexer, num_vertices, vertices, faces, building_index);
      for (const auto& wall : m_building_outer_walls)
        wall.output_for_object( 
      indexer, num_vertices, vertices, faces, building_index);
      for (const auto& wall : m_building_inner_walls)
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
    const std::vector<Point_3>& m_better_cluster;
    Building& m_building;
    
    bool m_empty;
    Points_3 m_cluster;
    std::vector<Indices> m_roof_points_3;
    Indices m_unclassified_points_3;
    Approximate_face m_building_ground;
    std::vector<Approximate_face> m_building_outer_walls;
    std::vector<Approximate_face> m_building_inner_walls;
    std::vector<Approximate_face> m_building_roofs;
    Partition_3 m_partition_3;

    std::shared_ptr<Generic_simplifier> m_simplifier_ptr;
    
    std::vector<Segment_2> m_partitioning_constraints_2;
    std::vector< std::vector<Segment_2> > m_inner_wall_contours;
    std::vector< std::vector<Segment_2> > m_inner_roof_contours;
    std::vector<Plane_3> m_roof_planes;

    Partition_2 m_partition_2;

    void create_input_cluster_3(
      const FT region_growing_scale_3,
      const FT region_growing_angle_3) {
      
      if (empty()) return;

      const Building_roofs_creator creator(m_data.point_map_3);
      creator.create_cluster(
        m_input,
        region_growing_scale_3,
        region_growing_angle_3,
        m_cluster);
    }

    void extract_roof_regions_3(
      const FT region_growing_scale_3,
      const FT region_growing_noise_level_3,
      const FT region_growing_angle_3,
      const FT region_growing_min_area_3,
      const FT region_growing_distance_to_line_3,
      const FT alpha_shape_size_2) {
        
      if (empty()) return;
      
      const Building_roofs_creator creator(m_data.point_map_3);
      creator.create_roof_regions(
        m_cluster,
        region_growing_scale_3,
        region_growing_noise_level_3,
        region_growing_angle_3,
        region_growing_min_area_3,
        region_growing_distance_to_line_3,
        alpha_shape_size_2,
        m_roof_points_3,
        m_unclassified_points_3);
    }

    void make_approximate_bounds(
      const FT grid_cell_width_2,
      const FT alpha_shape_size_2,
      const FT imagecut_beta_2,
      const FT max_height_difference,
      const FT image_noise_2,
      const FT min_length_2,
      const FT angle_bound_2,
      const FT ordinate_bound_2) {
        
      // Roofs.
      bool success = add_approximate_roofs();
      if (!success) return;

      // Ground.
      success = add_approximate_ground();
      if (!success) return;

      // Walls.
      FT top_z = m_building.top_z;
      const FT bottom_z = m_building.bottom_z;
      CGAL_assertion(top_z > bottom_z);
      top_z -= (top_z - bottom_z) / FT(2);

      const Building_walls_estimator westimator(
        m_building.edges1,
        bottom_z, 
        top_z);

      success = add_outer_walls(westimator);
      if (!success) return;
      success = add_inner_walls_new_version(
        grid_cell_width_2,
        alpha_shape_size_2,
        imagecut_beta_2,
        max_height_difference,
        image_noise_2,
        min_length_2,
        angle_bound_2,
        ordinate_bound_2,
        westimator);
      if (!success) return;
    }

    bool add_approximate_roofs() {
      
      const Building_roofs_estimator restimator(
        m_cluster,
        m_data.point_map_3,
        m_roof_points_3);
      restimator.estimate(m_building_roofs);

      if (m_building_roofs.empty()) {
        m_empty = true; return false;
      }
      return true;
    }

    bool add_approximate_ground() {

      const FT bottom_z = m_building.bottom_z;
      const Building_ground_estimator gestimator(
        m_building.base1.triangulation.delaunay,
        bottom_z);
      gestimator.estimate(m_building_ground);
      return true;
    }

    bool add_outer_walls(
      const Building_walls_estimator& westimator) {
     
      westimator.estimate(m_building_outer_walls);

      // Safety feature to protect buildings with insufficient number of walls.
      if (m_building_outer_walls.size() < 3) {
        m_building_outer_walls.clear();
        m_building_outer_walls.reserve(m_building.edges1.size());

        Approximate_face wall;
        for (const auto& edge : m_building.edges1) {
          westimator.estimate_wall(edge, wall.polygon);
          m_building_outer_walls.push_back(wall);
        }
        CGAL_assertion(m_building_outer_walls.size() == m_building.edges1.size());
        return false;
      }
      return true;
    }

    void extract_partitioning_constraints_2(
      const FT grid_cell_width_2,
      const FT alpha_shape_size_2,
      const FT imagecut_beta_2,
      const FT max_height_difference,
      const FT image_noise_2,
      const FT min_length_2,
      const FT angle_bound_2,
      const FT ordinate_bound_2) {

      // Create image.
      create_image(
        grid_cell_width_2,
        alpha_shape_size_2,
        imagecut_beta_2,
        max_height_difference,
        image_noise_2,
        min_length_2);

      // Create wall contours.
      create_inner_contours(true, "wall", m_inner_wall_contours);

      // Create and regularize wall segments.
      std::vector<Segment_2> wall_segments;
      regularize_inner_contours(
        true, true,
        min_length_2, angle_bound_2, ordinate_bound_2,
        "wall",
        m_inner_wall_contours,
        wall_segments);

      // Create roof contours.
      create_inner_contours(false, "roof", m_inner_roof_contours);

      // Create and regularize roof segments.
      std::vector<Segment_2> roof_segments;
      regularize_inner_contours(
        false, false,
        min_length_2, angle_bound_2, ordinate_bound_2,
        "roof",
        m_inner_roof_contours,
        roof_segments);

      // Merge all segments.
      m_partitioning_constraints_2.clear();
      for (const auto& segment : wall_segments)
        m_partitioning_constraints_2.push_back(segment);
      for (const auto& segment : roof_segments)
        m_partitioning_constraints_2.push_back(segment);
      for (const auto& edge : m_building.edges1)
        m_partitioning_constraints_2.push_back(edge.segment);

      save_polylines(m_partitioning_constraints_2, 
      "/Users/monet/Documents/lod/logs/buildings/tmp/interior-edges-kinetic");
    }

    void regularize_inner_contours(
      const bool regularize_angles, 
      const bool regularize_ordinates,
      const FT min_length_2,
      const FT angle_bound_2,
      const FT ordinate_bound_2,
      const std::string name,
      std::vector< std::vector<Segment_2> >& contours,
      std::vector<Segment_2>& segments) {

      if (contours.empty()) return;

      Segment_regularizer regularizer(
        min_length_2, angle_bound_2, ordinate_bound_2);
      
      if (regularize_angles) {

        std::vector<Segment_2> outer_segments;
        for (const auto& edge : m_building.edges1)
          outer_segments.push_back(edge.segment);

        regularizer.compute_multiple_directions(
          outer_segments, contours);
        regularizer.regularize_contours(contours);
      }

      save_contours(contours,
      "/Users/monet/Documents/lod/logs/buildings/tmp/interior-" 
      + name + "-edges-after");

      segments.clear();
      for (const auto& contour : contours)
        for (const auto& segment : contour)
          segments.push_back(segment);

      if (regularize_ordinates)
        regularizer.merge_closest(segments);

      save_polylines(segments,
      "/Users/monet/Documents/lod/logs/buildings/tmp/interior-" 
      + name + "-edges-merged");
    }

    void create_inner_contours(
      const bool height_based,
      const std::string name,
      std::vector< std::vector<Segment_2> >& contours) {
      
      m_simplifier_ptr->create_inner_contours(height_based);
      m_simplifier_ptr->get_contours(
        contours);

      save_contours(contours,
      "/Users/monet/Documents/lod/logs/buildings/tmp/interior-" 
      + name + "-edges-before");
    }

    void create_image(
      const FT grid_cell_width_2,
      const FT alpha_shape_size_2,
      const FT imagecut_beta_2,
      const FT max_height_difference,
      const FT image_noise_2,
      const FT min_length_2) {

      m_simplifier_ptr = std::make_shared<Generic_simplifier>(
        m_cluster, 
        m_data.point_map_3,
        grid_cell_width_2,
        alpha_shape_size_2,
        imagecut_beta_2,
        max_height_difference,
        image_noise_2,
        min_length_2);

      m_simplifier_ptr->create_cluster_from_regions(
        m_roof_points_3,
        m_unclassified_points_3);
      m_simplifier_ptr->transform_cluster();
      m_simplifier_ptr->create_grid();
      m_simplifier_ptr->create_image(
        m_building.base1.triangulation, true);
    }

    void save_contours(
      const std::vector< std::vector<Segment_2> >& contours,
      const std::string path) {

      m_partitioning_constraints_2.clear();
      for (const auto& contour : contours)
        for (const auto& segment : contour)
          m_partitioning_constraints_2.push_back(segment);
      save_polylines(m_partitioning_constraints_2, path);
    }

    void partition_2(
      const FT kinetic_min_face_width_2,
      const std::size_t kinetic_max_intersections_2) {

      if (m_partitioning_constraints_2.empty()) return;
      
			const Kinetic_partitioning_2 kinetic(
        kinetic_min_face_width_2,
        kinetic_max_intersections_2);
			kinetic.compute(
        m_partitioning_constraints_2,
        m_partition_2);

      save_partition_2(
        "/Users/monet/Documents/lod/logs/buildings/tmp/partition_2", false);
      std::cout << "partition finished" << std::endl;
    }

    void save_partition_2(
      const std::string path,
      const bool with_roof_colors) {

      const FT z = m_building.bottom_z;
      std::size_t num_vertices = 0;
      internal::Indexer<Point_3> indexer;

      std::vector<Point_3> vertices; 
      std::vector<Indices> faces; 
      std::vector<CGAL::Color> fcolors;

      Polygon_inserter<Traits> inserter(faces, fcolors);
      auto output_vertices = std::back_inserter(vertices);
      auto output_faces = boost::make_function_output_iterator(inserter);

      for (const auto& face : m_partition_2.faces) {
        if (!with_roof_colors) {
          face.output_for_visibility(
            indexer, num_vertices, output_vertices, output_faces, z);
        } else {
          face.output_with_label_color(
            indexer, num_vertices, output_vertices, output_faces, z);
        }
      }
      
      Saver<Traits> saver;
      saver.export_polygon_soup(vertices, faces, fcolors, path);
    }

    void compute_visibility_2() {
      
      if (m_partition_2.empty()) return;

      std::vector<Point_3> points;
      std::vector<Indices> updated_regions;
      m_simplifier_ptr->get_points_for_visibility_3(
        m_building.base1.triangulation,
        m_cluster,
        m_roof_points_3,
        points,
        updated_regions,
        m_roof_planes);

      using Roof_visibility_2 = internal::Roof_visibility_2<Traits>;
      Roof_visibility_2 visibility(
        points,
        m_building,
        updated_regions);
      visibility.compute(m_partition_2);

      save_partition_2(
        "/Users/monet/Documents/lod/logs/buildings/tmp/visibility_inout_2", false);
      save_partition_2(
        "/Users/monet/Documents/lod/logs/buildings/tmp/visibility_roofs_2", true);
      std::cout << "visibility finished" << std::endl;
    }

    void apply_graphcut_2(
      const FT graphcut_beta_2) {

      if (m_partition_2.empty()) return;

      using Roof_graphcut_2 = internal::Roof_graphcut_2<Traits>;
      const Roof_graphcut_2 graphcut(graphcut_beta_2);
      graphcut.apply(m_partition_2);

      save_partition_2(
        "/Users/monet/Documents/lod/logs/buildings/tmp/graphcut_roofs_2", true);
      std::cout << "graphcut finished" << std::endl;
    }

    bool add_inner_walls_new_version(
      const FT grid_cell_width_2,
      const FT alpha_shape_size_2,
      const FT imagecut_beta_2,
      const FT max_height_difference,
      const FT image_noise_2,
      const FT min_length_2,
      const FT angle_bound_2,
      const FT ordinate_bound_2,
      const Building_walls_estimator& westimator) {

      // Create image.
      create_image(
        grid_cell_width_2,
        alpha_shape_size_2,
        imagecut_beta_2,
        max_height_difference,
        image_noise_2,
        min_length_2);

      // Create wall contours.
      create_inner_contours(true, "wall", m_inner_wall_contours);

      // Regularize wall contours.
      std::vector<Segment_2> segments;
      regularize_inner_contours(
        true, true,
        min_length_2, angle_bound_2, ordinate_bound_2,
        "wall",
        m_inner_wall_contours,
        segments);

      // Create walls.
      create_inner_walls_from_segments(segments, westimator);
      return true;
    }

    void save_polylines(
      const std::vector<Segment_2>& segments,
      const std::string name) {
      
      CGAL_assertion(segments.size() > 0);
      std::vector< std::vector<Point_3> > polylines(segments.size());
      for (std::size_t i = 0; i < segments.size(); ++i) {
        const Point_2& s = segments[i].source();
        const Point_2& t = segments[i].target();
        
        polylines[i].push_back(Point_3(s.x(), s.y(), FT(0)));
        polylines[i].push_back(Point_3(t.x(), t.y(), FT(0)));
      }
      
      Saver<Traits> saver;
      saver.export_polylines(polylines, name);
    }

    bool add_inner_walls_stable_version(
      const FT grid_cell_width_2,
      const FT alpha_shape_size_2,
      const FT imagecut_beta_2,
      const FT max_height_difference,
      const FT image_noise_2,
      const FT min_length_2,
      const FT angle_bound_2,
      const Building_walls_estimator& westimator) {
      
      // Create image.
      create_image(
        grid_cell_width_2,
        alpha_shape_size_2,
        imagecut_beta_2,
        max_height_difference,
        image_noise_2,
        min_length_2);

      // Create inner points.
      std::vector<Point_2> points;
      m_simplifier_ptr->get_inner_boundary_points_2(points);

      // Create region growing based segments.
      std::vector<Segment_2> segments;
      create_segments_region_growing(points, segments); 

      // Regularize segments with the global regularizer.
      regularize_segments_global(angle_bound_2, segments);

      // Create walls.
      create_inner_walls_from_segments(segments, westimator);
      return true;
    }

    void create_inner_walls_from_segments(
      const std::vector<Segment_2>& segments,
      const Building_walls_estimator& westimator) {

      m_building_inner_walls.clear();
      m_building_inner_walls.reserve(segments.size());

      Approximate_face wall; Boundary boundary;
      for (const auto& segment : segments) {
        boundary.segment = segment;
        westimator.estimate_wall(boundary, wall.polygon);
        m_building_inner_walls.push_back(wall);
      }
    }

    void create_segments_region_growing(
      const std::vector<Point_2>& points, 
      std::vector<Segment_2>& segments) {

      Building_walls_creator creator(points);
      std::vector<Indices> regions;
      creator.create_wall_regions(
        m_data.parameters.buildings.region_growing_scale_2 / FT(2),
        m_data.parameters.buildings.region_growing_noise_level_2 / FT(2),
        m_data.parameters.buildings.region_growing_angle_2,
        m_data.parameters.buildings.region_growing_min_length_2,
        regions);

      creator.create_boundaries(regions, segments);
      save_polylines(segments, 
      "/Users/monet/Documents/lod/logs/buildings/tmp/interior-edges-before");
    }

    void regularize_segments_global(
      const FT angle_bound_2,
      std::vector<Segment_2>& segments) {

      CGAL_assertion(segments.size() >= 2);
      Regularization regularization;
      regularization.regularize_angles(
        segments,
        angle_bound_2);

      save_polylines(segments, 
      "/Users/monet/Documents/lod/logs/buildings/tmp/interior-edges-after");
    }

    void partition_3(
      const std::size_t kinetic_max_intersections_3) {

      Kinetic_partitioning_3 kinetic(
        m_building_outer_walls,
        m_building_inner_walls,
        m_building_roofs,
        m_building_ground,
        kinetic_max_intersections_3);
      kinetic.compute(m_partition_3);

      std::cout << "kinetic finished" << std::endl;
    }

    void compute_visibility_3() {

      std::vector<Point_3> points;
      std::vector<Indices> updated_regions;
      std::vector<Plane_3> planes;
      m_simplifier_ptr->get_points_for_visibility_3(
        m_building.base1.triangulation,
        m_cluster,
        m_roof_points_3,
        points,
        updated_regions,
        planes);
      
      using Identity_map_3 = CGAL::Identity_property_map<Point_3>;
      using Visibility_3 = internal::Visibility_3<Traits, std::vector<Point_3>, Identity_map_3>;
      Identity_map_3 identity_map_3;

      if (m_partition_3.empty()) return;
      Visibility_3 visibility(
        points,
        identity_map_3, 
        m_building,
        updated_regions);
      visibility.compute(m_partition_3);

      std::cout << "visibility finished" << std::endl;
    }

    void apply_graphcut_3(
      const FT graphcut_beta_3) {

      if (m_partition_3.empty()) return;
      const Graphcut_3 graphcut(graphcut_beta_3);
      graphcut.apply(m_partition_3);

      std::cout << "graphcut finished" << std::endl;
    }

    void compute_roofs_and_corresponding_walls_3(
      const FT max_height_difference) {

      if (m_partition_3.empty()) return;
      const FT height_threshold = max_height_difference / FT(4);
      const Building_builder_3 builder(m_partition_3, height_threshold);
      builder.add_lod2_from_partition_3(m_building);

      if (m_building.roofs2.empty() || m_building.walls2.empty())
        m_empty = true;

      std::cout << "builder finished" << std::endl;
    }

    void compute_roofs_and_corresponding_walls_2(
      const FT max_height_difference) {

      if (m_partition_2.empty()) return;
      const Building_builder_2 builder(m_partition_2, max_height_difference);
      builder.add_lod2_from_partition_2(
        m_roof_planes, m_building);

      if (m_building.roofs2.empty() || m_building.walls2.empty())
        m_empty = true;

      std::cout << "builder finished: " << 
      m_building.walls2.size() << " " << m_building.roofs2.size() << std::endl;
    }
  };

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_INTERNAL_BUILDING_ROOFS_H
