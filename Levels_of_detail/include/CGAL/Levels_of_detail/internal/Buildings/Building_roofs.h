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
#include <CGAL/Levels_of_detail/internal/Partitioning/Partition_23_adapter.h>
#include <CGAL/Levels_of_detail/internal/Partitioning/Kinetic_partitioning_2.h>
#include <CGAL/Levels_of_detail/internal/Partitioning/Kinetic_partitioning_3.h>
#include <CGAL/Levels_of_detail/internal/Deprecated/Partition_builder_from_image_2_v2.h>
#include <CGAL/Levels_of_detail/internal/Partitioning/Partition_builder_from_image_2.h>

// Visibility.
#include <CGAL/Levels_of_detail/internal/Visibility/Visibility_3.h>
#include <CGAL/Levels_of_detail/internal/Visibility/Roof_visibility_2.h>
#include <CGAL/Levels_of_detail/internal/Visibility/Roof_visibility_3.h>

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

// Simplification.
#include <CGAL/Levels_of_detail/internal/Simplification/Generic_simplifier.h>

// Regularization.
#include <CGAL/Levels_of_detail/internal/Regularization/Segment_merger.h>
#include <CGAL/Levels_of_detail/internal/Regularization/Segment_regularizer.h>
#include <CGAL/Levels_of_detail/internal/Regularization/Regularization.h>

// Reconstruction.
#include <CGAL/Levels_of_detail/internal/Reconstruction/LOD2_image_reconstruction.h>

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
    using Partition_23_adapter = internal::Partition_23_adapter<Traits>;

    using Building_builder_3 = internal::Building_builder<Traits, Partition_3, Points_3, Point_map_3>;

    using Generic_simplifier = internal::Generic_simplifier<Traits, Point_map_3>;
    using Regularization = internal::Regularization<Traits>;

    using Segment_merger = internal::Segment_merger<Traits>;
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
    m_empty(false),
    m_num_roofs(std::size_t(-1)),
    m_num_labels(std::size_t(-1)),
    m_is_image_created(false) {

      if (input.empty())
        m_empty = true;
    }

    void detect_roofs() {

      /* detect_roofs_2_v1(true); */

      /* detect_roofs_2_v2(); */

      /* detect_roofs_2_v3(); */

      /* detect_roofs_2_v4(); */

      detect_roofs_3();

      /* detect_roofs_23(); */
    }

    void detect_roofs_2_v1(const bool use_default) {

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
        m_data.parameters.buildings.alpha_shape_size_2,
        m_data.parameters.buildings.max_height_difference);

      const bool use_lidar = true;
      extract_partitioning_constraints_2(
        m_data.parameters.buildings.grid_cell_width_2,
        m_data.parameters.buildings.alpha_shape_size_2,
        m_data.parameters.buildings.imagecut_beta_2,
        m_data.parameters.buildings.max_height_difference,
        m_data.parameters.buildings.image_noise_2,
        m_data.parameters.buildings.regularization_min_length_2,
        m_data.parameters.buildings.regularization_angle_bound_2,
        m_data.parameters.buildings.regularization_ordinate_bound_2,
        use_lidar, use_default);
    }

    void detect_roofs_2_v2() {

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
        m_data.parameters.buildings.alpha_shape_size_2,
        m_data.parameters.buildings.max_height_difference);

      m_is_image_created = create_image(
        m_data.parameters.buildings.grid_cell_width_2,
        m_data.parameters.buildings.alpha_shape_size_2,
        m_data.parameters.buildings.imagecut_beta_2,
        m_data.parameters.buildings.max_height_difference,
        m_data.parameters.buildings.image_noise_2,
        m_data.parameters.buildings.regularization_min_length_2,
        m_data.parameters.noise_level,
        m_data.parameters.buildings.region_growing_scale_3,
        m_data.parameters.buildings.region_growing_angle_3,
        false);
    }

    void detect_roofs_2_v3() {

      detect_roofs_2_v2();
    }

    void detect_roofs_2_v4() {

      detect_roofs_2_v2();
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
        m_data.parameters.buildings.alpha_shape_size_2,
        m_data.parameters.buildings.max_height_difference);

      const bool use_lidar = true;
      make_approximate_bounds(
        m_data.parameters.buildings.grid_cell_width_2,
        m_data.parameters.buildings.alpha_shape_size_2,
        m_data.parameters.buildings.imagecut_beta_2,
        m_data.parameters.buildings.max_height_difference,
        m_data.parameters.buildings.image_noise_2,
        m_data.parameters.buildings.regularization_min_length_2,
        m_data.parameters.buildings.regularization_angle_bound_2,
        m_data.parameters.buildings.regularization_ordinate_bound_2,
        use_lidar);
    }

    void detect_roofs_23() {

      if (empty())
        return;

      const bool use_default = true;
      detect_roofs_2_v1(use_default);
      const bool construct = false;
      const bool with_gc = true;
      compute_roofs_2_v1(construct, with_gc);
    }

    void compute_roofs() {

      /* compute_roofs_2_v1(true, true); */

      /* compute_roofs_2_v2(); */

      /* compute_roofs_2_v3(); */

      /* compute_roofs_2_v4(); */

      compute_roofs_3(true);

      /* compute_roofs_23(); */
    }

    void compute_roofs_2_v1(
      const bool construct,
      const bool with_gc) {

      if (empty())
        return;

      partition_2(
        m_data.parameters.buildings.kinetic_min_face_width_2,
        m_data.parameters.buildings.kinetic_max_intersections_2);

      compute_visibility_2();

      if (with_gc)
        apply_graphcut_2(
         m_data.parameters.buildings.graphcut_beta_2);

      if (construct)
        compute_roofs_and_corresponding_walls_2(
          m_data.parameters.buildings.max_height_difference);
    }

    void compute_roofs_2_v2() {

      if (empty())
        return;

      partition_2_from_image_v2(
        m_data.parameters.buildings.regularization_min_length_2,
        m_data.parameters.buildings.regularization_angle_bound_2,
        m_data.parameters.buildings.regularization_ordinate_bound_2);

      compute_roofs_and_corresponding_walls_2(
        m_data.parameters.buildings.max_height_difference);
    }

    void compute_roofs_2_v3() {

      if (empty())
        return;

      partition_2_from_image_v3(
        m_data.parameters.noise_level,
        m_data.parameters.buildings.regularization_min_length_2,
        m_data.parameters.buildings.regularization_angle_bound_2,
        m_data.parameters.buildings.regularization_ordinate_bound_2);

      compute_roofs_and_corresponding_walls_2(
        m_data.parameters.buildings.max_height_difference);
    }

    void compute_roofs_2_v4() {

      if (empty())
        return;

      partition_2_from_image_v4(
        m_data.parameters.noise_level / FT(2),
        m_data.parameters.buildings.regularization_min_length_2,
        m_data.parameters.buildings.regularization_angle_bound_2,
        m_data.parameters.buildings.regularization_ordinate_bound_2,
        m_data.parameters.buildings.max_height_difference,
        m_data.parameters.buildings.graphcut_beta_2);
    }

    void compute_roofs_3(const bool use_image) {

      if (empty())
        return;

      partition_3(
        m_data.parameters.buildings.kinetic_max_intersections_3);

      compute_visibility_3(use_image);

      apply_graphcut_3(
        m_data.parameters.buildings.graphcut_beta_3);

      compute_roofs_and_corresponding_walls_3(
        m_data.parameters.buildings.max_height_difference);
    }

    void compute_roofs_23() {

      if (empty())
        return;

      compute_partition_data_23(
        m_data.parameters.buildings.regularization_ordinate_bound_2,
        m_data.parameters.buildings.max_height_difference);

      const bool use_image = false;
      compute_roofs_3(use_image);
    }

    void set_flat_roofs() {

      m_building.edges2 = m_building.edges1;
      m_building.base2  = m_building.base1;
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
    std::shared_ptr<Partition_23_adapter> m_partion_23_adapter_ptr;

    std::vector<Segment_2> m_partitioning_constraints_2;
    std::vector< std::vector<Segment_2> > m_inner_wall_contours;
    std::vector< std::vector<Segment_2> > m_inner_roof_contours;
    std::vector<Plane_3> m_roof_planes;

    Partition_2 m_partition_2;

    std::size_t m_num_roofs;
    std::size_t m_num_labels;
    bool m_is_image_created;

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
      const FT alpha_shape_size_2,
      const FT max_height_difference) {

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
        max_height_difference,
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
      const FT ordinate_bound_2,
      const bool use_lidar) {

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
        use_lidar,
        westimator);
      if (!success) return;

      merge_inner_with_outer_walls(
        ordinate_bound_2,
        westimator);
    }

    void merge_inner_with_outer_walls(
      const FT ordinate_bound_2,
      const Building_walls_estimator& westimator) {

      std::vector<Segment_2> segments;
      for (const auto& face : m_building_inner_walls) {
        const auto& p0 = face.polygon[0];
        const auto& p1 = face.polygon[1];

        const Point_2 s = Point_2(p0.x(), p0.y());
        const Point_2 t = Point_2(p1.x(), p1.y());
        segments.push_back(Segment_2(s, t));
      }

      std::vector<Segment_2> outer;
      for (const auto& face : m_building_outer_walls) {
        const auto& p0 = face.polygon[0];
        const auto& p1 = face.polygon[1];

        const Point_2 s = Point_2(p0.x(), p0.y());
        const Point_2 t = Point_2(p1.x(), p1.y());
        outer.push_back(Segment_2(s, t));
      }

      Segment_merger merger(ordinate_bound_2);
      merger.merge_segments_with_outer_boundary(outer, segments);

      m_building_outer_walls.clear();
      create_inner_walls_from_segments(segments, westimator);
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
      const FT ordinate_bound_2,
      const bool use_lidar,
      const bool use_default) {

      // Create image.
      const bool success = create_image(
        grid_cell_width_2,
        alpha_shape_size_2,
        imagecut_beta_2,
        max_height_difference,
        image_noise_2,
        min_length_2,
        FT(0), FT(0), FT(0),
        use_lidar);

      if (!success)
        return;

      // Create wall contours.
      create_inner_contours(true, "wall", m_inner_wall_contours);

      // Create and regularize wall segments.
      std::vector<Segment_2> wall_segments;
      regularize_inner_contours(
        true, true, true,
        min_length_2, angle_bound_2, ordinate_bound_2,
        "wall",
        m_inner_wall_contours,
        wall_segments);

      // Create roof contours.
      create_inner_contours(false, "roof", m_inner_roof_contours);
      /* create_roof_contours_test(m_inner_roof_contours); */

      // Create and regularize roof segments.
      std::vector<Segment_2> roof_segments;

      bool st1 = false, st2 = false, st3 = false;
      if (!use_default) {
        st1 = true; st2 = true; st3 = true;
      }

      regularize_inner_contours(
        st1, st2, st3,
        min_length_2, angle_bound_2, ordinate_bound_2,
        "roof",
        m_inner_roof_contours,
        roof_segments);

      // Merge all segments.
      Segment_merger merger(ordinate_bound_2);

      std::vector<Segment_2> segments;
      for (const auto& segment : wall_segments)
        segments.push_back(segment);

      if (!use_default) {
        for (const auto& segment : roof_segments)
          segments.push_back(segment);
        merger.merge_segments(segments);
      }

      std::vector<Segment_2> outer;
      for (const auto& edge : m_building.edges1)
        outer.push_back(edge.segment);

      merger.merge_segments_with_outer_boundary(
        outer, segments);

      if (use_default) {
        for (const auto& segment : roof_segments)
          segments.push_back(segment);
      }

      m_partitioning_constraints_2.clear();
      m_partitioning_constraints_2 = segments;

      for (auto& segment : m_partitioning_constraints_2) {
        const auto& p = segment.source();
        const auto& q = segment.target();

        const FT b1 = FT(1)  / FT(100);
        const FT b2 = FT(99) / FT(100);

        const FT x1 = b1 * p.x() + b2 * q.x();
        const FT x2 = b1 * q.x() + b2 * p.x();

        const FT y1 = b1 * p.y() + b2 * q.y();
        const FT y2 = b1 * q.y() + b2 * p.y();

        const Point_2 s = Point_2(x1, y1);
        const Point_2 t = Point_2(x2, y2);

        segment = Segment_2(s, t);
      }

      Saver<Traits> saver;
      saver.save_polylines(m_partitioning_constraints_2,
      "/Users/monet/Documents/gf/lod/logs/buildings/tmp/interior-edges-kinetic");
    }

    void regularize_inner_contours(
      const bool regularize_angles,
      const bool regularize_ordinates,
      const bool snap,
      const FT min_length_2,
      const FT angle_bound_2,
      const FT ordinate_bound_2,
      const std::string name,
      std::vector< std::vector<Segment_2> >& contours,
      std::vector<Segment_2>& segments) {

      if (contours.empty()) return;

      Segment_regularizer regularizer(
        min_length_2, angle_bound_2);

      std::vector<Segment_2> outer_segments;
      for (const auto& edge : m_building.edges1)
        outer_segments.push_back(edge.segment);

      if (regularize_angles) {
        regularizer.compute_multiple_directions(
          outer_segments, contours);
        regularizer.regularize_contours(contours);
      }

      save_contours(contours,
      "/Users/monet/Documents/gf/lod/logs/buildings/tmp/interior-"
      + name + "-edges-after");

      segments.clear();
      for (const auto& contour : contours)
        for (const auto& segment : contour)
          segments.push_back(segment);

      Segment_merger merger(ordinate_bound_2);

      if (regularize_ordinates)
        merger.merge_segments(segments);

      Saver<Traits> saver;
      saver.save_polylines(segments,
        "/Users/monet/Documents/gf/lod/logs/buildings/tmp/interior-"
      + name + "-edges-merged");

      if (snap)
        merger.snap_segments(outer_segments, segments);

      saver.save_polylines(segments,
      "/Users/monet/Documents/gf/lod/logs/buildings/tmp/interior-"
      + name + "-edges-snapped");
    }

    void create_inner_contours(
      const bool height_based,
      const std::string name,
      std::vector< std::vector<Segment_2> >& contours) {

      m_simplifier_ptr->create_inner_contours(height_based);
      m_simplifier_ptr->get_contours(
        contours);

      save_contours(contours,
      "/Users/monet/Documents/gf/lod/logs/buildings/tmp/interior-"
      + name + "-edges-before");
    }

    void create_roof_contours_test(
      std::vector< std::vector<Segment_2> >& contours) {

      std::vector<Segment_2> segments;
      for (const auto& edge : m_building.edges1)
        segments.push_back(edge.segment);

      m_simplifier_ptr->create_roof_contours_test(segments);
      m_simplifier_ptr->get_contours(
        contours);

      save_contours(contours,
      "/Users/monet/Documents/gf/lod/logs/buildings/tmp/interior-roof-edges-before");
    }

    bool create_image(
      const FT grid_cell_width_2,
      const FT alpha_shape_size_2,
      const FT imagecut_beta_2,
      const FT max_height_difference,
      const FT image_noise_2,
      const FT min_length_2,
      const FT noise_level,
      const FT region_growing_scale_3,
      const FT region_growing_angle_3,
      const bool use_lidar) {

      std::vector<Segment_2> segments;
      for (const auto& edge : m_building.edges1)
        segments.push_back(edge.segment);

      Saver<Traits> saver;
      saver.save_polylines(segments,
      "/Users/monet/Documents/gf/lod/logs/buildings/tmp/boundary-edges");

      if (m_cluster.empty() || m_roof_points_3.empty())
        return false;

      m_simplifier_ptr = std::make_shared<Generic_simplifier>(
        m_cluster,
        m_data.point_map_3,
        grid_cell_width_2,
        alpha_shape_size_2,
        imagecut_beta_2,
        max_height_difference,
        image_noise_2,
        min_length_2,
        noise_level,
        region_growing_scale_3,
        region_growing_angle_3,
        use_lidar);

      m_simplifier_ptr->create_cluster_from_regions(
        m_roof_points_3,
        m_unclassified_points_3);
      m_simplifier_ptr->transform_cluster();
      m_simplifier_ptr->create_grid();
      m_simplifier_ptr->create_image(
        m_building.base1.triangulation, true);

      return true;
    }

    void save_contours(
      const std::vector< std::vector<Segment_2> >& contours,
      const std::string path) {

      m_partitioning_constraints_2.clear();
      for (const auto& contour : contours)
        for (const auto& segment : contour)
          m_partitioning_constraints_2.push_back(segment);

      Saver<Traits> saver;
      saver.save_polylines(m_partitioning_constraints_2, path);
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
        "/Users/monet/Documents/gf/lod/logs/buildings/tmp/partition_2", false);
      std::cout << "partition finished" << std::endl;
    }

    void partition_2_from_image_v1() {

    }

    void partition_2_from_image_v2(
      const FT min_length_2,
      const FT angle_bound_2,
      const FT ordinate_bound_2) {

      if (!m_is_image_created) {
        m_partition_2.clear(); return;
      }

      std::vector<Segment_2> boundary;
      for (const auto& edge : m_building.edges1)
        boundary.push_back(edge.segment);

      m_partition_2.clear();
      using Partition_builder_from_image_2 = internal::Partition_builder_from_image_2_v2<
      Traits, std::shared_ptr<Generic_simplifier> >;

      Partition_builder_from_image_2 builder(
        boundary, m_building.base0.triangulation,
        m_simplifier_ptr, m_partition_2,
        min_length_2, angle_bound_2, ordinate_bound_2);

      builder.build();
      builder.add_constraints();
      builder.compute_visibility();
      builder.label_faces();
      builder.optimize();

      builder.get_roof_planes(m_roof_planes);
      std::cout << "partition finished" << std::endl;
    }

    void partition_2_from_image_v3(
      const FT noise_level,
      const FT min_length_2,
      const FT angle_bound_2,
      const FT ordinate_bound_2) {

      if (!m_is_image_created) {
        m_partition_2.clear(); return;
      }

      std::vector<Segment_2> boundary;
      for (const auto& edge : m_building.edges1)
        boundary.push_back(edge.segment);

      m_partition_2.clear();
      using Partition_builder_from_image_2 = internal::Partition_builder_from_image_2<
      Traits, std::shared_ptr<Generic_simplifier> >;

      Partition_builder_from_image_2 builder(
        boundary, m_building.base0.triangulation,
        m_simplifier_ptr, m_partition_2,
        noise_level, min_length_2, angle_bound_2, ordinate_bound_2);

      builder.build();
      builder.get_roof_planes(m_roof_planes);

      builder.create_triangulation();
      builder.compute_visibility();
      builder.label_faces();

      std::cout << "partition finished" << std::endl;
    }

    void partition_2_from_image_v4(
      const FT noise_level_2,
      const FT min_length_2,
      const FT angle_bound_2,
      const FT ordinate_bound_2,
      const FT max_height_difference,
      const FT beta) {

      if (!m_is_image_created) {
        m_partition_2.clear(); return;
      }

      std::vector<Segment_2> boundary;
      for (const auto& edge : m_building.edges1)
        boundary.push_back(edge.segment);

      m_partition_2.clear();
      using LOD2_image_reconstruction = internal::LOD2_image_reconstruction<
        Traits, std::shared_ptr<Generic_simplifier>, Point_map_3>;

      LOD2_image_reconstruction builder(
        m_input, m_data.point_map_3,
        boundary, m_building.directions,
        m_building.base0.triangulation,
        m_simplifier_ptr, m_partition_2,
        noise_level_2,
        min_length_2, angle_bound_2, ordinate_bound_2,
        max_height_difference, beta, m_building.top_z);

      builder.build();
      builder.simplify();
      builder.create_tree();

      builder.regularize();
      builder.get_roof_planes(m_roof_planes);
      builder.get_lod2(m_building);
      std::cout << "partition finished" << std::endl;
    }

    void save_partition_2(
      const std::string path,
      const bool with_roof_colors) {

      const FT z = FT(0);
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
        points,
        updated_regions,
        m_roof_planes);

      using Roof_visibility_2 = internal::Roof_visibility_2<Traits>;
      Roof_visibility_2 visibility(
        points,
        m_building,
        updated_regions);
      visibility.compute(m_partition_2);

      m_num_roofs  = updated_regions.size();
      m_num_labels = visibility.number_of_actual_roofs(m_partition_2);

      std::cout << "Num roofs/labels: "
        << m_num_roofs << "/" << m_num_labels << std::endl;

      if (m_num_roofs == 0 || m_num_labels == 0) {
        m_partition_2.clear(); return;
      }

      save_partition_2(
        "/Users/monet/Documents/gf/lod/logs/buildings/tmp/visibility_inout_2", false);
      save_partition_2(
        "/Users/monet/Documents/gf/lod/logs/buildings/tmp/visibility_roofs_2", true);
      std::cout << "visibility finished" << std::endl;
    }

    void apply_graphcut_2(
      const FT graphcut_beta_2) {

      if (m_partition_2.empty()) return;

      using Roof_graphcut_2 = internal::Roof_graphcut_2<Traits>;
      const Roof_graphcut_2 graphcut(
        m_num_roofs, m_num_labels, graphcut_beta_2);
      const bool success = graphcut.apply(m_partition_2);

      if (!success) {
        m_partition_2.clear(); return;
      }

      save_partition_2(
        "/Users/monet/Documents/gf/lod/logs/buildings/tmp/graphcut_roofs_2", true);
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
      const bool use_lidar,
      const Building_walls_estimator& westimator) {

      // Create image.
     const bool success = create_image(
        grid_cell_width_2,
        alpha_shape_size_2,
        imagecut_beta_2,
        max_height_difference,
        image_noise_2,
        min_length_2,
        FT(0), FT(0), FT(0),
        use_lidar);

      if (!success)
        return false;

      // Create wall contours.
      create_inner_contours(true, "wall", m_inner_wall_contours);

      // Regularize wall contours.
      std::vector<Segment_2> segments;
      regularize_inner_contours(
        true, true, true,
        min_length_2, angle_bound_2, ordinate_bound_2,
        "wall",
        m_inner_wall_contours,
        segments);

      // Create walls.
      create_inner_walls_from_segments(segments, westimator);
      return true;
    }

    bool add_inner_walls_stable_version(
      const FT grid_cell_width_2,
      const FT alpha_shape_size_2,
      const FT imagecut_beta_2,
      const FT max_height_difference,
      const FT image_noise_2,
      const FT min_length_2,
      const FT angle_bound_2,
      const bool use_lidar,
      const Building_walls_estimator& westimator) {

      // Create image.
      const bool success = create_image(
        grid_cell_width_2,
        alpha_shape_size_2,
        imagecut_beta_2,
        max_height_difference,
        image_noise_2,
        min_length_2,
        FT(0), FT(0), FT(0),
        use_lidar);

      if (!success)
        return false;

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
      Saver<Traits> saver;
      saver.save_polylines(segments,
      "/Users/monet/Documents/gf/lod/logs/buildings/tmp/interior-edges-before");
    }

    void regularize_segments_global(
      const FT angle_bound_2,
      std::vector<Segment_2>& segments) {

      CGAL_assertion(segments.size() >= 2);
      Regularization regularization;
      regularization.regularize_angles(
        segments,
        angle_bound_2);

      Saver<Traits> saver;
      saver.save_polylines(segments,
      "/Users/monet/Documents/gf/lod/logs/buildings/tmp/interior-edges-after");
    }

    void compute_partition_data_23(
      const FT ordinate_bound_2,
      const FT max_height_difference) {

      if (m_partition_2.empty()) {
        m_building_roofs.clear();
        m_building_inner_walls.clear();
        m_building_outer_walls.clear();
        return;
      }

      std::map<std::size_t, Plane_3> plane_map;
      m_simplifier_ptr->get_plane_map(plane_map);

      FT top_z = m_building.top_z;
      const FT bottom_z = m_building.bottom_z;
      CGAL_assertion(top_z > bottom_z);
      top_z -= (top_z - bottom_z) / FT(2);

      m_partion_23_adapter_ptr =
        std::make_shared<Partition_23_adapter>(
          plane_map,
          bottom_z, top_z,
          ordinate_bound_2,
          max_height_difference,
          m_partition_2);

      // Ground.
      bool success = add_approximate_ground();
      if (!success) return;

      // Roofs.
      success = m_partion_23_adapter_ptr->get_approximate_roofs(
        m_building_roofs);
      if (!success) return;

      // Outer walls.
      const Building_walls_estimator westimator(
        m_building.edges1,
        bottom_z,
        top_z);

      m_building_outer_walls.clear();
      success = add_outer_walls(westimator);
      if (!success) return;

      // Inner walls.
      success = m_partion_23_adapter_ptr->get_approximate_inner_walls(
        m_building_inner_walls);
      if (!success) return;

      merge_inner_with_outer_walls(
        ordinate_bound_2,
        westimator);
    }

    void partition_3(
      const std::size_t kinetic_max_intersections_3) {

      if (
        m_building_roofs.empty() &&
        m_building_inner_walls.empty() &&
        m_building_outer_walls.empty())
      return;

      Kinetic_partitioning_3 kinetic(
        m_building_outer_walls,
        m_building_inner_walls,
        m_building_roofs,
        m_building_ground,
        kinetic_max_intersections_3);
      kinetic.compute(m_partition_3);

      std::cout << "kinetic finished" << std::endl;
    }

    void compute_visibility_3(const bool use_image) {

      if (m_partition_3.empty()) return;
      std::vector<Point_3> points;
      std::vector<Indices> updated_regions;
      std::vector<Plane_3> stub;

      if (use_image)
        m_simplifier_ptr->get_points_for_visibility_3(
          m_building.base1.triangulation,
          points,
          updated_regions,
          stub);
      else {

        /*
        m_partion_23_adapter_ptr->get_points_for_visibility_3(
          points,
          updated_regions,
          stub); */

        compute_roof_visibility_3();
        return;
      }

      using Identity_map_3 = CGAL::Identity_property_map<Point_3>;
      using Visibility_3 = internal::Visibility_3<Traits, std::vector<Point_3>, Identity_map_3>;
      Identity_map_3 identity_map_3;

      Visibility_3 visibility(
        points,
        identity_map_3,
        m_building,
        updated_regions);
      visibility.compute(m_partition_3);

      std::cout << "visibility finished" << std::endl;
    }

    void compute_roof_visibility_3() {

      if (m_partition_2.empty() || m_partition_3.empty())
        return;

      using Visibility_3 = internal::Roof_visibility_3<Traits>;
      Visibility_3 visibility(
        m_partition_2, m_building.bottom_z, m_building.top_z);
      visibility.compute(m_partition_3);

      std::cout << "visibility finished" << std::endl;
    }

    void apply_graphcut_3(
      const FT graphcut_beta_3) {

      if (m_partition_3.empty()) return;
      const Graphcut_3 graphcut(
        graphcut_beta_3);
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
