#ifndef CGAL_LEVELS_OF_DETAIL_BUILDINGS_H
#define CGAL_LEVELS_OF_DETAIL_BUILDINGS_H

// STL includes.
#include <vector>
#include <algorithm>

// Boost includes.
#include <boost/iterator/transform_iterator.hpp>

// CGAL includes.
#include <CGAL/array.h>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utilities.h>

// Simplification.
#include <CGAL/Levels_of_detail/internal/Simplification/Thinning_2.h>
#include <CGAL/Levels_of_detail/internal/Simplification/Grid_based_filtering_2.h>
#include <CGAL/Levels_of_detail/internal/Simplification/Alpha_shapes_filtering_2.h>

// Shape detection.
#include <CGAL/Levels_of_detail/internal/Shape_detection/Region_growing.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Estimate_normals_2.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Points_2_fuzzy_sphere_connectivity.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Points_2_k_nearest_neighbors_connectivity.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Points_2_least_squares_line_fit_conditions.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Polygon_faces_2_stored_connectivity.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Polygon_faces_2_visibility_conditions.h>

// Partitioning.
#include <CGAL/Levels_of_detail/internal/Partitioning/Kinetic_partitioning_2.h>

// Visibility.
#include <CGAL/Levels_of_detail/internal/Visibility/K_nearest_neighbors_search_2.h>
#include <CGAL/Levels_of_detail/internal/Visibility/Visibility_2.h>

// Buildings.
#include <CGAL/Levels_of_detail/internal/Buildings/Building_footprints_2.h>
#include <CGAL/Levels_of_detail/internal/Buildings/Building_boundaries_2.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename DataStructure>
  class Buildings {

  public:
    using Data_structure = DataStructure;

    using Traits = typename Data_structure::Traits;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Line_2 = typename Traits::Line_2;
    using Segment_2 = typename Traits::Segment_2;

    using Building = typename Data_structure::Building;

    using Grid_based_filtering_2 = Grid_based_filtering_2<Traits>;
    using Alpha_shapes_filtering_2 = Alpha_shapes_filtering_2<Traits>;

    using Points_connectivity_2 = 
    Points_2_k_nearest_neighbors_connectivity<Traits>;
    using Normals_estimator_2 = 
    Estimate_normals_2<Traits, Points_connectivity_2>;
    using Points_conditions_2 = 
    Points_2_least_squares_line_fit_conditions<Traits>;
    using Points_region_growing_2 = 
    Region_growing<Points_connectivity_2, Points_conditions_2>;
    
    using Kinetic_partitioning_2 = Kinetic_partitioning_2<Traits>;

    using Input_range = typename Data_structure::Input_range;
    using Point_map = typename Data_structure::Point_map;
    using Visibility_map = typename Data_structure::Visibility_map;

    using Knn_2 =
    K_nearest_neighbors_search_2<Traits, Input_range, Point_map>;
    using Visibility_2 = 
    Visibility_2<Traits, Input_range, Knn_2, Visibility_map>;
    
    using Polygon_faces_connectivity_2 =
    Polygon_faces_2_stored_connectivity<Traits>;
    using Polygon_faces_conditions_2 =
    Polygon_faces_2_visibility_conditions<Traits>;
    using Polygon_faces_region_growing_2 = 
    Region_growing<Polygon_faces_connectivity_2, Polygon_faces_conditions_2>;

    using Building_footprints_2 = Building_footprints_2<Traits>;
    using Building_boundaries_2 = Building_boundaries_2<Traits>;

    Buildings(Data_structure& data_structure) :
    m_data(data_structure),
    m_has_exact_boundaries(false)
    { }
    
    // PROCESSING

    void detect_boundaries(
      const FT alpha_shape_size,
      const FT grid_cell_width,
      const FT region_growing_search_size,
      const FT region_growing_noise_level,
      const FT region_growing_angle,
      const FT region_growing_min_length) {

      if (m_data.verbose) 
        std::cout << std::endl << "- Detecting building boundaries" 
        << std::endl;

      extract_boundary_points_2(
        alpha_shape_size, 
        grid_cell_width);

      extract_wall_points_2(
        region_growing_search_size,
        region_growing_noise_level,
        region_growing_angle,
        region_growing_min_length);
    }

    void detect_footprints(
      const FT kinetic_min_face_width,
      const std::size_t kinetic_max_intersections,
      const std::size_t min_faces_per_building) {

      if (m_data.verbose) 
        std::cout << std::endl << "- Detecting building footprints" 
        << std::endl;

      partition_2(
        kinetic_min_face_width, 
        kinetic_max_intersections);

      compute_visibility_2();

      detect_building_footprints_2();

      finilize_buildings(min_faces_per_building);
    }

    // FLAGS.

    const bool has_exact_boundaries() const {
      return m_has_exact_boundaries;
    }

    // OUTPUT

    template<typename OutputIterator>
    void return_boundary_points(OutputIterator output) const {

      const auto& points = m_data.building_boundary_points_2;
      const auto& plane = m_data.ground_plane;

      CGAL_precondition(!points.empty());
      std::copy(
        boost::make_transform_iterator(
          points.begin(),
          internal::Point_3_from_point_2_and_plane<Traits>(plane)),
        boost::make_transform_iterator(
          points.end(),
          internal::Point_3_from_point_2_and_plane<Traits>(plane)),
        output);
    }

    template<typename OutputIterator>
    void return_wall_points(OutputIterator output) const {

      const auto& points = m_data.building_boundary_points_2;
      const auto& indices = m_data.building_boundary_indices_2;
      const auto& plane = m_data.ground_plane;

      std::vector<std::pair<Point_3, int> > data(points.size());
      for (std::size_t i = 0; i < points.size(); ++i)
        data[i] = std::make_pair(
                  internal::position_on_plane_3(
                    points[i], 
                    plane), -1);

        for(std::size_t i = 0; i < indices.size(); ++i)
          for(std::size_t j = 0; j < indices[i].size(); ++j)
            data[indices[i][j]].second = i;

        std::copy(data.begin(), data.end(), output);
    }

    template<typename OutputIterator>
    void return_approximate_boundary_edges(OutputIterator output) const {

      const auto& points = m_data.building_boundary_points_2;
      const auto& indices = m_data.building_boundary_indices_2;
      const auto& plane = m_data.ground_plane;

      std::copy(
        boost::make_transform_iterator(
          indices.begin(),
          internal::Segment_3_from_points_and_plane<Traits>(
            points, plane)),
        boost::make_transform_iterator(
          indices.end(),
          internal::Segment_3_from_points_and_plane<Traits>(
            points, plane)),
        output);
    }

    template<typename OutputIterator>
    void return_exact_boundary_edges(OutputIterator output) const {

      const auto& buildings = m_data.buildings;
      const auto& plane = m_data.ground_plane;

      for (std::size_t i = 0; i < buildings.size(); ++i) {
        const auto& edges = buildings[i].boundaries;

        for (std::size_t j = 0; j < edges.size(); ++j)
          *(output++) = 
          internal::segment_3_from_segment_2_and_plane(edges[j], plane);
      }
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void return_partitioning(
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces) const {

      const auto& faces = m_data.building_polygon_faces_2;
      const auto& plane = m_data.ground_plane;

      internal::Indexer<Point_2> indexer;
      std::size_t num_vertices = 0;

      for (std::size_t i = 0; i < faces.size(); ++i) {
        const auto& vertices = faces[i].vertices;

        std::vector<std::size_t> face;
        for (std::size_t j = 0; j < vertices.size(); ++j) {

          const std::size_t idx = indexer(vertices[j]);
          if (idx == num_vertices) {
              
            *(output_vertices++) = 
            internal::position_on_plane_3(vertices[j], plane);
            ++num_vertices;
          }
          face.push_back(idx);
        }
        *(output_faces++) = std::make_pair(face, faces[i].visibility);
      }
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void return_footprints(
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces) const {

      const auto& buildings = m_data.buildings;
      const auto& plane = m_data.ground_plane;
      
      internal::Indexer<Point_2> indexer;
      std::size_t num_vertices = 0;

      for (std::size_t i = 0; i < buildings.size(); ++i) {
        const auto& triangles = buildings[i].footprint;

        for (std::size_t j = 0; j < triangles.size(); ++j) {
          cpp11::array<std::size_t, 3> face;
          
          for (std::size_t k = 0; k < 3; ++k) {
            const auto& point = triangles[j][k];

            const std::size_t idx = indexer(point);
            if (idx == num_vertices) {

              *(output_vertices++) = 
              internal::position_on_plane_3(point, plane);
              ++num_vertices;
            }
            face[k] = idx;
          }
          *(output_faces++) = std::make_pair(face, i);
        }
      }
    }

  private:
    Data_structure& m_data;
    bool m_has_exact_boundaries;

    // Boundaries.
    void extract_boundary_points_2(
      const FT alpha_shape_size, 
      const FT grid_cell_width) {

      m_data.building_boundary_points_2.clear();
      const FT sampling = grid_cell_width;

      apply_alpha_shapes_filtering_2(alpha_shape_size, sampling);
      apply_grid_based_filtering_2(grid_cell_width);

      if (m_data.verbose)
        std::cout << "-> " << m_data.building_boundary_points_2.size()
        << " boundary point(s) extracted" << std::endl;
    }

    void apply_alpha_shapes_filtering_2(
      const FT alpha_shape_size, 
      const FT sampling) {

      if (m_data.verbose) 
        std::cout << "* alpha shapes filtering" 
        << std::endl;

      const std::size_t numb_pts = m_data.building_boundary_points().size();
      const std::size_t numi_pts = m_data.building_interior_points().size();

      CGAL_precondition(numb_pts >= 3 || numi_pts >= 3);
      Alpha_shapes_filtering_2 filtering(alpha_shape_size);

      const auto& boundary_points = m_data.building_boundary_points();
      const auto& interior_points = m_data.building_interior_points();

      if (numb_pts >= 3)
        filtering.add_points(
          boundary_points,
          m_data.point_map);

      if (numi_pts >= 3)
        filtering.add_points(
          interior_points,
          m_data.point_map);

      filtering.get_filtered_points(
        sampling, m_data.building_boundary_points_2);
    }

    void apply_grid_based_filtering_2(const FT grid_cell_width) {

        if (m_data.verbose) 
        std::cout << "* grid-based filtering" 
        << std::endl;

        const Grid_based_filtering_2 filtering(grid_cell_width);
				filtering.apply(m_data.building_boundary_points_2);
    }

    // Walls.
    void extract_wall_points_2(
      const FT region_growing_search_size,
      const FT region_growing_noise_level,
      const FT region_growing_angle,
      const FT region_growing_min_length) {

      if (m_data.verbose) 
        std::cout << "* region growing" 
        << std::endl;

      m_data.building_boundary_indices_2.clear();
      const auto& points = m_data.building_boundary_points_2;

      Points_connectivity_2 connectivity(
        points, 
        region_growing_search_size);

      Normals_estimator_2 estimator(
        points, 
        connectivity);

      Points_conditions_2 conditions(
        points, 
        estimator.normals(),
        region_growing_noise_level,
        region_growing_angle,
        region_growing_min_length);

      std::vector<std::size_t> indices(points.size());
      for (std::size_t i = 0; i < points.size(); ++i)
        indices[i] = i;
      std::stable_sort(
        indices.begin(), indices.end(), estimator.sorter());

      Points_region_growing_2 region_growing(
        indices,
        connectivity,
        conditions);

      region_growing.detect(
        m_data.building_boundary_indices_2);

      if (m_data.verbose) 
        std::cout << "-> " << m_data.building_boundary_indices_2.size()
        << " wall(s) extracted" 
        << std::endl;
    }

    // Partitioning.
    void partition_2(
      const FT kinetic_min_face_width,
      const std::size_t kinetic_max_intersections) {

      if (m_data.verbose) 
        std::cout << "* kinetic partitioning" 
      << std::endl;

      m_data.building_polygon_faces_2.clear();
      const auto& points = m_data.building_boundary_points_2;
      const auto& indices = m_data.building_boundary_indices_2;

      Line_2 line;
      Point_2 p, q;
      const std::size_t num_segments = indices.size();
      
      std::vector<Segment_2> segments(num_segments);
      for (std::size_t i = 0; i < num_segments; ++i) {

        line_from_points_2(points, indices[i], line);
        boundary_points_on_line_2(points, indices[i], line, p, q);

        segments[i] = Segment_2(p, q);
      }

			const Kinetic_partitioning_2 kinetic_partitioning_2(
        kinetic_max_intersections,
        kinetic_min_face_width);

			kinetic_partitioning_2.compute(
        segments,
        m_data.building_polygon_faces_2);
        
      if (m_data.verbose)
        std::cout << "-> " << m_data.building_polygon_faces_2.size()
        << " polygon face(s) created" 
        << std::endl;
    }

    // Visibility.
    void compute_visibility_2() {

      if (m_data.verbose) 
        std::cout << "* computing visibility" 
      << std::endl;

      Knn_2 connectivity(
        m_data.input_range, 
        m_data.point_map, 
        1);

      const Visibility_2 visibility(
        m_data.input_range, 
        connectivity, 
        m_data.visibility_map);

      visibility.compute(m_data.building_polygon_faces_2);
    }

    // Footprints.
    void detect_building_footprints_2() {

      if (m_data.verbose) 
        std::cout << "* detecting footprints" 
      << std::endl;

      m_data.building_footprints_2.clear();
      const auto& faces = m_data.building_polygon_faces_2;

      Polygon_faces_connectivity_2 connectivity(faces);
      Polygon_faces_conditions_2 conditions(faces);

      std::vector<std::size_t> indices(faces.size());
      for (std::size_t i = 0; i < faces.size(); ++i)
        indices[i] = i;

      Polygon_faces_region_growing_2 region_growing(
        indices,
        connectivity,
        conditions);

      region_growing.detect(
        m_data.building_footprints_2);

      if (m_data.verbose)
        std::cout << "-> " << m_data.building_footprints_2.size()
        << " building footprint(s) detected" 
        << std::endl;
    }

    // Final buildings.
    void finilize_buildings(const std::size_t min_faces_per_building) {

      m_data.buildings.clear();

      const auto& faces = m_data.building_polygon_faces_2;
      const auto& footprints = m_data.building_footprints_2;

      Building building;
      auto& triangles = building.footprint;
      auto& segments = building.boundaries;

      Building_footprints_2 footprints_extractor;
      Building_boundaries_2 boundaries_extractor;

      for (std::size_t i = 0; i < footprints.size(); ++i) {  
        const auto& indices = footprints[i];

        footprints_extractor.create_footprint_triangles(
          faces, 
          indices, 
          triangles);

        boundaries_extractor.create_boundary_segments(
          faces,
          indices,
          segments);

        CGAL_precondition(min_faces_per_building >= 1);
        if (triangles.size() >= min_faces_per_building && segments.size() >= 3)
          m_data.buildings.push_back(building);
      }
      m_has_exact_boundaries = true;
    }

  }; // Buildings

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_BUILDINGS_H
