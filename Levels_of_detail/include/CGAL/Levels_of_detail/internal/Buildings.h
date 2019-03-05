#ifndef CGAL_LEVELS_OF_DETAIL_BUILDINGS_H
#define CGAL_LEVELS_OF_DETAIL_BUILDINGS_H

// STL includes.
#include <tuple>
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
#include <CGAL/Levels_of_detail/internal/Shape_detection/Estimate_normals_3.h>

#include <CGAL/Levels_of_detail/internal/Shape_detection/Points_2_fuzzy_sphere_connectivity.h>

#include <CGAL/Levels_of_detail/internal/Shape_detection/Points_2_k_nearest_neighbors_connectivity.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Points_3_k_nearest_neighbors_connectivity.h>

#include <CGAL/Levels_of_detail/internal/Shape_detection/Points_2_least_squares_line_fit_conditions.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Points_3_least_squares_plane_fit_conditions.h>

#include <CGAL/Levels_of_detail/internal/Shape_detection/Polygon_faces_2_stored_connectivity.h>
#include <CGAL/Levels_of_detail/internal/Shape_detection/Polygon_faces_2_visibility_conditions.h>

#include <CGAL/Levels_of_detail/internal/Graphcut/Weight_quality_estimator_3.h>
#include <CGAL/Levels_of_detail/internal/Graphcut/Graphcut_3.h>

// Partitioning.
#include <CGAL/Levels_of_detail/internal/Partitioning/Kinetic_partitioning_2.h>
#include <CGAL/Levels_of_detail/internal/Partitioning/Kinetic_partitioning_3.h>

// Visibility.
#include <CGAL/Levels_of_detail/internal/Visibility/K_nearest_neighbors_search_2.h>
#include <CGAL/Levels_of_detail/internal/Visibility/Visibility_2.h>
#include <CGAL/Levels_of_detail/internal/Visibility/Visibility_3.h>

// Buildings.
#include <CGAL/Levels_of_detail/internal/Buildings/Building_footprints_2.h>
#include <CGAL/Levels_of_detail/internal/Buildings/Building_boundaries_2.h>
#include <CGAL/Levels_of_detail/internal/Buildings/Buildings_clustering.h>
#include <CGAL/Levels_of_detail/internal/Buildings/Building_height_estimator.h>
#include <CGAL/Levels_of_detail/internal/Buildings/Roof_cleaner.h>
#include <CGAL/Levels_of_detail/internal/Buildings/Roof_estimator.h>
#include <CGAL/Levels_of_detail/internal/Buildings/Wall_estimator.h>
#include <CGAL/Levels_of_detail/internal/Buildings/Building_ground_estimator.h>
#include <CGAL/Levels_of_detail/internal/Buildings/Roof_wall_extractor.h>

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
    using Filtered_range = typename Data_structure::Filtered_range;
    using Filtered_range_iterator = typename Data_structure::Filtered_range_iterator;
    using Cluster = std::vector<Filtered_range_iterator>;

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

    using Buildings_clustering = 
    Buildings_clustering<Traits, Filtered_range, Point_map>;

    using Building_height_estimator =
    Building_height_estimator<Traits, Filtered_range, Point_map>;

    struct Dereference_map {

      Point_map m_point_map;

      Dereference_map(const Point_map point_map) : 
      m_point_map(point_map)
      { }

      using key_type = Filtered_range_iterator;
      using value_type = typename Point_map::value_type;
      using reference = const value_type&;
      using category = boost::lvalue_property_map_tag;

      friend reference get(
      const Dereference_map& dereference_map, 
      const key_type& key) {
        
        return get(dereference_map.m_point_map, *key);
      }
    };

    using Points_connectivity_3 = 
    Points_3_k_nearest_neighbors_connectivity<Traits, Cluster, Dereference_map>;
    using Normals_estimator_3 = 
    Estimate_normals_3<Traits, Cluster, Dereference_map, Points_connectivity_3>;
    using Points_conditions_3 = 
    Points_3_least_squares_plane_fit_conditions<Traits, Cluster, Dereference_map>;
    using Points_region_growing_3 = 
    Region_growing<Points_connectivity_3, Points_conditions_3>;

    using Roof_cleaner = Roof_cleaner<Traits, Cluster, Dereference_map>;
    using Roof_estimator = Roof_estimator<Traits, Cluster, Dereference_map>;
    using Wall_estimator = Wall_estimator<Traits>;
    using Building_ground_estimator = Building_ground_estimator<Traits>;

    using Kinetic_partitioning_3 = Kinetic_partitioning_3<Traits>;
    using Visibility_3 = Visibility_3<Traits, Cluster, Dereference_map>;

    using Weight_quality_estimator_3 = Weight_quality_estimator_3<Traits>;
    using Graphcut_3 = Graphcut_3<Traits>;
    using Roof_wall_extractor = Roof_wall_extractor<Traits>;

    Buildings(Data_structure& data_structure) :
    m_data(data_structure),
    m_has_exact_boundaries(false),
    m_has_exact_roofs(false)
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

    void compute_footprints(
      const FT kinetic_min_face_width,
      const std::size_t kinetic_max_intersections,
      const std::size_t min_faces_per_building) {

      if (m_data.verbose) 
        std::cout << std::endl << "- Computing building footprints" 
        << std::endl;

      partition_2(
        kinetic_min_face_width, 
        kinetic_max_intersections);

      compute_visibility_2();

      compute_building_footprints_2();

      finilize_buildings(min_faces_per_building);
    }

    void extrude_footprints(const Extrusion_type extrusion_type) {
      
      if (m_data.verbose) 
        std::cout << std::endl << "- Extruding building footprints" 
        << std::endl;

      cluster_building_points();

      compute_building_heights(
        extrusion_type);
    }

    void detect_roofs(
      const FT region_growing_search_size,
      const FT region_growing_noise_level,
      const FT region_growing_angle,
      const FT region_growing_min_area,
      const FT min_size) {

      if (m_data.verbose) 
        std::cout << std::endl << "- Detecting building roofs" 
        << std::endl;

      extract_roof_regions_3(
        region_growing_search_size,
        region_growing_noise_level,
        region_growing_angle,
        region_growing_min_area);

      clean_roof_regions_3(min_size);

      create_approximate_roofs();
    }

    void compute_roofs(
      const std::size_t kinetic_max_intersections,
      const FT graph_cut_beta_3) {

      if (m_data.verbose) 
        std::cout << std::endl << "- Computing building roofs" 
        << std::endl;

      create_partitioning_input_3();

      create_partitioning_output_3(kinetic_max_intersections);

      create_building_bounds(graph_cut_beta_3);

      create_exact_roofs_and_walls();
    }

    // FLAGS.

    const bool has_exact_boundaries() const {
      return m_has_exact_boundaries;
    }

    const bool has_exact_roofs() const {
      return m_has_exact_roofs;
    }

    // OUTPUT

    template<typename OutputIterator>
    void return_boundary_points(OutputIterator output) const {

      const auto& points = m_data.building_boundary_points_2;
      const auto& plane = m_data.planar_ground.plane;

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
      const auto& plane = m_data.planar_ground.plane;

      std::vector< std::pair<Point_3, long> > data(points.size());
      for (std::size_t i = 0; i < points.size(); ++i)
        data[i] = std::make_pair(
                  internal::position_on_plane_3(
                    points[i], 
                    plane), -1);

        for(std::size_t i = 0; i < indices.size(); ++i)
          for(std::size_t j = 0; j < indices[i].size(); ++j)
            data[indices[i][j]].second = static_cast<long>(i);

        std::copy(data.begin(), data.end(), output);
    }

    template<typename OutputIterator>
    void return_approximate_boundary_edges(OutputIterator output) const {

      const auto& points = m_data.building_boundary_points_2;
      const auto& indices = m_data.building_boundary_indices_2;
      const auto& plane = m_data.planar_ground.plane;

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
      const auto& plane = m_data.planar_ground.plane;

      for (std::size_t i = 0; i < buildings.size(); ++i) {
        const auto& edges = buildings[i].edges;

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
      const auto& plane = m_data.planar_ground.plane;

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
      FacesOutputIterator output_faces,
      const bool extruded) const {

      const auto& buildings = m_data.buildings;
      const auto& plane = m_data.planar_ground.plane;
      
      internal::Indexer<Point_2> indexer;
      std::size_t num_vertices = 0;

      for (std::size_t i = 0; i < buildings.size(); ++i) {
        const auto& triangles = buildings[i].triangles;

        for (std::size_t j = 0; j < triangles.size(); ++j) {
          cpp11::array<std::size_t, 3> face;
          
          for (std::size_t k = 0; k < 3; ++k) {
            const auto& point = triangles[j][k];

            const std::size_t idx = indexer(point);
            if (idx == num_vertices) {

              if (!extruded)
                *(output_vertices++) = 
                internal::position_on_plane_3(point, plane);
              else 
                *(output_vertices++) = 
                Point_3(point.x(), point.y(), buildings[i].height);

              ++num_vertices;
            }
            face[k] = idx;
          }
          *(output_faces++) = std::make_pair(face, i);
        }
      }
    }

    template<typename OutputIterator>
    void return_roof_points(OutputIterator output) const {

      const auto& clusters = m_data.building_clusters;
      const auto& buildings = m_data.buildings;

      for (std::size_t i = 0; i < buildings.size(); ++i) {  
        const auto& cluster = clusters[buildings[i].cluster_index];

        std::vector< std::tuple<Point_3, long, long> > data(cluster.size());
        for (std::size_t j = 0; j < cluster.size(); ++j)
          data[j] = std::make_tuple(
            get(m_data.point_map, *(cluster[j])), 
            static_cast<long>(i), 
            -1);

        const auto& indices = buildings[i].roof_indices;
        for (std::size_t j = 0; j < indices.size(); ++j)
          for (std::size_t k = 0; k < indices[j].size(); ++k)
            std::get<2>(data[indices[j][k]]) = static_cast<long>(j);
        
        for (std::size_t j = 0; j < data.size(); ++j)
          *(output++) = data[j];
      }
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void return_approximate_roofs(
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces) const {

      const auto& buildings = m_data.buildings;
      
      internal::Indexer<Point_3> indexer;
      std::size_t num_vertices = 0;

      for (std::size_t i = 0; i < buildings.size(); ++i) {
        const auto& roofs = buildings[i].approximate_roofs;

        for (std::size_t j = 0; j < roofs.size(); ++j) {
          std::vector<std::size_t> face(roofs[j].size());
          
          for (std::size_t k = 0; k < roofs[j].size(); ++k) {
            const auto& point = roofs[j][k];

            const std::size_t idx = indexer(point);
            if (idx == num_vertices) {

              *(output_vertices++) = point;
              ++num_vertices;
            }
            face[k] = idx;
          }
          *(output_faces++) = std::make_tuple(face, i, j);
        }
      }
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void return_exact_roofs(
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces) const {

      const auto& buildings = m_data.buildings;
      
      internal::Indexer<Point_3> indexer;
      std::size_t num_vertices = 0;

      for (std::size_t i = 0; i < buildings.size(); ++i) {
        const auto& roofs = buildings[i].roofs;

        for (std::size_t j = 0; j < roofs.size(); ++j) {
          std::vector<std::size_t> face(roofs[j].vertices.size());
          
          for (std::size_t k = 0; k < roofs[j].vertices.size(); ++k) {
            const auto& point = roofs[j].vertices[k];

            const std::size_t idx = indexer(point);
            if (idx == num_vertices) {

              *(output_vertices++) = point;
              ++num_vertices;
            }
            face[k] = idx;
          }
          *(output_faces++) = std::make_tuple(face, i, j);
        }
      }
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void return_partitioning_input_3(
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces) const {

      const auto& buildings = m_data.buildings;
      
      internal::Indexer<Point_3> indexer;
      std::size_t num_vertices = 0;

      for (std::size_t i = 0; i < buildings.size(); ++i) {
        
        const auto& roofs = buildings[i].approximate_roofs;
        const auto& walls = buildings[i].approximate_walls;
        const auto& ground = buildings[i].approximate_ground;

        // Roofs.
        for (std::size_t j = 0; j < roofs.size(); ++j) {
          std::vector<std::size_t> face(roofs[j].size());
          
          for (std::size_t k = 0; k < roofs[j].size(); ++k) {
            const auto& point = roofs[j][k];

            const std::size_t idx = indexer(point);
            if (idx == num_vertices) {

              *(output_vertices++) = point;
              ++num_vertices;
            }
            face[k] = idx;
          }
          *(output_faces++) = std::make_pair(face, i);
        }

        // Walls.
        for (std::size_t j = 0; j < walls.size(); ++j) {
          std::vector<std::size_t> face(walls[j].size());
          
          for (std::size_t k = 0; k < walls[j].size(); ++k) {
            const auto& point = walls[j][k];

            const std::size_t idx = indexer(point);
            if (idx == num_vertices) {

              *(output_vertices++) = point;
              ++num_vertices;
            }
            face[k] = idx;
          }
          *(output_faces++) = std::make_pair(face, i);
        }

        // Ground.
        std::vector<std::size_t> face(ground.size());
        for (std::size_t k = 0; k < ground.size(); ++k) {
          const auto& point = ground[k];

          const std::size_t idx = indexer(point);
          if (idx == num_vertices) {

            *(output_vertices++) = point;
            ++num_vertices;
          }
          face[k] = idx;
        }
        if (!ground.empty())
          *(output_faces++) = std::make_pair(face, i);
      }
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void return_partitioning_output_3(
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces,
      const bool with_visibility) const {

      const auto& buildings = m_data.buildings;
      internal::Indexer<Point_3> indexer;
      
      std::size_t num_vertices = 0;
      for (std::size_t i = 0; i < buildings.size(); ++i) {
        
        const auto& polyhedrons = buildings[i].polyhedrons;
        for (std::size_t j = 0; j < polyhedrons.size(); ++j) {

          if (with_visibility &&
            polyhedrons[j].visibility == Visibility_label::OUTSIDE)
            continue;
          
          const auto& faces = polyhedrons[j].faces;
          const auto& vertices = polyhedrons[j].vertices;
          
          for (std::size_t k = 0; k < faces.size(); ++k) {
            
            std::vector<std::size_t> face(faces[k].size());
            for (std::size_t l = 0; l < faces[k].size(); ++l) {
              
              const auto& point = vertices[faces[k][l]];

              const std::size_t idx = indexer(point);
              if (idx == num_vertices) {

                *(output_vertices++) = point;
                ++num_vertices;
              }
              face[l] = idx;
            }
            *(output_faces++) = std::make_pair(face, i);
          }
        }
      } 
    }

  private:
    Data_structure& m_data;
    
    bool m_has_exact_boundaries;
    bool m_has_exact_roofs;

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
      std::sort(indices.begin(), indices.end(), estimator.sorter());

      Points_region_growing_2 region_growing(
        indices,
        connectivity,
        conditions);

      region_growing.detect(
        m_data.building_boundary_indices_2);

      if (m_data.verbose) 
        std::cout << "-> " << m_data.building_boundary_indices_2.size()
        << " wall(s) detected" 
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
    void compute_building_footprints_2() {

      if (m_data.verbose) 
        std::cout << "* computing footprints" 
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
        << " building footprint(s) computed" 
        << std::endl;
    }

    // Final buildings.
    void finilize_buildings(const std::size_t min_faces_per_building) {

      m_data.buildings.clear();

      const auto& faces = m_data.building_polygon_faces_2;
      const auto& footprints = m_data.building_footprints_2;

      Building building;
      auto& triangles = building.triangles;
      auto& segments = building.edges;

      Building_footprints_2 footprints_extractor;
      Building_boundaries_2 boundaries_extractor;

      std::size_t building_count = 0;
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
        if (triangles.size() >= min_faces_per_building && segments.size() >= 3) {

          building.cluster_index = building_count;
          m_data.buildings.push_back(building);
          ++building_count;
        }
      }
      m_has_exact_boundaries = true;
    }

    // Clustering.
    void cluster_building_points() {

      if (m_data.verbose) 
        std::cout << "* clustering building points"
        << std::endl;

      m_data.building_clusters.clear();
      const auto& points = m_data.building_interior_points();

      const Buildings_clustering clustering(
        points, 
        m_data.point_map);
      
      clustering.create(
        m_data.buildings,
        m_data.building_clusters);
    }

    // Compute building heights.
    void compute_building_heights(const Extrusion_type extrusion_type) {

      if (m_data.verbose) 
        std::cout << "* computing heights"
        << std::endl;

      const auto& clusters = m_data.building_clusters;
      
      const Building_height_estimator estimator(clusters, m_data.point_map);
      estimator.compute_heights(extrusion_type, m_data.buildings);
    }

    // Roofs.
    void extract_roof_regions_3(
      const FT region_growing_search_size,
      const FT region_growing_noise_level,
      const FT region_growing_angle,
      const FT region_growing_min_area) {

      if (m_data.verbose) 
        std::cout << "* region growing" 
        << std::endl;

      auto& buildings = m_data.buildings;
      const auto& clusters = m_data.building_clusters;

      std::size_t num_regions = 0;
      for (std::size_t i = 0; i < buildings.size(); ++i) {
        
        Building& building = buildings[i];
        const Cluster& cluster = clusters[building.cluster_index];

        num_regions += apply_region_growing_3(
          region_growing_search_size,
          region_growing_noise_level,
          region_growing_angle,
          region_growing_min_area,
          cluster,
          building);
      }

      if (m_data.verbose)
        std::cout << "-> " << num_regions
        << " regions(s) found" 
        << std::endl;
    }

    std::size_t apply_region_growing_3(
      const FT region_growing_search_size,
      const FT region_growing_noise_level,
      const FT region_growing_angle,
      const FT region_growing_min_area,
      const Cluster& cluster,
      Building& building) {

      building.roof_indices.clear();

      Points_connectivity_3 connectivity(
        cluster,
        m_data.point_map, 
        region_growing_search_size);

      Normals_estimator_3 estimator(
        cluster, 
        m_data.point_map,
        connectivity);

      Points_conditions_3 conditions(
        cluster,
        m_data.point_map, 
        estimator.normals(),
        region_growing_noise_level,
        region_growing_angle,
        region_growing_min_area);

      std::vector<std::size_t> indices(cluster.size());
      for (std::size_t i = 0; i < cluster.size(); ++i)
        indices[i] = i;
      std::sort(indices.begin(), indices.end(), estimator.sorter());

      Points_region_growing_3 region_growing(
        indices,
        connectivity,
        conditions);

      region_growing.detect(building.roof_indices);
      building.normals = estimator.normals();

      return building.roof_indices.size();
    }

    void clean_roof_regions_3(const FT scale) {

      if (m_data.verbose) 
        std::cout << "* cleaning" 
        << std::endl;

      auto& buildings = m_data.buildings;
      const auto& clusters = m_data.building_clusters;

      std::size_t num_roofs = 0;
      for (std::size_t i = 0; i < buildings.size(); ++i) {
        
        Building& building = buildings[i];
        const Cluster& cluster = clusters[building.cluster_index];

        const Roof_cleaner cleaner(
          cluster, 
          m_data.point_map,
          building.normals,
          scale);

        cleaner.clean(building.roof_indices);
        num_roofs += building.roof_indices.size();
      }

      if (m_data.verbose)
        std::cout << "-> " << num_roofs
        << " roof(s) detected" 
        << std::endl;
    }

    void create_approximate_roofs() {

      if (m_data.verbose) 
        std::cout << "* creating approximate roofs" 
        << std::endl;

      auto& buildings = m_data.buildings;
      const auto& clusters = m_data.building_clusters;

      for (std::size_t i = 0; i < buildings.size(); ++i) {
        
        const Roof_estimator estimator(
          clusters[buildings[i].cluster_index], 
          m_data.point_map,
          buildings[i].roof_indices);

        estimator.estimate(buildings[i].approximate_roofs);
      }
    }

    void create_partitioning_input_3() {

      if (m_data.verbose) 
        std::cout << "* creating partitioning input" 
        << std::endl;

      auto& buildings = m_data.buildings;
      for (std::size_t i = 0; i < buildings.size(); ++i) {
        
        const Wall_estimator westimator(
          buildings[i].edges, m_data.planar_ground.plane, buildings[i].height);
        westimator.estimate(buildings[i].approximate_walls);

        const Building_ground_estimator gestimator(
          buildings[i].triangles, m_data.planar_ground.plane);
        gestimator.estimate(buildings[i].approximate_ground);
      }
    }

    void create_partitioning_output_3(
      const std::size_t kinetic_max_intersections) {

      if (m_data.verbose) 
        std::cout << "* creating partitioning output" 
        << std::endl;

      std::size_t num_polyhedrons = 0;
      auto& buildings = m_data.buildings;

      for (std::size_t i = 0; i < buildings.size(); ++i) {
        auto& building = buildings[i];

        const Kinetic_partitioning_3 kinetic(
          building.approximate_walls,
          building.approximate_roofs,
          building.approximate_ground,
          kinetic_max_intersections);

        kinetic.compute(building.polyhedrons);
        num_polyhedrons += building.polyhedrons.size();
      }

      if (m_data.verbose)
        std::cout << "-> " << num_polyhedrons
        << " polyhedron facet(s) created" 
        << std::endl;
    }

    void create_building_bounds(const FT graph_cut_beta_3) {

      compute_visibility_3();
      apply_graph_cut_3(graph_cut_beta_3);
    }

    void compute_visibility_3() {

      if (m_data.verbose) 
        std::cout << "* computing visibility" 
        << std::endl;

      std::size_t num_polyhedrons = 0;
      auto& buildings = m_data.buildings;

      for (std::size_t i = 0; i < buildings.size(); ++i) {
        auto& building = buildings[i];

        const Visibility_3 visibility(
          m_data.building_clusters[building.cluster_index],
          m_data.point_map,
          building);

        visibility.compute(building.polyhedrons);

        for (std::size_t j = 0; j < building.polyhedrons.size(); ++j)
          if (building.polyhedrons[j].visibility == Visibility_label::INSIDE)
            ++num_polyhedrons;
      }

      if (m_data.verbose)
        std::cout << "-> " << num_polyhedrons
        << " polyhedron facet(s) after visibility" 
        << std::endl;
    }

    void apply_graph_cut_3(const FT graph_cut_beta_3) {

      if (m_data.verbose) 
        std::cout << "* applying graphcut" 
        << std::endl;

      std::size_t num_polyhedrons = 0;
      auto& buildings = m_data.buildings;

      for (std::size_t i = 0; i < buildings.size(); ++i) {
        auto& building = buildings[i];

        const Weight_quality_estimator_3 estimator(building.polyhedrons);
        estimator.estimate(building.graphcut_faces);

        const Graphcut_3 graphcut(graph_cut_beta_3);
        graphcut.apply(
          building.graphcut_faces,
          building.polyhedrons);

        for (auto& polyhedron : building.polyhedrons)
          if (polyhedron.visibility == Visibility_label::INSIDE)
            ++num_polyhedrons;
      }

      if (m_data.verbose)
        std::cout << "-> " << num_polyhedrons
        << " polyhedron facet(s) after graphcut" 
        << std::endl;
    }

    void create_exact_roofs_and_walls() {

      m_has_exact_roofs = true;

      if (m_data.verbose) 
        std::cout << "* extracting final roofs" 
        << std::endl;

      std::size_t num_roofs = 0;
      for (auto& building : m_data.buildings) {

        const Roof_wall_extractor extractor(
          building.polyhedrons, m_data.planar_ground.plane);

        extractor.extract(building.roofs, building.walls);
        num_roofs += building.roofs.size();
      }

      if (m_data.verbose)
        std::cout << "-> " << num_roofs
        << " final roof(s) extracted" 
        << std::endl;
    }

  }; // Buildings

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_BUILDINGS_H
