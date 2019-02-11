#ifndef CGAL_LEVELS_OF_DETAIL_VEGETATION_H
#define CGAL_LEVELS_OF_DETAIL_VEGETATION_H

// STL includes.
#include <string>
#include <vector>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utilities.h>

// Vegetation.
#include <CGAL/Levels_of_detail/internal/Vegetation/Trees.h>
#include <CGAL/Levels_of_detail/internal/Vegetation/Tree_footprints_2.h>
#include <CGAL/Levels_of_detail/internal/Vegetation/Tree_boundaries_2.h>
#include <CGAL/Levels_of_detail/internal/Vegetation/Vegetation_clustering.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename DataStructure>
  class Vegetation {

  public:
    using Data_structure = DataStructure;

    using Traits = typename DataStructure::Traits;
    using Filtered_range = typename Data_structure::Filtered_range;
    using Point_map = typename Data_structure::Point_map;

    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;

    using Tree = typename Data_structure::Tree;

    using Vegetation_clustering = 
    Vegetation_clustering<Traits, Filtered_range, Point_map>;
    
    using Trees = 
    Trees<Traits, Filtered_range, Point_map>;

    using Tree_footprints_2 = Tree_footprints_2<Traits>;
    using Tree_boundaries_2 = Tree_boundaries_2<Traits>;

    Vegetation(Data_structure& data_structure) :
    m_data(data_structure)
    { }
    
    // PROCESSING

    void detect_footprints(
      const FT grid_cell_width, 
      const FT min_height, 
      const FT min_radius,
      const std::size_t min_faces_per_tree) {

      if (m_data.verbose) 
        std::cout << std::endl << "- Detecting tree footprints" 
        << std::endl;

      cluster_vegetation_points(
        grid_cell_width, 
        min_height);

      estimate_trees(
        min_radius);

      finilize_trees(
        min_faces_per_tree);
    }

    void extrude_footprints(const Extrusion_type extrusion_type) {
      
      if (m_data.verbose) 
        std::cout << std::endl << "- Extruding tree footprints" 
        << std::endl;

      compute_tree_heights(
        extrusion_type);
    }

    // OUTPUT

    template<typename OutputIterator>
    void return_clustered_points(OutputIterator output) const {

      const auto& points = m_data.vegetation_points();
      const auto& clusters = m_data.vegetation_clusters;

      for (std::size_t i = 0; i < clusters.size(); ++i) {
        for (std::size_t j = 0; j < clusters[i].size(); ++j) {

          const auto& it = clusters[i][j];
          const auto& point = get(m_data.point_map, *it);

          *(output++) = std::make_pair(point, i);
        }
      }
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void return_footprints(
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces,
      const bool extruded = false) const {

      const auto& trees = m_data.trees;
      const auto& plane = m_data.ground_plane;
      
      internal::Indexer<Point_2> indexer;
      std::size_t num_vertices = 0;

      for (std::size_t i = 0; i < trees.size(); ++i) {
        const auto& triangles = trees[i].footprint;

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
                Point_3(point.x(), point.y(), trees[i].height);

              ++num_vertices;
            }
            face[k] = idx;
          }
          *(output_faces++) = std::make_pair(face, i);
        }
      }
    }

    template<typename OutputIterator>
    void return_boundary_edges(OutputIterator output) const {

      const auto& trees = m_data.trees;
      const auto& plane = m_data.ground_plane;

      for (std::size_t i = 0; i < trees.size(); ++i) {
        const auto& edges = trees[i].boundaries;

        for (std::size_t j = 0; j < edges.size(); ++j)
          *(output++) = 
          internal::segment_3_from_segment_2_and_plane(edges[j], plane);
      }
    }

  private:
    Data_structure& m_data;

    // Clustering.
    void cluster_vegetation_points(
      const FT grid_cell_width, 
      const FT min_height) {

      if (m_data.verbose) 
        std::cout << "* clustering vegetation points" 
        << std::endl;

      m_data.vegetation_clusters.clear();
      const auto& points = m_data.vegetation_points();

      Vegetation_clustering clustering(
        points, 
        m_data.point_map);
      
      clustering.detect(
        grid_cell_width, 
        min_height,
        m_data.vegetation_clusters);

      if (m_data.verbose) 
        std::cout << "-> " << m_data.vegetation_clusters.size()
        << " cluster(s) found"
        << std::endl;
    }

    // Estimating trees.
    void estimate_trees(const FT min_radius) {

      if (m_data.verbose) 
        std::cout << "* detecting footprints" 
        << std::endl;

      m_data.trees.clear();
      const auto& clusters = m_data.vegetation_clusters;

      const Trees trees(clusters, m_data.point_map);
      trees.estimate(min_radius, m_data.trees);

      if (m_data.verbose) 
        std::cout << "-> " << m_data.trees.size()
        << " trees footprint(s) detected"
        << std::endl;
    }

    // Final trees.
    void finilize_trees(const std::size_t min_faces_per_tree) {

      auto& trees = m_data.trees;
      Tree_footprints_2 footprints_extractor;
      Tree_boundaries_2 boundaries_extractor;

      for (std::size_t i = 0; i < trees.size(); ++i) {  
        auto& tree = trees[i];

        footprints_extractor.create_footprint_triangles(
          tree.center, 
          tree.radius,
          min_faces_per_tree,
          tree.footprint);

        boundaries_extractor.create_boundary_segments(
          tree.center,
          tree.radius,
          min_faces_per_tree,
          tree.boundaries);
      }
    }

    // Compute tree heights.
    void compute_tree_heights(const Extrusion_type extrusion_type) {

      if (m_data.verbose) 
        std::cout << "* computing heights"
        << std::endl;

      const auto& clusters = m_data.vegetation_clusters;
      
      const Trees trees(clusters, m_data.point_map);
      trees.compute_heights(extrusion_type, m_data.trees);
    }

  }; // Vegetation

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_VEGETATION_H
