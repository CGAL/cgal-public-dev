#ifndef CGAL_LEVELS_OF_DETAIL_VEGETATION_H
#define CGAL_LEVELS_OF_DETAIL_VEGETATION_H

// STL includes.
#include <string>
#include <vector>

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utilities.h>

// Vegetation.
// #include <CGAL/Levels_of_detail/internal/Vegetation/Trees.h>
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
    using Tree = typename Data_structure::Tree;

    using Vegetation_clustering = 
    Vegetation_clustering<Traits, Filtered_range, Point_map>;

    // using Trees = Trees<Traits>;

    Vegetation(Data_structure& data_structure) :
    m_data(data_structure)
    { }
    
    // PROCESSING

    void detect_footprints(
      const FT grid_cell_width, 
      const FT min_height, 
      const FT min_radius) {

      if (m_data.verbose) 
        std::cout << std::endl << "- Detecting tree footprints" 
        << std::endl;

      cluster_vegetation_points(
        grid_cell_width, 
        min_height);

      estimate_trees(min_radius);
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
      FacesOutputIterator output_faces) const {

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
        std::cout << "* estimating trees" 
        << std::endl;

      m_data.trees.clear();

      const auto& clusters = m_data.vegetation_clusters;

      // Trees trees(clusters);
      // trees.estimate(min_radius, m_data.trees);

      if (m_data.verbose) 
        std::cout << "-> " << m_data.trees.size()
        << " trees(s) estimated"
        << std::endl;
    }

  }; // Vegetation

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_VEGETATION_H
