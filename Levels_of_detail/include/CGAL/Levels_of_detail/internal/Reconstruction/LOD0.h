#ifndef CGAL_LEVELS_OF_DETAIL_0_H
#define CGAL_LEVELS_OF_DETAIL_0_H

// Internal includes.
#include <CGAL/Levels_of_detail/internal/structures.h>
#include <CGAL/Levels_of_detail/internal/Ground/Planar_ground_estimator.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename DataStructure>
  class LOD0 {

  public:
    using Data_structure = DataStructure;
    using Traits = typename Data_structure::Traits;
    using Planar_ground_estimator = Planar_ground_estimator<Traits>;

    LOD0(Data_structure& data_structure) :
    m_data(data_structure),
    m_planar_ground_estimator(m_data.planar_ground.plane)
    { }

    void reconstruct() {

      const auto& ground = m_data.planar_ground;
      m_planar_ground_estimator.initialize(ground.bounding_box);

      auto& buildings = m_data.buildings;
      auto& trees = m_data.trees;

      set_object_indices(buildings);
      set_object_indices(trees);

      for (const auto& building : buildings)
        m_planar_ground_estimator.add_urban_object(building);
      for (const auto& tree : trees)
        m_planar_ground_estimator.add_urban_object(tree);

      for (const auto& building : buildings)
        m_planar_ground_estimator.tag_faces(building);
      for (const auto& tree : trees)
        m_planar_ground_estimator.tag_faces(tree);

      m_planar_ground_estimator.finilize();
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void output_as_triangle_soup(
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces) {

      m_planar_ground_estimator.output_as_triangle_soup(
        output_vertices, output_faces);
    }

  private:
    Data_structure& m_data;
    Planar_ground_estimator m_planar_ground_estimator;

    template<typename Urban_object>
    void set_object_indices(std::vector<Urban_object>& objects) const {
      for (std::size_t i = 0; i < objects.size(); ++i)
        objects[i].object_index = i;
    }

  }; // LOD0

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_0_H
