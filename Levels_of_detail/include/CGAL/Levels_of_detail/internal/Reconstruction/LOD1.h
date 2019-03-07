#ifndef CGAL_LEVELS_OF_DETAIL_1_H
#define CGAL_LEVELS_OF_DETAIL_1_H

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utilities.h>
#include <CGAL/Levels_of_detail/internal/structures.h>
#include <CGAL/Levels_of_detail/internal/Ground/Smooth_ground_estimator.h>
#include <CGAL/Levels_of_detail/internal/Reconstruction/Knn_search.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename DataStructure>
  class LOD1 {

  public:
    using Data_structure = DataStructure;
    
    using Traits = typename Data_structure::Traits;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Segment_2 = typename Traits::Segment_2;
    using Plane_3 = typename Traits::Plane_3;
    using Triangle_2 = typename Traits::Triangle_2;
    using Triangle_3 = typename Traits::Triangle_3;

    using Filtered_range = typename Data_structure::Filtered_range;
    using Filter_iterator = typename Data_structure::Filter_iterator;
    using Point_map = typename Data_structure::Point_map;

    using Point_2_from_iterator = 
    Point_2_from_iterator_map<Filter_iterator, Point_2, Point_map>;
    using Knn_search = 
    Knn_search<Traits, Filter_iterator, Point_2_from_iterator>;
    using Smooth_ground_estimator = 
    Smooth_ground_estimator<Traits, Knn_search, Filtered_range, Point_map>;

    using Building = typename Data_structure::Building;

    LOD1(
      Data_structure& data_structure,
      const FT ground_precision) :
    m_data(data_structure),
    m_ground_plane(m_data.planar_ground.plane),
    m_ground_points(m_data.ground_points()),
    m_point_map(m_data.point_map),
    m_knn_search(
      m_ground_points,
      Point_2_from_iterator(m_point_map),
      6),
    m_smooth_ground_estimator(
      m_ground_plane,
      m_knn_search,
      m_ground_points,
      m_point_map,
      ground_precision)
    { }

    void reconstruct() {

      const auto& ground = m_data.planar_ground;
      m_smooth_ground_estimator.initialize(ground.bounding_box);

      auto& buildings = m_data.buildings;
      auto& trees = m_data.trees;

      set_object_indices(buildings);
      set_object_indices(trees);

      for (const auto& building : buildings)
        m_smooth_ground_estimator.add_urban_object(building);
      for (const auto& tree : trees)
        m_smooth_ground_estimator.add_urban_object(tree);

      m_smooth_ground_estimator.finilize();
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void output_as_triangle_soup(
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces) {

      std::size_t num_vertices = 0;
      internal::Indexer<Point_3> indexer;

      m_smooth_ground_estimator.output_as_triangle_soup(
        output_vertices, output_faces, num_vertices, indexer);

      // Buildings.
      add_urban_objects(
        m_data.buildings, Urban_object_type::BUILDING,
        output_vertices, output_faces,
        num_vertices, indexer);

      // Trees.
      add_urban_objects(
        m_data.trees, Urban_object_type::TREE,
        output_vertices, output_faces,
        num_vertices, indexer);
    }

  private:
    Data_structure& m_data;

    const Plane_3& m_ground_plane;
    const Filtered_range m_ground_points;
    const Point_map m_point_map;

    Knn_search m_knn_search;
    Smooth_ground_estimator m_smooth_ground_estimator;

    template<typename Urban_object>
    void set_object_indices(std::vector<Urban_object>& objects) const {
      for (std::size_t i = 0; i < objects.size(); ++i)
        objects[i].object_index = i;
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator,
    typename Urban_object>
    void add_urban_objects(
      const std::vector<Urban_object>& objects,
      const Urban_object_type urban_tag,
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces,
      std::size_t& num_vertices,
      internal::Indexer<Point_3>& indexer) const {

      for (const auto& object : objects) {
        add_walls(
          object, urban_tag,
          output_vertices, output_faces,
          num_vertices, indexer);
        add_roofs(
          object, urban_tag,
          output_vertices, output_faces,
          num_vertices, indexer);
      }
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator,
    typename Urban_object>
    void add_walls(
      const Urban_object& object,
      const Urban_object_type urban_tag,
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces,
      std::size_t& num_vertices,
      internal::Indexer<Point_3>& indexer) const {

      for (const Segment_2& segment : object.boundaries()) {
        const Point_2& s = segment.source();
        const Point_2& t = segment.target();

        const FT z1 = m_smooth_ground_estimator.get_z(s);
        const FT z2 = m_smooth_ground_estimator.get_z(t);
        const FT ztop = object.height;

        const Point_3 p1 = Point_3(s.x(), s.y(), z1);
        const Point_3 p2 = Point_3(t.x(), t.y(), z2);
        const Point_3 p3 = Point_3(t.x(), t.y(), ztop);
        const Point_3 p4 = Point_3(s.x(), s.y(), ztop);

        const Triangle_3 tri1(p1, p2, p3);
        const Triangle_3 tri2(p3, p4, p1);

        add_face(
          tri1, urban_tag,
          output_vertices, output_faces, 
          num_vertices, indexer);
        add_face(
          tri2, urban_tag,
          output_vertices, output_faces, 
          num_vertices, indexer);
      }
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator,
    typename Urban_object>
    void add_roofs(
      const Urban_object& object,
      const Urban_object_type urban_tag,
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces,
      std::size_t& num_vertices,
      internal::Indexer<Point_3>& indexer) const {

      for (const Triangle_2& triangle : object.footprint()) { 
        const Point_2& a = triangle[0];
        const Point_2& b = triangle[1];
        const Point_2& c = triangle[2];

        const FT ztop = object.height;

        const Point_3 p1 = Point_3(a.x(), a.y(), ztop);
        const Point_3 p2 = Point_3(b.x(), b.y(), ztop);
        const Point_3 p3 = Point_3(c.x(), c.y(), ztop);

        const Triangle_3 tri(p1, p2, p3);
        add_face(
          tri, urban_tag,
          output_vertices, output_faces, 
          num_vertices, indexer);
      }
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void add_face(
      const Triangle_3& triangle,
      const Urban_object_type urban_tag,
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces,
      std::size_t& num_vertices,
      internal::Indexer<Point_3>& indexer) const {

      cpp11::array<std::size_t, 3> face;
      for (std::size_t k = 0; k < 3; ++k) {
        const Point_3& p = triangle[k];

        const std::size_t idx = indexer(p);
        if (idx == num_vertices) {

          *(output_vertices++) = p;
          ++num_vertices;
        }
        face[k] = idx;
      }
      *(output_faces++) = 
        std::make_pair(face, urban_tag);
    }

  }; // LOD1

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_1_H
