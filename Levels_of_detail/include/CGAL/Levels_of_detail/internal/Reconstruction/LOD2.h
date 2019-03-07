#ifndef CGAL_LEVELS_OF_DETAIL_2_H
#define CGAL_LEVELS_OF_DETAIL_2_H

// Internal includes.
#include <CGAL/Levels_of_detail/internal/utilities.h>
#include <CGAL/Levels_of_detail/internal/structures.h>
#include <CGAL/Levels_of_detail/internal/Ground/Smooth_ground_estimator.h>
#include <CGAL/Levels_of_detail/internal/Reconstruction/Knn_search.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<typename DataStructure>
  class LOD2 {

  public:
    using Data_structure = DataStructure;
    
    using Traits = typename Data_structure::Traits;
    using FT = typename Traits::FT;
    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Segment_3 = typename Traits::Segment_3;
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
    using Tree = typename Data_structure::Tree;

    LOD2(
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
      add_buildings(
        m_data.buildings,
        output_vertices, output_faces,
        num_vertices, indexer);

      // Trees.
      add_trees(
        m_data.trees,
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
    typename FacesOutputIterator>
    void add_buildings(
      const std::vector<Building>& buildings,
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces,
      std::size_t& num_vertices,
      internal::Indexer<Point_3>& indexer) const {

      for (const auto& building : buildings) {
        add_walls(
          building,
          output_vertices, output_faces,
          num_vertices, indexer);
        add_roofs(
          building, 
          output_vertices, output_faces,
          num_vertices, indexer);
      }
    }

    bool are_equal(const Point_3& p, const Point_3& q) const {

      const FT tolerance = FT(1) / FT(1000000);
      return 
        CGAL::abs(p.x() - q.x()) < tolerance &&
        CGAL::abs(p.y() - q.y()) < tolerance &&
        CGAL::abs(p.z() - q.z()) < tolerance;
    }

    bool are_equal(const Segment_3& s1, const Segment_3& s2) const {

      return 
        (are_equal(s1.source(), s2.source()) && 
        are_equal(s1.target(), s2.target())) ||
        (are_equal(s1.source(), s2.target()) && 
        are_equal(s1.target(), s2.source()));
    }

    bool is_boundary(
      const Segment_3& query, 
      const Building& building) const {

      const auto& roofs = building.roofs;

      std::size_t count = 0;
      for (const auto& roof : roofs) {
        const auto& vertices = roof.vertices;

        for (std::size_t i = 0; i < vertices.size(); ++i) {
          const std::size_t ip = (i + 1) % vertices.size();
          const Segment_3 segment = Segment_3(vertices[i], vertices[ip]);

          if (are_equal(query, segment))
            ++count;

          if (count > 1) 
            return false;
        }
      }
      return true;
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void add_walls(
      const Building& building,
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces,
      std::size_t& num_vertices,
      internal::Indexer<Point_3>& indexer) const {

      std::vector<Segment_3> segments;
      for (const auto& roof : building.roofs) {

        const auto& vertices = roof.vertices;
        for (std::size_t i = 0; i < vertices.size(); ++i) {
          const std::size_t ip = (i + 1) % vertices.size();

          const Segment_3 segment = Segment_3(vertices[i], vertices[ip]);
          if (is_boundary(segment, building))
            segments.push_back(segment);
        }
      }

      for (const Segment_3& segment : segments) {

        const Point_3& s = segment.source();
        const Point_3& t = segment.target();

        const FT z3 = m_smooth_ground_estimator.get_z(s);
        const FT z4 = m_smooth_ground_estimator.get_z(t);

        const Point_3 p1 = s;
        const Point_3 p2 = t;
        const Point_3 p3 = Point_3(t.x(), t.y(), z4);
        const Point_3 p4 = Point_3(s.x(), s.y(), z3);

        const Triangle_3 tri1(p1, p2, p3);
        const Triangle_3 tri2(p3, p4, p1);

        add_face(
          tri1, Urban_object_type::BUILDING,
          output_vertices, output_faces, 
          num_vertices, indexer);
        add_face(
          tri2, Urban_object_type::BUILDING,
          output_vertices, output_faces, 
          num_vertices, indexer);
      }
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void add_roofs(
      const Building& building,
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces,
      std::size_t& num_vertices,
      internal::Indexer<Point_3>& indexer) const {

      for (const auto& roof : building.roofs) {

        const auto& vertices = roof.vertices;
        const Point_3& p = vertices[0];

        for (std::size_t i = 1; i < vertices.size() - 1; ++i) {
          const std::size_t ip = i + 1;
          const Triangle_3 tri = Triangle_3(p, vertices[i], vertices[ip]);
          add_face(
            tri, Urban_object_type::BUILDING,
            output_vertices, output_faces, 
            num_vertices, indexer);
        }
      }
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void add_trees(
      const std::vector<Tree>& trees,
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces,
      std::size_t& num_vertices,
      internal::Indexer<Point_3>& indexer) const {

      for (const auto& tree : trees) {
        add_tree(
          tree, 
          output_vertices, output_faces,
          num_vertices, indexer);
      }
    }

    template<
    typename VerticesOutputIterator,
    typename FacesOutputIterator>
    void add_tree(
      const Tree& tree,
      VerticesOutputIterator output_vertices,
      FacesOutputIterator output_faces,
      std::size_t& num_vertices,
      internal::Indexer<Point_3>& indexer) const {
        
      const auto& vertices = tree.vertices;
      const auto& faces = tree.faces;
      const auto& bottom = tree.bottom;

      for (std::size_t i = 0; i < faces.size(); ++i) {
        const Point_3& p1 = vertices[faces[i][0]];
        const Point_3& p2 = vertices[faces[i][1]];
        const Point_3& p3 = vertices[faces[i][2]];

        Point_3 a, b, c;
        if (!bottom[faces[i][0]]) a = p1;
        else a = Point_3(p1.x(), p1.y(), m_smooth_ground_estimator.get_z(p1));
        if (!bottom[faces[i][1]]) b = p2;
        else b = Point_3(p2.x(), p2.y(), m_smooth_ground_estimator.get_z(p2));
        if (!bottom[faces[i][2]]) c = p3;
        else c = Point_3(p3.x(), p3.y(), m_smooth_ground_estimator.get_z(p3));

        const Triangle_3 tri = Triangle_3(a, b, c);
        add_face(
          tri, Urban_object_type::TREE,
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

  }; // LOD2

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_2_H
