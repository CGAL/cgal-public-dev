#ifndef CGAL_LEVELS_OF_DETAIL_DATA_STRUCTURE_H
#define CGAL_LEVELS_OF_DETAIL_DATA_STRUCTURE_H

// STL includes.
#include <vector>

// Boost includes.
#include <boost/iterator/filter_iterator.hpp>

// CGAL includes.
#include <CGAL/Iterator_range.h>

// Internal includes.
#include <CGAL/Levels_of_detail/enumerations.h>
#include <CGAL/Levels_of_detail/internal/structures.h>

namespace CGAL {
namespace Levels_of_detail {
namespace internal {

  template<
  typename GeomTraits, 
  typename InputRange, 
  typename PointMap,
  typename SemanticMap, 
  typename VisibilityMap>
  struct Data_structure {

  public:
        
    // Types.
    using Traits = GeomTraits;
    using Input_range = InputRange;
    using Point_map = PointMap;
    using Semantic_map = SemanticMap;
    using Visibility_map = VisibilityMap;

    using Point_2 = typename Traits::Point_2;
    using Point_3 = typename Traits::Point_3;
    using Plane_3 = typename Traits::Plane_3;

    struct Filter_points_by_label {

    public:
      Semantic_label label;
      SemanticMap semantic_map;

      Filter_points_by_label() { }
          
      Filter_points_by_label(
        Semantic_label label, 
        SemanticMap semantic_map) : 
      label(label), 
      semantic_map(semantic_map) 
      { }

      bool operator()(const typename SemanticMap::key_type& key) const {
        return get(semantic_map, key) == label;
      }
    };

    using Iterator = typename Input_range::const_iterator;
    using Filter_iterator = 
    boost::filter_iterator<Filter_points_by_label, Iterator>;
    using Filtered_range = Iterator_range<Filter_iterator>;
    using Filtered_range_iterator = typename Filtered_range::const_iterator;

    using Polygon_face_2 = Polygon_face_2<Traits>;
    using Building = Building<Traits>;
    using Tree = Tree<Traits>;
    using Smooth_ground = Smooth_ground<Traits>;

    // Input.
    const Input_range& input_range;
    const bool verbose;

    // Property maps.
    Point_map point_map;
    Semantic_map semantic_map;
    Visibility_map visibility_map;

    // Access containers.
    Plane_3 ground_plane;
    std::vector<Point_3> planar_ground;
    std::vector<Point_2> building_boundary_points_2;
    std::vector< std::vector<std::size_t> > building_boundary_indices_2;
    std::vector<Polygon_face_2> building_polygon_faces_2;
    std::vector< std::vector<std::size_t> > building_footprints_2;
    std::vector<Building> buildings;
    std::vector< std::vector<Filtered_range_iterator> > vegetation_clusters;
    std::vector<Tree> trees;
    Smooth_ground smooth_ground;
    std::vector< std::vector<Filtered_range_iterator> > building_clusters;

    // Constructor.
    Data_structure(
      const Input_range& input_range_, 
      Point_map point_map_,
      Semantic_map semantic_map_, 
      Visibility_map visibility_map_,
      const bool verbose_ = false) : 
    input_range(input_range_),
    point_map(point_map_),
    semantic_map(semantic_map_),
    visibility_map(visibility_map_),
    verbose(verbose_) 
    { }

    ~Data_structure() 
    { }

    // Access functions.
    inline Filtered_range ground_points() const {
      return make_range(
        boost::make_filter_iterator(
          Filter_points_by_label(
            Semantic_label::GROUND, semantic_map),
            input_range.begin(), input_range.end()),
        boost::make_filter_iterator(
          Filter_points_by_label(
            Semantic_label::GROUND, semantic_map),
            input_range.end(), input_range.end()));
    }

    inline Filtered_range building_boundary_points() const {
      return make_range(
        boost::make_filter_iterator(
          Filter_points_by_label(
            Semantic_label::BUILDING_BOUNDARY, semantic_map),
            input_range.begin(), input_range.end()),
        boost::make_filter_iterator(
          Filter_points_by_label(
            Semantic_label::BUILDING_BOUNDARY, semantic_map),
            input_range.end(), input_range.end()));
    }

    inline Filtered_range building_interior_points() const {
      return make_range(
        boost::make_filter_iterator(
          Filter_points_by_label(
            Semantic_label::BUILDING_INTERIOR, semantic_map),
            input_range.begin(), input_range.end()),
        boost::make_filter_iterator(
          Filter_points_by_label(
            Semantic_label::BUILDING_INTERIOR, semantic_map),
            input_range.end(), input_range.end()));
    }

    inline Filtered_range vegetation_points() const {
      return make_range(
        boost::make_filter_iterator(
          Filter_points_by_label(
            Semantic_label::VEGETATION, semantic_map),
            input_range.begin(), input_range.end()),
        boost::make_filter_iterator(
          Filter_points_by_label(
            Semantic_label::VEGETATION, semantic_map),
            input_range.end(), input_range.end()));
    }

  }; // Data_structure

} // internal
} // Levels_of_detail
} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_DATA_STRUCTURE_H
