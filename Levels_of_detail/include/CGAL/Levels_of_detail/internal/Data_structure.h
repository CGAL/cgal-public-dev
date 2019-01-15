#ifndef CGAL_LEVELS_OF_DETAIL_DATA_STRUCTURE_H
#define CGAL_LEVELS_OF_DETAIL_DATA_STRUCTURE_H

// Boost includes.
#include <boost/iterator/filter_iterator.hpp>

// CGAL includes.
#include <CGAL/Iterator_range.h>

namespace CGAL {

  namespace Levels_of_detail {

    namespace internal {

      template<
      typename GeometricTraits, 
      typename InputRange, 
      typename PointMap,
      typename SemanticMap, 
      typename VisibilityMap>
      struct Data_structure {

      public:
        
        using Kernel = GeometricTraits;
        using Input_range = InputRange;
        using Point_map = PointMap;
        using Semantic_map = SemanticMap;
        using Visibility_map = VisibilityMap;

        struct Filter_points_by_label {
          
          Semantic_label label;
          SemanticMap semantic_map;

          Filter_points_by_label() { }
          
          Filter_points_by_label(
            Semantic_label label, 
            SemanticMap semantic_map) : 
          label(label), 
          semantic_map(semantic_map) 
          { }

          bool operator()(const typename SemanticMap::key_type &key) const {
            return get(semantic_map, key) == label;
          }
        };

        using iterator = typename Input_range::const_iterator;
        using Filtered_iterator = 
        boost::filter_iterator<Filter_points_by_label, iterator>;
        using Filtered_range = Iterator_range<Filtered_iterator>;

        const Input_range &input_range;

        Point_map point_map;
        Semantic_map semantic_map;
        Visibility_map visibility_map;

        Data_structure(
          const Input_range &input_range_, 
          Point_map point_map_,
          Semantic_map semantic_map_, 
          Visibility_map visibility_map_) : 
        input_range(input_range_),
        point_map(point_map_),
        semantic_map(semantic_map_),
        visibility_map(visibility_map_)
        { }

        inline Filtered_range ground_points() const {
          return make_range(
            boost::make_filter_iterator(
              Filter_points_by_label(Semantic_label::GROUND, semantic_map),
              input_range.begin(), input_range.end()),
            boost::make_filter_iterator(
              Filter_points_by_label(Semantic_label::GROUND, semantic_map),
              input_range.end(), input_range.end()));
        }

        inline Filtered_range building_boundary_points() const {
          return make_range(
            boost::make_filter_iterator(
              Filter_points_by_label(Semantic_label::BUILDING_BOUNDARY, semantic_map),
              input_range.begin(), input_range.end()),
            boost::make_filter_iterator(
              Filter_points_by_label(Semantic_label::BUILDING_BOUNDARY, semantic_map),
              input_range.end(), input_range.end()));
        }

        inline Filtered_range building_interior_points() const {
          return make_range(
            boost::make_filter_iterator(
              Filter_points_by_label(Semantic_label::BUILDING_INTERIOR, semantic_map),
              input_range.begin(), input_range.end()),
            boost::make_filter_iterator(
              Filter_points_by_label(Semantic_label::BUILDING_INTERIOR, semantic_map),
              input_range.end(), input_range.end()));
        }

        inline Filtered_range vegetation_points() const {
          return make_range(
            boost::make_filter_iterator(
              Filter_points_by_label(Semantic_label::VEGETATION, semantic_map),
              input_range.begin(), input_range.end()),
            boost::make_filter_iterator(
              Filter_points_by_label(Semantic_label::VEGETATION, semantic_map),
              input_range.end(), input_range.end()));
        }

      }; // Data_structure

    } // internal

  } // Levels_of_detail

} // CGAL

#endif // CGAL_LEVELS_OF_DETAIL_DATA_STRUCTURE_H
