#ifndef CGAL_LEVEL_OF_DETAIL_CLASSIFICATION_BARYCENTRE_VISIBILITY_STRATEGY_2_H
#define CGAL_LEVEL_OF_DETAIL_CLASSIFICATION_BARYCENTRE_VISIBILITY_STRATEGY_2_H

// STL includes.
#include <map>
#include <cmath>
#include <vector>
#include <memory>
#include <utility>
#include <cassert>

// Boost includes.
#include <boost/tuple/tuple.hpp>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/property_map.h>

#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_adapter.h>

// New CGAL includes.
#include <CGAL/Visibility_2/Polygon_based_visibility_2/Level_of_detail_polygon_data_estimator_2.h>
#include <CGAL/Visibility_2/Polygon_based_visibility_2/Level_of_detail_classification_labels_matcher_2.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputContainer, class DataStructure, class VisibilityOutput>
		class Level_of_detail_classification_barycentre_visibility_strategy_2 {

        public:
            typedef KernelTraits     Kernel;
            typedef InputContainer   Input_container;
            typedef DataStructure    Data_structure;
            typedef VisibilityOutput Visibility_output;

            typename Kernel::Compute_squared_distance_2 squared_distance;

            using FT      = typename Kernel::FT;
            using Point_2 = typename Kernel::Point_2;
            using Line_2  = typename Kernel::Line_2;

            using Point_label      = int;
            using Point_with_label = typename std::pair<Point_2, Point_label>;
			using Point_map        = typename CGAL::First_of_pair_property_map<Point_with_label>;
            
			using Search_traits_2 = CGAL::Search_traits_2<Kernel>;
			using Search_traits   = CGAL::Search_traits_adapter<Point_with_label, Point_map, Search_traits_2>;
			using Fuzzy_circle    = CGAL::Fuzzy_sphere<Search_traits>;
			using Fuzzy_tree      = CGAL::Kd_tree<Search_traits>;

            using Container  = typename Data_structure::Container;
            using Containers = typename Data_structure::Containers;

            using Polygon = typename Container::Polygon;

            using Points     = Input_container;
            using Visibility = Visibility_output;

            using Labels_matcher         = CGAL::LOD::Level_of_detail_classification_labels_matcher_2<Kernel, Visibility>;
            using Polygon_data_estimator = CGAL::LOD::Level_of_detail_polygon_data_estimator_2<Kernel, Polygon>;

            Level_of_detail_classification_barycentre_visibility_strategy_2(const Points &points, const Data_structure &data_structure) : 
            m_points(points), m_data_structure(data_structure) { }

            void estimate(Visibility &visibility) {
                assert(visibility.size() > 0);

                const Containers &containers = m_data_structure.containers();
				assert(containers.size() > 0 && containers.size() == visibility.size());

                assert(m_points.size() > 0);
                Fuzzy_tree tree(m_points.begin(), m_points.end());

                for (size_t i = 0; i < containers.size(); ++i) {
                    
                    const Polygon &polygon = containers[i].polygon;
                    estimate_polygon_visibility(tree, polygon, i, visibility);
                }
            }

        private:
            const Points         &m_points;
            const Data_structure &m_data_structure;

            Labels_matcher m_labels_matcher;
            std::shared_ptr<Polygon_data_estimator> m_polygon_data_estimator;

            void estimate_polygon_visibility(const Fuzzy_tree &tree, const Polygon &polygon, const size_t container_index, Visibility &visibility) {

                // Some assertions.
                if (polygon.size() == 0) {
                    
                    visibility[container_index] = std::make_pair(FT(0), FT(0));
                    return;
                }

                // Compute some preliminary data.
                Point_2 polygon_barycentre;
                m_polygon_data_estimator = std::make_shared<Polygon_data_estimator>(polygon);

                m_polygon_data_estimator->compute_barycentre(polygon_barycentre);
                m_polygon_data_estimator->compute_distances_to_edges(polygon_barycentre);

                const FT min_distance = m_polygon_data_estimator->compute_minimum_distance_to_edges();

                const FT mean = m_polygon_data_estimator->compute_mean_distance_to_edges();
                const FT stde = m_polygon_data_estimator->compute_standard_deviation_from_distances_to_edges(mean);

                // Set disc for searching natural neighbours.
                assert(min_distance > FT(0));
                FT radius = FT(0);
                
                const FT half_mean = mean / FT(2);
                if (stde < half_mean) radius = min_distance * FT(3) / FT(5);
                else radius = min_distance * FT(7) / FT(5);

                const Fuzzy_circle circle(polygon_barycentre, radius);

                // Search for natural neighbours.
                Points found_points;
				tree.search(std::back_inserter(found_points), circle);

                // Compute local visibility.
                compute_visibility(found_points, container_index, visibility);
            }

            void compute_visibility(const Points &points, const size_t index, Visibility &result) const {
                const size_t container_index = 0;

                Visibility visibility;
                visibility[container_index] = std::make_pair(FT(0), FT(0));

                assert(points.size() > 0);
                for (size_t i = 0; i < points.size(); ++i) {

                    const Point_label point_label = points[i].second;
                    m_labels_matcher.add_visibility(container_index, point_label, visibility);
                }

                result[index] = visibility.at(container_index);
            }
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_CLASSIFICATION_BARYCENTRE_VISIBILITY_STRATEGY_2_H