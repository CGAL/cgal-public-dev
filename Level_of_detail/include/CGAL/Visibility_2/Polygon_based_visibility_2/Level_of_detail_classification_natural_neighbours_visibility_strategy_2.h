#ifndef CGAL_LEVEL_OF_DETAIL_CLASSIFICATION_NATURAL_NEIGHBOURS_VISIBILITY_STRATEGY_2_H
#define CGAL_LEVEL_OF_DETAIL_CLASSIFICATION_NATURAL_NEIGHBOURS_VISIBILITY_STRATEGY_2_H

// STL includes.
#include <map>
#include <cmath>
#include <vector>
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
#include <CGAL/Visibility_2/Polygon_based_visibility_2/Level_of_detail_classification_labels_matcher_2.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputContainer, class DataStructure, class VisibilityOutput>
		class Level_of_detail_classification_natural_neighbours_visibility_strategy_2 {

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

            using Polygon                 = typename Container::Polygon;
			using Polygon_vertex_iterator = typename Polygon::Vertex_const_iterator;

            using Points     = Input_container;
            using Visibility = Visibility_output;

            using Labels_matcher = CGAL::LOD::Level_of_detail_classification_labels_matcher_2<Kernel, Visibility>;

            Level_of_detail_classification_natural_neighbours_visibility_strategy_2(const Points &points, const Data_structure &data_structure) : 
            m_points(points), m_data_structure(data_structure), m_scale(-FT(1)) { }

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

            void set_scale(const FT new_value) {
                
                assert(new_value > FT(0));
                m_scale = new_value;
            }

        private:
            const Points         &m_points;
            const Data_structure &m_data_structure;

            FT             m_scale;
            Labels_matcher m_labels_matcher;

            void estimate_polygon_visibility(const Fuzzy_tree &tree, const Polygon &polygon, const size_t container_index, Visibility &visibility) {
                

            }
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_CLASSIFICATION_NATURAL_NEIGHBOURS_VISIBILITY_STRATEGY_2_H