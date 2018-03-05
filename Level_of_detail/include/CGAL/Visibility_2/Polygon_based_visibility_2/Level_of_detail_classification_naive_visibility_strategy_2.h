#ifndef CGAL_LEVEL_OF_DETAIL_CLASSIFICATION_NAIVE_VISIBILITY_STRATEGY_2_H
#define CGAL_LEVEL_OF_DETAIL_CLASSIFICATION_NAIVE_VISIBILITY_STRATEGY_2_H

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

// New CGAL includes.
#include <CGAL/Visibility_2/Polygon_based_visibility_2/Level_of_detail_classification_labels_matcher_2.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputContainer, class DataStructure, class VisibilityOutput>
		class Level_of_detail_classification_naive_visibility_strategy_2 {

        public:
            typedef KernelTraits     Kernel;
            typedef InputContainer   Input_container;
            typedef DataStructure    Data_structure;
            typedef VisibilityOutput Visibility_output;

            using FT      = typename Kernel::FT;
            using Point_2 = typename Kernel::Point_2;

            using Point_label = int;
            
            using Points     = Input_container;
            using Visibility = Visibility_output;

            using Labels_matcher = CGAL::LOD::Level_of_detail_classification_labels_matcher_2<Kernel, Visibility>;

            Level_of_detail_classification_naive_visibility_strategy_2(const Points &points, const Data_structure &data_structure) : 
            m_points(points), m_data_structure(data_structure), m_scale(-FT(1)) { }

            void estimate(Visibility &visibility) {
                
                assert(m_points.size()   > 0);
                assert(visibility.size() > 0);
                
				for (size_t i = 0; i < m_points.size(); ++i) {
                    const Point_2 &point = m_points[i].first;

                    const int container_index = m_data_structure.locate(point);
                    if (container_index < 0) continue;

                    const Point_label point_label = m_points[i].second;
                    m_labels_matcher.add_visibility(static_cast<size_t>(container_index), point_label, visibility);
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
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_CLASSIFICATION_NAIVE_VISIBILITY_STRATEGY_2_H