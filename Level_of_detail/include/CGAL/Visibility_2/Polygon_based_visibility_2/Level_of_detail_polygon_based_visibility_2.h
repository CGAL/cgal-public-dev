#ifndef CGAL_LEVEL_OF_DETAIL_POLYGON_BASED_VISIBILITY_2_H
#define CGAL_LEVEL_OF_DETAIL_POLYGON_BASED_VISIBILITY_2_H

#if defined(WIN32) || defined(_WIN32) 
#define PSR "\\" 
#else 
#define PSR "/" 
#endif

// STL includes.
#include <map>
#include <cmath>
#include <memory>
#include <vector>
#include <utility>
#include <cassert>
#include <iostream>

// Boost includes.
#include <boost/tuple/tuple.hpp>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/IO/Color.h>
#include <CGAL/number_utils.h>
#include <CGAL/property_map.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputContainer, class DataStructure, class VisibilityStrategy>
		class Level_of_detail_polygon_based_visibility_2 {

        public:
            typedef KernelTraits       Kernel;
            typedef InputContainer     Input;
            typedef DataStructure      Data_structure;
            typedef VisibilityStrategy Visibility_strategy;

            using FT  = typename Kernel::FT;
            using Log = CGAL::LOD::Mylog;

            using Point_2 = typename Kernel::Point_2;
            using Point_3 = typename Kernel::Point_3;

            using Point_label     = int;
            using Point_label_map = typename Input:: template Property_map<Point_label>;

			using Point_with_label = typename std::pair<Point_2, Point_label>;
			using Point_map        = typename CGAL::First_of_pair_property_map<Point_with_label>;
            using Points           = std::vector<Point_with_label>;

            using Container  = typename Data_structure::Container;
            using Containers = typename Data_structure::Containers;
				
			using Polygon                 = typename Container::Polygon;
			using Polygon_vertex_iterator = typename Polygon::Vertex_const_iterator;

            using Visibility_data = std::pair<FT, FT>;
            using Visibility      = std::map<size_t, Visibility_data>;

            using Point_iterator = typename Input::const_iterator;
            using Colour         = typename Container::Colour;

            Level_of_detail_polygon_based_visibility_2(const Input &input, Data_structure &data_structure) :
            m_silent(false), m_input(input), m_data_structure(data_structure) {
                
                set_point_labels();
            }

            void compute() {
                
                set_points();
                Visibility visibility;

                set_initial_visibility(visibility);
                estimate_visibility(visibility);
                
                set_visibility_to_containers(visibility);
                if (!m_silent) save_data_structure();
            }

            void make_silent(const bool new_state) {
                m_silent = new_state;
            }

        private:
            bool            m_silent;
            Points          m_points;
            Point_label_map m_point_labels;
            
            const Input    &m_input;
            Data_structure &m_data_structure;

            std::shared_ptr<Visibility_strategy> m_visibility_strategy;
            
			inline void set_point_labels() {
				boost::tie(m_point_labels, boost::tuples::ignore) = m_input.template property_map<Point_label>("label");
			}

            void set_points() {

                m_points.clear();
                m_points.resize(m_input.number_of_points());

                Point_2 point; size_t i = 0;
                for (Point_iterator pit = m_input.begin(); pit != m_input.end(); ++pit, ++i) {

                    const Point_3 &original = m_input.point(*pit);
                    point = Point_2(original.x(), original.y());

                    m_points[i] = std::make_pair(point, get_point_label(pit));
                }
                assert(i == m_input.number_of_points());
            }

            Point_label get_point_label(const Point_iterator &pit) {
                return m_point_labels[*pit];
            }

            void set_initial_visibility(Visibility &visibility) {
                
                const Containers &containers = m_data_structure.containers();
				assert(containers.size() > 0);

                visibility.clear();
                for (size_t i = 0; i < containers.size(); ++i)
                    visibility[i] = std::make_pair(FT(0), FT(0));
            }

            inline void estimate_visibility(Visibility &visibility) {
                
                assert(m_points.size() > 0);
                
                m_visibility_strategy = std::make_shared<Visibility_strategy>(m_points, m_data_structure);
                m_visibility_strategy->estimate(visibility);
            }

            void set_visibility_to_containers(const Visibility &visibility) {

                Containers &containers = m_data_structure.containers();
				assert(containers.size() > 0);

                assert(visibility.size() == containers.size());
                for (size_t i = 0; i < containers.size(); ++i) {

                    containers[i].inside = get_final_visibility_value(visibility.at(i));
                    containers[i].colour = get_final_visibility_colour(containers[i].inside);
                }
            }

            FT get_final_visibility_value(const Visibility_data &visibility_data) {

                const FT inside  = visibility_data.first;
                const FT outside = visibility_data.second;

                if (inside > outside) return FT(1);
                return FT(0);
            }

            Colour get_final_visibility_colour(const FT visibility_value) {
                assert(visibility_value == FT(0) || visibility_value == FT(1));

				if (visibility_value == FT(1)) return Colour(55, 255, 55);
                return Colour(255, 55, 55);
            }

            void save_data_structure() {
                
                const Containers &containers = m_data_structure.containers();
                assert(containers.size() > 0);

                Log exporter;
                exporter.save_polygons<Containers, Polygon, Kernel>(containers, "tmp" + std::string(PSR) + "visibility", true);
            }
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_POLYGON_BASED_VISIBILITY_2_H