#ifndef CGAL_LEVEL_OF_DETAIL_POLYGON_BASED_VISIBILITY_2_H
#define CGAL_LEVEL_OF_DETAIL_POLYGON_BASED_VISIBILITY_2_H

// STL includes.
#include <utility>
#include <cassert>
#include <vector>
#include <map>
#include <cmath>
#include <iostream>

// Boost includes.
#include <boost/tuple/tuple.hpp>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/IO/Color.h>
#include <CGAL/number_utils.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputContainer, class DataStructure>
		class Level_of_detail_polygon_based_visibility_2 {

        public:
            typedef KernelTraits   Kernel;
            typedef InputContainer Input;
            typedef DataStructure  Data_structure;

            using FT = typename Kernel::FT;
            using Log = CGAL::LOD::Mylog;
            
            using Point_label = int;
            using Point_label_map = typename Input:: template Property_map<Point_label>;

            using Point_2 = typename Kernel::Point_2;
            using Point_3 = typename Kernel::Point_3;

            using Containers = typename Data_structure::Containers;
			using Container  = typename Data_structure::Container;
				
			using Polygon                 = typename Container::Polygon;
			using Polygon_vertex_iterator = typename Polygon::Vertex_const_iterator;

            using Visibility_data = std::pair<FT, FT>;
            using Visibility      = std::map<size_t, Visibility_data>;

            using Point_iterator = typename Input::const_iterator;
            using Colour = typename Container::Colour;

            Level_of_detail_polygon_based_visibility_2(const Input &input, Data_structure &data_structure) :
            m_silent(false), m_input(input), m_data_structure(data_structure) { 
                
                set_point_labels();
            }

            void compute() {

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
            bool m_silent;
            Point_label_map m_point_labels;

            const Input    &m_input;
            Data_structure &m_data_structure;

			inline void set_point_labels() {
				boost::tie(m_point_labels, boost::tuples::ignore) = m_input.template property_map<Point_label>("label");
			}

            void set_initial_visibility(Visibility &visibility) {
                
                const Containers &containers = m_data_structure.containers();
				assert(containers.size() > 0);

                visibility.clear();
                for (size_t i = 0; i < containers.size(); ++i)
                    visibility[i] = std::make_pair(FT(0), FT(0));
            }

            void estimate_visibility(Visibility &visibility) {

                Point_2 point;
				for (Point_iterator pit = m_input.begin(); pit != m_input.end(); ++pit) {
					
                    const Point_3 &original = m_input.point(*pit);
                    point = Point_2(original.x(), original.y());

                    // Change it to more precise location type, for edge - add value to both faces, for vertex to all adjacent faces.
                    const int container_index = m_data_structure.locate(point);
                    if (container_index < 0) continue; 

                    const Point_label point_label = get_point_label(pit);
                    add_visibility(static_cast<size_t>(container_index), point_label, visibility);
                }
            }

            Point_label get_point_label(const Point_iterator &pit) {
                assert(m_point_labels.size() == m_input.number_of_points());
                return m_point_labels[*pit];
            }

            void add_visibility(const size_t container_index, const Point_label point_label, Visibility &visibility) {

				const Point_label ground     = 0;
				const Point_label facade     = 1;
				const Point_label roof       = 2;
				const Point_label vegetation = 3;

				switch (point_label) {

					case ground:
						set_outside(container_index, visibility);
						break;

					case facade:
						set_unknown(container_index, visibility);
						break;

					case roof:
						set_inside(container_index, visibility);
						break;

					case vegetation:
						set_outside(container_index, visibility);
						break;

					default:
                        set_unknown(container_index, visibility);
                        break;
				}
			}

            inline void set_inside(const size_t container_index, Visibility &visibility) {
				visibility[container_index].first += FT(1);
			}

			inline void set_outside(const size_t container_index, Visibility &visibility) {
				visibility[container_index].second += FT(1);
			}

			void set_unknown(const size_t container_index, Visibility &visibility) {
				visibility[container_index].first  += FT(1) / FT(2);
                visibility[container_index].second += FT(1) / FT(2);
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

                if (inside >= outside) return FT(1);
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
                exporter.save_polygons<Containers, Polygon, Kernel>(containers, "tmp" + std::string(PS) + "visibility", true);
            }
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_POLYGON_BASED_VISIBILITY_2_H