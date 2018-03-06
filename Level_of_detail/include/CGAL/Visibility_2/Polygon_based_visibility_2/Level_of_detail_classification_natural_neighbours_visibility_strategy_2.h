#ifndef CGAL_LEVEL_OF_DETAIL_CLASSIFICATION_NATURAL_NEIGHBOURS_VISIBILITY_STRATEGY_2_H
#define CGAL_LEVEL_OF_DETAIL_CLASSIFICATION_NATURAL_NEIGHBOURS_VISIBILITY_STRATEGY_2_H

#if defined(WIN32) || defined(_WIN32) 
#define PSR "\\" 
#else 
#define PSR "/" 
#endif 

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

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Interpolation_traits_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Visibility_2/Polygon_based_visibility_2/Level_of_detail_polygon_sampler_2.h>
#include <CGAL/Visibility_2/Polygon_based_visibility_2/Level_of_detail_polygon_data_estimator_2.h>
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

            using Labels_matcher         = CGAL::LOD::Level_of_detail_classification_labels_matcher_2<Kernel, Visibility>;
            using Polygon_sampler        = CGAL::LOD::Level_of_detail_polygon_sampler_2<Kernel, Polygon>;
            using Polygon_data_estimator = CGAL::LOD::Level_of_detail_polygon_data_estimator_2<Kernel, Polygon>;

            using Delaunay_triangulation = CGAL::Delaunay_triangulation_2<Kernel>;
			using Interpolation_traits   = CGAL::Interpolation_traits_2<Kernel>;
			
			using Function_type = std::map<Point_2, FT, typename Kernel::Less_xy_2>;
			using Value_access  = CGAL::Data_access<Function_type>;

            using Log = CGAL::LOD::Mylog;

            Level_of_detail_classification_natural_neighbours_visibility_strategy_2(const Points &points, const Data_structure &data_structure) : 
            m_points(points), m_data_structure(data_structure), m_norm_threshold(FT(1000)), m_debug(false) { }

            void estimate(Visibility &visibility) {
                
                assert(visibility.size() > 0);
                if (m_debug) m_debug_samples.clear();

                const Containers &containers = m_data_structure.containers();
				assert(containers.size() > 0 && containers.size() == visibility.size());

				assert(m_points.size() > 0);
                set_delunay_and_function_values();
                
                for (size_t i = 0; i < containers.size(); ++i) {
                    
                    const Polygon &polygon = containers[i].polygon;
                    estimate_polygon_visibility(polygon, i, visibility);
                }

                if (m_debug) {
                    Log log; log.export_points("tmp" + std::string(PSR) + "samples_shepard_visibility", m_debug_samples);
                }
            }

        private:
            const Points         &m_points;
            const Data_structure &m_data_structure;

            Labels_matcher m_labels_matcher;
            FT             m_norm_threshold;

            std::shared_ptr<Polygon_data_estimator> m_polygon_data_estimator;
            std::shared_ptr<Polygon_sampler>        m_polygon_sampler;

            bool m_debug;
            std::vector<Point_2> m_debug_samples;

            Delaunay_triangulation m_dt;
			Function_type          m_function_values;

			void set_delunay_and_function_values() {

                assert(m_points.size() > 0);
                for (size_t i = 0; i < m_points.size(); ++i) {
                    
                    const Point_label point_label = m_points[i].second;
                    const FT inside = m_labels_matcher.match_label(point_label);

                    m_dt.insert(m_points[i].first);
					m_function_values.insert(std::make_pair(m_points[i].first, inside));
                }
			}

            void estimate_polygon_visibility(const Polygon &polygon, const size_t container_index, Visibility &visibility) {
                
                visibility[container_index] = std::make_pair(FT(0), FT(0));

                // Some assertions.
                if (polygon.size() == 0) return;

                // Compute some preliminary data.
                Point_2 polygon_barycentre;

                m_polygon_data_estimator = std::make_shared<Polygon_data_estimator>(polygon);
                m_polygon_data_estimator->compute_barycentre(polygon_barycentre);

                // Sample polygon.
                std::vector<Point_2> samples;
                m_polygon_sampler = std::make_shared<Polygon_sampler>(polygon);

                m_polygon_sampler->set_number_of_subdivision_steps(2);
                m_polygon_sampler->create_samples(samples);

                samples.push_back(polygon_barycentre);

                if (m_debug) {
                    for (size_t i = 0; i < samples.size(); ++i) 
                        m_debug_samples.push_back(samples[i]);
                }

                // Compute visibility.
                for (size_t i = 0; i < samples.size(); ++i) {
				    
                    const Point_2 &query = samples[i];
                    add_visibility(query, container_index, visibility);
                }
            }

            void add_visibility(const Point_2 &query, const size_t index, Visibility &result) {

                FT inside = FT(0), outside = FT(0);
                const FT intp_value = interpolate(query);

                if (intp_value > FT(1) / FT(2)) inside += FT(1);
                else outside += FT(1);

                result[index].first  += inside;
                result[index].second += outside;
            }

            FT interpolate(const Point_2 &query) {
                
               	std::vector<std::pair<Point_2, FT> > coords;
				const auto triple = CGAL::natural_neighbor_coordinates_2(m_dt, query, std::back_inserter(coords));

				const bool success = triple.third;
				const FT norm      = static_cast<FT>(triple.second);

				if (!success) 			   return FT(0);
				if (is_invalid_norm(norm)) return FT(0);

				assert(norm > FT(0));
				const FT intp = CGAL::linear_interpolation(coords.begin(), coords.end(), norm, Value_access(m_function_values));
				return intp;
            }

            bool is_invalid_norm(const FT norm) {
				return (!std::isfinite(CGAL::to_double(norm)) || norm <= FT(0) || norm > m_norm_threshold);
			}
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_CLASSIFICATION_NATURAL_NEIGHBOURS_VISIBILITY_STRATEGY_2_H