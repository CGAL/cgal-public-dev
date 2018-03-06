#ifndef CGAL_LEVEL_OF_DETAIL_CLASSIFICATION_SHEPARD_VISIBILITY_STRATEGY_2_H
#define CGAL_LEVEL_OF_DETAIL_CLASSIFICATION_SHEPARD_VISIBILITY_STRATEGY_2_H

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

#include <CGAL/Kd_tree.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_adapter.h>

// New CGAL includes.
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Visibility_2/Polygon_based_visibility_2/Level_of_detail_polygon_data_estimator_2.h>
#include <CGAL/Visibility_2/Polygon_based_visibility_2/Level_of_detail_classification_labels_matcher_2.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputContainer, class DataStructure, class VisibilityOutput>
		class Level_of_detail_classification_shepard_visibility_strategy_2 {

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
            using Polygon_data_estimator = CGAL::LOD::Level_of_detail_polygon_data_estimator_2<Kernel, Polygon>;

            using Log = CGAL::LOD::Mylog;

            Level_of_detail_classification_shepard_visibility_strategy_2(const Points &points, const Data_structure &data_structure) : 
            m_points(points), m_data_structure(data_structure), m_debug(false){ }

            void estimate(Visibility &visibility) {
                
                assert(visibility.size() > 0);
                if (m_debug) m_debug_samples.clear();

                const Containers &containers = m_data_structure.containers();
				assert(containers.size() > 0 && containers.size() == visibility.size());

                assert(m_points.size() > 0);
                Fuzzy_tree tree(m_points.begin(), m_points.end());

                for (size_t i = 0; i < containers.size(); ++i) {
                    
                    const Polygon &polygon = containers[i].polygon;
                    estimate_polygon_visibility(tree, polygon, i, visibility);
                }

                if (m_debug) {
                    Log log; log.export_points("tmp" + std::string(PSR) + "samples_shepard_visibility", m_debug_samples);
                }
            }

        private:
            const Points         &m_points;
            const Data_structure &m_data_structure;

            Labels_matcher                          m_labels_matcher;
            std::shared_ptr<Polygon_data_estimator> m_polygon_data_estimator;

            bool m_debug, m_extra_sample;
            std::vector<Point_2> m_debug_samples;

            void estimate_polygon_visibility(const Fuzzy_tree &tree, const Polygon &polygon, const size_t container_index, Visibility &visibility) {
                
                visibility[container_index] = std::make_pair(FT(0), FT(0));

                // Some assertions.
                if (polygon.size() == 0) return;

                // Compute some preliminary data.
                Point_2 polygon_barycentre;
                m_polygon_data_estimator = std::make_shared<Polygon_data_estimator>(polygon);

                m_polygon_data_estimator->compute_barycentre(polygon_barycentre);
                m_polygon_data_estimator->compute_distances_to_edges(polygon_barycentre);

                const FT min_distance = m_polygon_data_estimator->compute_minimum_distance_to_edges();

                const FT mean = m_polygon_data_estimator->compute_mean_distance_to_edges();
                const FT stde = m_polygon_data_estimator->compute_standard_deviation_from_distances_to_edges(mean);

                // Set radius of the disc for searching natural neighbours.
                assert(min_distance > FT(0));
                FT radius = FT(0);
                
                const FT half_mean = mean / FT(2);
                if (stde < half_mean) radius = min_distance * FT(3) / FT(5);
                else radius = min_distance * FT(7) / FT(5);

                // Sample polygon.
                std::vector<Point_2> samples;
                samples.push_back(polygon_barycentre);

                if (m_debug) {
                    for (size_t i = 0; i < samples.size(); ++i) 
                        m_debug_samples.push_back(samples[i]);
                }

                // Compute visibility.
                Points found_points;
                for (size_t i = 0; i < samples.size(); ++i) {
				    
                    const Point_2 &query = samples[i];
                    const Fuzzy_circle circle(query, radius);

                    found_points.clear();
                    tree.search(std::back_inserter(found_points), circle);

                    add_visibility(query, found_points, container_index, visibility);
                }
            }

            void add_visibility(const Point_2 &query, const Points &points, const size_t index, Visibility &result) {

                assert(points.size() > 0);
                std::vector<FT> values(points.size());
                
                for (size_t i = 0; i < points.size(); ++i) {
                    
                    const Point_label point_label = points[i].second;
                    values[i] = m_labels_matcher.match_label(point_label);
                }

                const FT intp_value = interpolate(query, points, values);

                FT inside = FT(0), outside = FT(0);
                if (intp_value > FT(1) / FT(2)) inside += FT(1);
                else outside += FT(1);

                result[index].first  += inside;
                result[index].second += outside;
            }

            FT interpolate(const Point_2 &query, const Points &points, const std::vector<FT> &values) {

                std::vector<FT> weights;
                compute_shepard_weights(query, points, weights);

                FT result = FT(0);
				for (size_t i = 0; i < values.size(); ++i) {
					
					assert(values[i] >= FT(0) && values[i] <= FT(1));
					result += values[i] * weights[i];
				}

                if (result < FT(0)) result = FT(0);
                if (result > FT(1)) result = FT(1);

				return result;
            }

            void compute_shepard_weights(const Point_2 &query, const Points &points, std::vector<FT> &weights) {

                weights.clear();
                weights.resize(points.size(), FT(0));

                FT sum = FT(0);
                for (size_t i = 0; i < points.size(); ++i) {
                    const FT distance = static_cast<FT>(CGAL::sqrt(CGAL::to_double(CGAL::squared_distance(query, points[i].first))));

                    if (distance == FT(0)) {
                     
                        weights.clear();
                        weights.resize(points.size(), FT(0));
                     
                        weights[i] = FT(1);
                        return;
                    }

                    weights[i] = FT(1) / distance;
                    sum += weights[i];
                }

                for (size_t i = 0; i < weights.size(); ++i) weights[i] /= sum;
            }
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_CLASSIFICATION_SHEPARD_VISIBILITY_STRATEGY_2_H