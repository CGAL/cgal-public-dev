#ifndef CGAL_LEVEL_OF_DETAIL_FACET_VISIBILITY_ESTIMATOR_H
#define CGAL_LEVEL_OF_DETAIL_FACET_VISIBILITY_ESTIMATOR_H

// STL includes.
#include <vector>
#include <memory>

// LOD includes.
#include <CGAL/Level_of_detail/Enumerations.h>
#include <CGAL/Level_of_detail/internal/Data/Polygon_data_estimator.h>
#include <CGAL/Level_of_detail/internal/Data/Kd_tree_with_data_creator.h>
#include <CGAL/Level_of_detail/internal/Property_maps/Point_from_value_property_map_2.h>
#include <CGAL/Level_of_detail/internal/Partitioning/Partition_element.h>

namespace CGAL {

	namespace Level_of_detail {

        namespace LOD = CGAL::Level_of_detail;

		template<class InputKernel, class InputRange, class PointMap, class VisibilityMap>
		class Facet_visibility_estimator {
			
        public:
            using Kernel      = InputKernel;
            using Input_range = InputRange;
            using Point_map_3 = PointMap;
            using Visibility_map   = VisibilityMap;

            using FT      = typename Kernel::FT;
            using Point_2 = typename Kernel::Point_2;
            using Vector_2 = typename Kernel::Vector_2;
            using Triangle_2 = typename Kernel::Triangle_2;

            typename Kernel::Compute_squared_distance_2 squared_distance_2;

            using Elements           = Input_range;
            using Element_identifier = typename Elements::Index;
            using Point_map_2        = LOD::Point_from_value_property_map_2<Element_identifier, Point_2, Point_map_3>;
            using Points_tree_2      = LOD::Kd_tree_with_data_creator<Kernel, Element_identifier, Elements, Point_map_2>;

            using Neighbours                = typename Points_tree_2::Neighbours;
            using Const_neighbours_iterator = typename Neighbours::const_iterator;

            using Function_values = std::vector<FT>;

            using Self     = Facet_visibility_estimator<Kernel, Input_range,
                                                        Point_map_3, Visibility_map>;

            using Facet = Partition_element<Kernel, CGAL::Polygon_2<Kernel> >;

            using Facet_data_estimator = LOD::Polygon_data_estimator<Kernel, Facet>;
            using Visibility_label     = LOD::Visibility_label;

            Facet_visibility_estimator(
                const Input_range &input_range, 
                const Point_map_3 &point_map_3, 
                const Visibility_map &visibility_map) :
            m_input_range(input_range),
            m_point_map_3(point_map_3),
            m_visibility_map(visibility_map) { 

                const FT local_search_radius = FT(0);

                m_point_map_2 = std::make_shared<Point_map_2>(m_point_map_3);
                m_tree        = std::make_shared<Points_tree_2>(m_input_range, *m_point_map_2.get(), local_search_radius);
            }

            void estimate_visibility(Facet& facet) {
              
//                compute_shepard_based_label(facet);
              compute_monte_carlo_label(facet, 100);
            }

        private:
            const Input_range &m_input_range;
            const Point_map_3 &m_point_map_3;
            const Visibility_map   &m_visibility_map;

            std::shared_ptr<Point_map_2>   m_point_map_2;
            std::shared_ptr<Points_tree_2> m_tree;

            void compute_monte_carlo_label(Facet& facet, std::size_t nb_probes) const {

              std::vector<std::pair<Triangle_2, double> > proba;
              double area = 0.;
              for (std::size_t i = 1; i < facet.size() - 1; ++ i)
              {
                Triangle_2 tri(facet[0], facet[i], facet[i+1]);
                proba.push_back (std::make_pair(tri, area));
                area += CGAL::abs (tri.area());
              }
              
              proba.push_back (std::make_pair(Triangle_2(), area));
              
              FT mean_value = FT(0.);
              for (std::size_t i = 0; i < nb_probes; ++ i)
              {
                Neighbours facet_barycentre_neighbours;
                m_tree->search_knn_2(random_point_in_triangles(proba), facet_barycentre_neighbours);
                mean_value +=
                  get_function_value(get(m_tree->element_identifier_map(), *(facet_barycentre_neighbours.begin())));
              }
              mean_value /= nb_probes;
              
              if (mean_value >= FT(1) / FT(2))
                facet.visibility_label() = Visibility_label::INSIDE;
              else
                facet.visibility_label() = Visibility_label::OUTSIDE;
            }


            Point_2 random_point_in_triangles (const std::vector<std::pair<Triangle_2, double> >& proba) const {

              double key = proba.back().second * (rand() / double(RAND_MAX));
              for (std::size_t j = 0; j < proba.size () - 1; ++ j)
                if (proba[j].second < key && key <= proba[j+1].second)
                  return random_point_in_triangle (proba[j].first);
              std::cerr << "Error! probability out of range" << std::endl;
              return Point_2();
            }
      
            Point_2 random_point_in_triangle (const Triangle_2& t) const {
              Vector_2 v01 (t[0], t[1]);
              Vector_2 v02 (t[0], t[2]);
              Point_2 out = t[0];

              double r01 = rand () / (double)RAND_MAX;
              double r02 = rand () / (double)RAND_MAX;

              if (r01 + r02 > 1.0)
              {
                r01 = 1. - r01;
                r02 = 1. - r02;
              }

              out = out + r01 * v01;
              out = out + r02 * v02;
              return out;

            }
              
            void compute_shepard_based_label(Facet &facet) const {

                Point_2 facet_barycentre;
                const FT local_search_radius = compute_local_search_radius(facet, facet_barycentre);

                Neighbours facet_barycentre_neighbours;
                m_tree->set_local_search_radius(local_search_radius);
                m_tree->search_2(facet_barycentre, facet_barycentre_neighbours);

                estimate_visibility(facet_barycentre, facet_barycentre_neighbours, facet);
            }

            FT compute_local_search_radius(const Facet &facet, Point_2 &facet_barycentre) const {
                
                Facet_data_estimator facet_data_estimator(facet);

                facet_barycentre = facet_data_estimator.compute_barycentre();
                facet_data_estimator.compute_distances_to_boundaries(facet_barycentre);

                const FT min_distance = facet_data_estimator.compute_minimum_distance_to_boundaries();

                const FT mean = facet_data_estimator.compute_mean_distance_to_boundaries();
                const FT stde = facet_data_estimator.compute_standard_deviation_on_distances_to_boundaries(mean);

                CGAL_precondition(min_distance > FT(0));
                FT local_search_radius = FT(0);
                
                const FT half_mean = mean / FT(2);
                if (stde < half_mean) local_search_radius = min_distance * FT(3) / FT(5);
                else local_search_radius = min_distance * FT(7) / FT(5);

                return local_search_radius;
            }

            void estimate_visibility(const Point_2 &query, const Neighbours &query_neighbours, Facet &facet) const {
              static std::ofstream f("no_point.xyz");
              f.precision(18);
                if (query_neighbours.size() == 0) {

                  f << query << " " << 0. << std::endl;
                    facet.visibility_label() = Visibility_label::OUTSIDE; return;
                }

                Function_values function_values;
                compute_function_values(query_neighbours, function_values);

                const FT interpolated_value = interpolate(query, query_neighbours, function_values);

                if (interpolated_value >= FT(1) / FT(2)) facet.visibility_label() = Visibility_label::INSIDE;
                else facet.visibility_label() = Visibility_label::OUTSIDE;
            }

            void compute_function_values(const Neighbours &points, Function_values &function_values) const {
                
                function_values.clear();
                function_values.resize(points.size(), FT(0));
                
                size_t i = 0;
                for (Const_neighbours_iterator cn_it = points.begin(); cn_it != points.end(); ++cn_it, ++i)
                    function_values[i] = get_function_value(get(m_tree->element_identifier_map(), *cn_it));
            }

            FT get_function_value(const Element_identifier &element_id) const {
                return get(m_visibility_map, element_id);
            }

            FT interpolate(const Point_2 &query, const Neighbours &query_neighbours, const Function_values &function_values) const {
                
                Function_values shepard_weights;
                compute_shepard_weights(query, query_neighbours, shepard_weights);

                FT result = FT(0);
				for (size_t i = 0; i < function_values.size(); ++i) {
					
					CGAL_precondition(function_values[i] >= FT(0) && function_values[i] <= FT(1));
					result += function_values[i] * shepard_weights[i];
				}

                if (result < FT(0)) result = FT(0);
                if (result > FT(1)) result = FT(1);

				return result;
            }

            void compute_shepard_weights(const Point_2 &query, const Neighbours &query_neighbours, Function_values &weights) const {

                weights.clear();
                weights.resize(query_neighbours.size(), FT(0));

                FT sum_weights = FT(0); size_t i = 0;
                for (Const_neighbours_iterator cn_it = query_neighbours.begin(); cn_it != query_neighbours.end(); ++cn_it, ++i) {
                    
                    const Point_2 &point = get(m_tree->point_map(), *cn_it);
                    const FT distance    = static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_distance_2(query, point))));

                    if (distance == FT(0)) {
                        weights.clear();

                        weights.resize(query_neighbours.size(), FT(0));
                        weights[i] = FT(1); return;
                    }

                    weights[i] = FT(1) / distance;
                    sum_weights += weights[i];
                }
                for (size_t i = 0; i < weights.size(); ++i) weights[i] /= sum_weights;
            }
		};

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_FACET_VISIBILITY_ESTIMATOR_H
