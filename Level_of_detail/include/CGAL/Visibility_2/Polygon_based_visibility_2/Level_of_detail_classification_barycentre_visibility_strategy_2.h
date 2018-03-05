#ifndef CGAL_LEVEL_OF_DETAIL_CLASSIFICATION_BARYCENTRE_VISIBILITY_STRATEGY_2_H
#define CGAL_LEVEL_OF_DETAIL_CLASSIFICATION_BARYCENTRE_VISIBILITY_STRATEGY_2_H

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

            using Polygon                 = typename Container::Polygon;
			using Polygon_vertex_iterator = typename Polygon::Vertex_const_iterator;

            using Points     = Input_container;
            using Visibility = Visibility_output;

            using Labels_matcher = CGAL::LOD::Level_of_detail_classification_labels_matcher_2<Kernel, Visibility>;

            Level_of_detail_classification_barycentre_visibility_strategy_2(const Points &points, const Data_structure &data_structure) : 
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
                
                
                // Compute the number of vertices in the given polygon.
                const size_t num_vertices = std::distance(polygon.vertices_begin(), polygon.vertices_end());
                std::vector<Point_2> vertices(num_vertices);


                // Compute barycentre of the polygon and map its vertices to a temporal data structure.
                FT x = FT(0), y = FT(0); size_t i = 0;
                for (Polygon_vertex_iterator vit = polygon.vertices_begin(); vit != polygon.vertices_end(); ++vit, ++i) {
                    const Point_2 &vertex = *vit;

                    x += vertex.x();
                    y += vertex.y();

                    vertices[i] = Point_2(vertex.x(), vertex.y());
                }

                if (num_vertices == 0) {
                    visibility[container_index] = std::make_pair(FT(0), FT(0));
                    return;
                }

                assert(num_vertices > 0);

                x /= static_cast<FT>(num_vertices);
                y /= static_cast<FT>(num_vertices);

                const Point_2 barycentre(x, y);


                // Compute mean to all the polygon's edges.
                FT min_distance = FT(1000000000000);

                FT mean = FT(0);
                std::vector<FT> distances(num_vertices);

                for (i = 0; i < num_vertices; ++i) {
                    const size_t ip = (i + 1) % num_vertices;

                    const Point_2 &p1 = vertices[i];
                    const Point_2 &p2 = vertices[ip];

                    const Line_2 line(p1, p2);
                    const Point_2 projected = line.projection(barycentre);

                    distances[i] = static_cast<FT>(CGAL::sqrt(CGAL::to_double(CGAL::squared_distance(barycentre, projected))));
                    min_distance = CGAL::min(min_distance, distances[i]);

                    mean += distances[i];
                }
                mean /= static_cast<FT>(distances.size());


                // Compute standard deviation wrt the mean above.
                FT stde = FT(0);
                for (size_t i = 0; i < distances.size(); ++i) stde += (distances[i] - mean) * (distances[i] - mean);
                
                stde /= static_cast<FT>(distances.size());
                stde  = static_cast<FT>(CGAL::sqrt(CGAL::to_double(stde)));


                // Set radius of the disc for searching natural neighbours.
                assert(min_distance > FT(0));
                FT radius = FT(0);
                
                const FT half_mean = mean / FT(2);
                if (stde < half_mean) radius = min_distance * FT(3) / FT(5);
                else radius = min_distance * FT(7) / FT(5);

                const Fuzzy_circle circle(barycentre, radius);


                // Search for natural neighbours.
                Points found_points;
				tree.search(std::back_inserter(found_points), circle);


                // Compute local visibility.
                const size_t local_container_index = 0;

                Visibility local_visibility;
                local_visibility[local_container_index] = std::make_pair(FT(0), FT(0));

                assert(found_points.size() > 0);
                for (i = 0; i < found_points.size(); ++i) {

                    const Point_label found_point_label = found_points[i].second;
                    m_labels_matcher.add_visibility(local_container_index, found_point_label, local_visibility);
                }

                visibility[container_index] = local_visibility.at(local_container_index);
            }
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_CLASSIFICATION_BARYCENTRE_VISIBILITY_STRATEGY_2_H