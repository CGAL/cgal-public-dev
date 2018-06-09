#ifndef CGAL_LEVEL_OF_DETAIL_POLYGON_DATA_ESTIMATOR_H
#define CGAL_LEVEL_OF_DETAIL_POLYGON_DATA_ESTIMATOR_H

// STL includes.
#include <map>
#include <cmath>
#include <vector>
#include <utility>
#include <cassert>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/property_map.h>

namespace CGAL {

	namespace Level_of_detail {

		template<class KernelTraits, class InputPolygon>
		class Polygon_data_estimator {

        public:
            typedef KernelTraits Kernel;
            typedef InputPolygon Polygon;

            typename Kernel::Compute_squared_distance_2 squared_distance;

            using FT      = typename Kernel::FT;
            using Point_2 = typename Kernel::Point_2;
            using Line_2  = typename Kernel::Line_2;

            using Polygon_vertex_iterator = typename Polygon::const_iterator;

            Polygon_data_estimator(const Polygon &polygon) : 
            m_polygon(polygon) 
            { }

            FT compute_mean_distance_to_boundaries() const {
                
                const size_t num_edges = m_distances_to_edges.size();
                assert(num_edges > 0 && num_edges == m_polygon.size());

                FT mean = FT(0);
                for (size_t i = 0; i < num_edges; ++i) mean += m_distances_to_edges[i];
                mean /= static_cast<FT>(m_distances_to_edges.size());

                return mean;
            }

            FT compute_standard_deviation_on_distances_to_boundaries(const FT mean) const {

                const size_t num_edges = m_distances_to_edges.size();
                assert(num_edges > 0 && num_edges == m_polygon.size());

                FT stde = FT(0);
                for (size_t i = 0; i < num_edges; ++i) stde += (m_distances_to_edges[i] - mean) * (m_distances_to_edges[i] - mean);
                
                stde /= static_cast<FT>(m_distances_to_edges.size());
                stde  = static_cast<FT>(CGAL::sqrt(CGAL::to_double(stde)));

                return stde;
            }

            FT compute_minimum_distance_to_boundaries() const {

                const size_t num_edges = m_distances_to_edges.size();
                assert(num_edges > 0 && num_edges == m_polygon.size());

                FT min_distance = FT(1000000000000);
                for (size_t i = 0; i < num_edges; ++i)
                    min_distance = CGAL::min(min_distance, m_distances_to_edges[i]);

                return min_distance;
            }

            void compute_distances_to_boundaries(const Point_2 &query) {

                const size_t num_vertices = m_vertices.size();
                assert(num_vertices > 0);

                m_distances_to_edges.clear();
                m_distances_to_edges.resize(num_vertices);

                for (size_t i = 0; i < num_vertices; ++i) {
                    const size_t ip = (i + 1) % num_vertices;

                    const Point_2 &p1 = m_vertices[i];
                    const Point_2 &p2 = m_vertices[ip];

                    const Line_2 line(p1, p2);
                    const Point_2 projected = line.projection(query);

                    m_distances_to_edges[i] = static_cast<FT>(CGAL::sqrt(CGAL::to_double(CGAL::squared_distance(query, projected))));
                }
            }

            void compute_barycentre(Point_2 &barycentre) const {

                const size_t num_vertices = m_vertices.size();
                assert(num_vertices > 0);

                FT x = FT(0), y = FT(0);
                for (size_t i = 0; i < num_vertices; ++i) {
                    const Point_2 &p = m_vertices[i];

                    x += p.x();
                    y += p.y();
                }

                x /= static_cast<FT>(num_vertices);
                y /= static_cast<FT>(num_vertices);

                barycentre = Point_2(x, y);
            }

        private:
            const Polygon &m_polygon;

            std::vector<Point_2> m_vertices;
            std::vector<FT>      m_distances_to_edges;
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_POLYGON_DATA_ESTIMATOR_H