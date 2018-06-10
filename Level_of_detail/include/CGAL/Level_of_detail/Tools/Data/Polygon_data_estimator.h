#ifndef CGAL_LEVEL_OF_DETAIL_POLYGON_DATA_ESTIMATOR_H
#define CGAL_LEVEL_OF_DETAIL_POLYGON_DATA_ESTIMATOR_H

// STL includes.
#include <vector>

// CGAL includes.
#include <CGAL/number_utils.h>
#include <CGAL/property_map.h>

// LOD includes.
#include <CGAL/Level_of_detail/Tools/Estimations/Barycentre_estimator.h>

namespace CGAL {

	namespace Level_of_detail {

        namespace LOD = CGAL::Level_of_detail;

		template<class InputKernel, class InputPolygon>
		class Polygon_data_estimator {

        public:
            using Kernel  = InputKernel;
            using Polygon = InputPolygon;

            typename Kernel::Compute_squared_distance_2 squared_distance_2;

            using FT      = typename Kernel::FT;
            using Point_2 = typename Kernel::Point_2;
            using Line_2  = typename Kernel::Line_2;

            using Const_vertex_iterator = typename Polygon::const_iterator;
            using Distances = std::vector<FT>;

            using Identity_point_map   = CGAL::Identity_property_map<Point_2>;
            using Barycentre_estimator = LOD::Barycentre_estimator<Kernel>;

            Polygon_data_estimator(const Polygon &polygon) : 
            m_polygon(polygon),
            m_big_value(FT(100000000000000))
            { }

            FT compute_mean_distance_to_boundaries() const {
                
                const size_t num_edges = m_distances_to_edges.size();
                CGAL_precondition(num_edges > 0 && num_edges == m_polygon.size());

                FT mean = FT(0);
                for (size_t i = 0; i < num_edges; ++i) mean += m_distances_to_edges[i];
                mean /= static_cast<FT>(m_distances_to_edges.size());

                return mean;
            }

            FT compute_standard_deviation_on_distances_to_boundaries(const FT mean) const {

                const size_t num_edges = m_distances_to_edges.size();
                CGAL_precondition(num_edges > 0 && num_edges == m_polygon.size());

                FT stde = FT(0);
                for (size_t i = 0; i < num_edges; ++i) 
                    stde += (m_distances_to_edges[i] - mean) * (m_distances_to_edges[i] - mean);
                
                stde /= static_cast<FT>(m_distances_to_edges.size());
                stde  = static_cast<FT>(CGAL::sqrt(CGAL::to_double(stde)));

                return stde;
            }

            FT compute_minimum_distance_to_boundaries() const {

                const size_t num_edges = m_distances_to_edges.size();
                CGAL_precondition(num_edges > 0 && num_edges == m_polygon.size());

                FT min_distance = m_big_value;
                for (size_t i = 0; i < num_edges; ++i)
                    min_distance = CGAL::min(min_distance, m_distances_to_edges[i]);

                return min_distance;
            }

            void compute_distances_to_boundaries(const Point_2 &query) {

                const size_t num_vertices = m_polygon.size();
                CGAL_precondition(num_vertices > 0);

                m_distances_to_edges.clear();
                m_distances_to_edges.resize(num_vertices);

                size_t i = 0;
                for (Const_vertex_iterator cv_it = m_polygon.begin(); cv_it != m_polygon.end(); ++i) {
                    const Point_2 &p1 = *cv_it;

                    const size_t ip = (i + 1) % num_vertices;
                    if (ip == 0) cv_it = m_polygon.begin();
                    else ++cv_it;
                    
                    const Point_2 &p2 = *cv_it;
                    if (ip == 0) cv_it = m_polygon.end();

                    const Line_2 line(p1, p2);
                    const Point_2 projected = line.projection(query);

                    m_distances_to_edges[i] = static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_distance_2(query, projected))));
                }
            }

            void compute_barycentre(Point_2 &barycentre) const {

                Identity_point_map identity_point_map;
                const Barycentre_estimator barycentre_estimator;
                barycentre_estimator.compute_barycentre_2(m_polygon, identity_point_map, barycentre);
            }

        private:
            const Polygon &m_polygon;
            
            Distances m_distances_to_edges;
            const FT  m_big_value;
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_POLYGON_DATA_ESTIMATOR_H