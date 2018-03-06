#ifndef CGAL_LEVEL_OF_DETAIL_POLYGON_SAMPLER_2_H
#define CGAL_LEVEL_OF_DETAIL_POLYGON_SAMPLER_2_H

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
#include <CGAL/Delaunay_triangulation_2.h>

// New CGAL includes.
#include <CGAL/Utils/Level_of_detail_uniform_sample_generator.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputPolygon>
		class Level_of_detail_polygon_sampler_2 {

        public:
            typedef KernelTraits Kernel;
            typedef InputPolygon Polygon;

            using FT         = typename Kernel::FT;
            using Point_2    = typename Kernel::Point_2;

            using Polygon_vertex_iterator = typename Polygon::Vertex_const_iterator;
            using Delaunay_triangulation  = CGAL::Delaunay_triangulation_2<Kernel>;

            using Faces_iterator   = typename Delaunay_triangulation::Finite_faces_iterator;
            using Sample_generator = CGAL::LOD::Uniform_sample_generator<Kernel>;

            Level_of_detail_polygon_sampler_2(const Polygon &polygon): m_polygon(polygon), m_num_sub_steps(0) { 
                m_dt.insert(m_polygon.vertices_begin(), m_polygon.vertices_end());
            }

            void create_samples(std::vector<Point_2> &samples) {

                samples.clear();
                m_sample_generator.set_number_of_samples(m_num_sub_steps);

                std::vector<Point_2> new_samples;
                for (Faces_iterator fit = m_dt.finite_faces_begin(); fit != m_dt.finite_faces_end(); ++fit) {

                    const Point_2 &p1 = fit->vertex(0)->point();
                    const Point_2 &p2 = fit->vertex(1)->point();
                    const Point_2 &p3 = fit->vertex(2)->point();

                    m_sample_generator.create_uniform_subdivision_samples(p1, p2, p3, new_samples);
                    for (size_t i = 0; i < new_samples.size(); ++i) samples.push_back(new_samples[i]);
                }
            }

            void set_number_of_subdivision_steps(const size_t new_value) {
                
                assert(new_value >= 0);
                m_num_sub_steps = new_value;
            }

        private:
            const Polygon &m_polygon;

            Delaunay_triangulation m_dt;
            Sample_generator       m_sample_generator;

            size_t m_num_sub_steps;
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_POLYGON_SAMPLER_2_H