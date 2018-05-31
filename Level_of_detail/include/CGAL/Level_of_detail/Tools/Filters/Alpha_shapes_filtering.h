#ifndef CGAL_LEVEL_OF_DETAIL_ALPHA_SHAPES_FILTERING_H
#define CGAL_LEVEL_OF_DETAIL_ALPHA_SHAPES_FILTERING_H

// CGAL includes.
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

namespace CGAL {

	namespace Level_of_detail {

		template<class InputKernel, class PointIdentifier>
		class Alpha_shapes_filtering {

        public:
            using Kernel           = InputKernel;
            using Point_identifier = PointIdentifier;
            
            using FT      = typename Kernel::FT;
            using Point_2 = typename Kernel::Point_2;

            using Vbi = CGAL::Triangulation_vertex_base_with_info_2<Point_identifier, Kernel>;
            using Vb  = CGAL::Alpha_shape_vertex_base_2<Kernel, Vbi>;
            using Fb  = CGAL::Alpha_shape_face_base_2<Kernel>;
            using Tds = CGAL::Triangulation_data_structure_2<Vb, Fb>;
            
            using Triangulation_2 = CGAL::Delaunay_triangulation_2<Kernel, Tds>;
            using Alpha_shape_2   = CGAL::Alpha_shape_2<Triangulation_2>;

            using Alpha_vertex_iterator = typename Alpha_shape_2::Alpha_shape_vertices_iterator;
            using Vertex_handle         = typename Triangulation_2::Vertex_handle;

            Alpha_shapes_filtering(const FT alpha) : 
            m_alpha(alpha) 
            { }

            template<class Elements, class Point_map, class Output>
            void add_points(const Elements &elements, const Point_map &point_map, Output &output) const {
                CGAL_precondition(elements.size() > 2);

				Triangulation_2 triangulation;
                create_triangulation(elements, point_map, triangulation);

				CGAL_precondition(m_alpha > FT(0));
				Alpha_shape_2 alpha_shape(triangulation, m_alpha, Alpha_shape_2::GENERAL);

				for (Alpha_vertex_iterator av_it = alpha_shape.alpha_shape_vertices_begin(); av_it != alpha_shape.alpha_shape_vertices_end(); ++av_it)
					output.push_back((*av_it)->info());
            }

        private:
            const FT m_alpha;

            template<class Elements, class Point_map>
            void create_triangulation(const Elements &elements, const Point_map &point_map, Triangulation_2 &triangulation) const {
                triangulation.clear();

                for (typename Elements::const_iterator element = elements.begin(); element != elements.end(); ++element) {
                    const Point_2 &point = get(point_map, *element);

                    Vertex_handle vertex_handle = triangulation.insert(point);
                    vertex_handle->info() = *element;
                }
            }
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_ALPHA_SHAPES_FILTERING_H