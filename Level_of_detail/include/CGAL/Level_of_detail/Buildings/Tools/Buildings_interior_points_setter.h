#ifndef CGAL_LEVEL_OF_DETAIL_BUILDINGS_INTERIOR_POINTS_SETTER_H
#define CGAL_LEVEL_OF_DETAIL_BUILDINGS_INTERIOR_POINTS_SETTER_H

namespace CGAL {

	namespace Level_of_detail {

        template<class InputKernel, class InputTriangulation>
        class Buildings_interior_points_setter {
        
        public:
            using Kernel        = InputKernel;
            using Triangulation = InputTriangulation;

            using Point_2 = typename Kernel::Point_2;
            
            using Triangulation_locate_type    = typename Triangulation::Locate_type;
            using Triangulation_face_handle    = typename Triangulation::Face_handle;
            using Triangulation_faces_iterator = typename Triangulation::Finite_faces_iterator;

            template<class Elements, class Point_map>
            void set_points(const Elements &elements, const Point_map &point_map, Triangulation &triangulation) const {
                
                clear_triangulation_elements(triangulation);
                set_elements(elements, point_map, triangulation);
            }

        private:
            void clear_triangulation_elements(Triangulation &triangulation) const {
                
                for (Triangulation_faces_iterator tf_it = triangulation.finite_faces_begin(); tf_it != triangulation.finite_faces_end(); ++tf_it)
                    tf_it->info().elements().clear();
            }

            template<class Elements, class Point_map>
            void set_elements(const Elements &elements, const Point_map &point_map, Triangulation &triangulation) const {

                CGAL_precondition(elements.size() > 0);
                using Const_elements_iterator = typename Elements::const_iterator;

                for (Const_elements_iterator ce_it = elements.begin(); ce_it != elements.end(); ++ce_it) {
                    const Point_2 &query = get(point_map, *ce_it);

                    Triangulation_locate_type locate_type; int locate_index_stub;
                    const Triangulation_face_handle face_handle = triangulation.locate(query, locate_type, locate_index_stub);

					if (locate_type == Triangulation::FACE || locate_type == Triangulation::EDGE || locate_type == Triangulation::VERTEX)
						face_handle->info().add_element(*ce_it);
                }
            }
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDINGS_INTERIOR_POINTS_SETTER_H