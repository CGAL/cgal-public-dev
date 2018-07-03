#ifndef CGAL_LEVEL_OF_DETAIL_CONSTRAINED_TRIANGULATION_CREATOR_H
#define CGAL_LEVEL_OF_DETAIL_CONSTRAINED_TRIANGULATION_CREATOR_H

// STL includes.
#include <vector>

// LOD includes.
#include <CGAL/Level_of_detail/internal/utils.h>

namespace CGAL {

	namespace Level_of_detail {

        namespace LOD = CGAL::Level_of_detail;

		template<class InputKernel, class InputTriangulation>
		class Constrained_triangulation_creator {

		public:
			using Kernel        = InputKernel;
            using Triangulation = InputTriangulation;

            using Point_2 = typename Kernel::Point_2;

            using Triangulation_vertex_handle  = typename Triangulation::Vertex_handle;
            using Triangulation_vertex_handles = std::vector< std::vector<Triangulation_vertex_handle> >;
            using Triangulation_faces_iterator = typename Triangulation::Finite_faces_iterator;

            template<class Faces_range, class Point_map>
            void make_triangulation_with_info(const Faces_range &faces_range, const Point_map &point_map, Triangulation &triangulation) const {

				// Insert points.
                Triangulation_vertex_handles triangulation_vertex_handles;
                insert_points(faces_range, point_map, triangulation, triangulation_vertex_handles);

                // Insert constraints.
                insert_constraints(triangulation_vertex_handles, triangulation);

                // Update info.
				update_info(faces_range, triangulation);
            }

        private:
            template<class Input_faces_range, class Point_map>
            void insert_points(const Input_faces_range &input_faces_range, const Point_map &point_map, Triangulation &triangulation, Triangulation_vertex_handles &triangulation_vertex_handles) const {

                using Input_faces_iterator = typename Input_faces_range::const_iterator;

                triangulation.clear();
                triangulation_vertex_handles.clear();
                triangulation_vertex_handles.resize(input_faces_range.size());

                size_t i = 0;
                for (Input_faces_iterator if_it = input_faces_range.begin(); if_it != input_faces_range.end(); ++if_it, ++i) {
					
                    const auto &vertices = *if_it;
                    triangulation_vertex_handles[i].resize(vertices.size());

                    size_t j = 0;
					for (auto cv_it = vertices.begin(); cv_it != vertices.end(); ++cv_it, ++j) {
						
                        const Point_2 &point = get(point_map, *cv_it);
                        triangulation_vertex_handles[i][j] = triangulation.insert(point);
                    }
				}
            }

            void insert_constraints(const Triangulation_vertex_handles &triangulation_vertex_handles, Triangulation &triangulation) const {
                
                const size_t size_i = triangulation_vertex_handles.size();
				for (size_t i = 0; i < size_i; ++i) {

                    const size_t size_j = triangulation_vertex_handles[i].size();
					for (size_t j = 0; j < size_j; ++j) {
						const size_t jp = (j + 1) % size_j;
						
						if (triangulation_vertex_handles[i][j] != triangulation_vertex_handles[i][jp])
							triangulation.insert_constraint(triangulation_vertex_handles[i][j], triangulation_vertex_handles[i][jp]);
					}
				}
            }

            template<class Input_faces_range>
			void update_info(const Input_faces_range &input_faces_range, Triangulation &triangulation) const {
                using Input_faces_iterator = typename Input_faces_range::const_iterator;
                
                for (Triangulation_faces_iterator tf_it = triangulation.finite_faces_begin(); tf_it != triangulation.finite_faces_end(); ++tf_it) {
                    Point_2 barycentre = internal::barycenter<Kernel> (tf_it);
  
                    for (Input_faces_iterator if_it = input_faces_range.begin(); if_it != input_faces_range.end(); ++if_it) {
                        const auto &polygon = *if_it;

                        if (polygon.has_on_bounded_side(barycentre)) {
							tf_it->info().visibility_label() = polygon.visibility_label(); break;
						}
                    }
                }
			}
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_CONSTRAINED_TRIANGULATION_CREATOR_H
