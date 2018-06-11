#ifndef CGAL_LEVEL_OF_DETAIL_BUILDINGS_SETTER_H
#define CGAL_LEVEL_OF_DETAIL_BUILDINGS_SETTER_H

namespace CGAL {

	namespace Level_of_detail {

        class Buildings_setter {
        
		public:
            template<class Face_to_building_map, class Triangulation>
			void set_buildings(const Face_to_building_map &face_to_building_map, Triangulation &triangulation) const {
				
				using Triangulation_faces_iterator = typename Triangulation::Finite_faces_iterator;
				using Triangulation_face_handle    = typename Triangulation::Face_handle;

				for (Triangulation_faces_iterator face = triangulation.finite_faces_begin(); face != triangulation.finite_faces_end(); ++face) {
					
					Triangulation_face_handle face_handle = static_cast<Triangulation_face_handle>(face);
					put(face_to_building_map, face_handle);
				}
			}
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDINGS_SETTER_H