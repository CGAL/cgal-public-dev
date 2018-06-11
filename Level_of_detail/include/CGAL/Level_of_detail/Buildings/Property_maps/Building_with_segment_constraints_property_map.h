#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_WITH_SEGMENT_CONSTRAINTS_PROPERTY_MAP_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_WITH_SEGMENT_CONSTRAINTS_PROPERTY_MAP_H

namespace CGAL {

	namespace Level_of_detail {

		template<class InputKernel, class InputTriangulation>
		class Building_with_segment_constraints_property_map {
			
        public:
            using Kernel        = InputKernel;
            using Triangulation = InputTriangulation;

            using FT = typename Kernel::FT;
            
            using Triangulation_face_handle = typename Triangulation::Face_handle; 

            template<class Elements, class Segment_map>
            Building_with_segment_constraints_property_map(
                const Triangulation &triangulation, 
                const Elements &elements,
                const Segment_map &segment_map,
                const FT tolerance) :
			m_triangulation(triangulation),
            m_tolerance(tolerance) { 

            }

            using key_type = Triangulation_face_handle;
            using Self     = Building_with_segment_constraints_property_map<Kernel, Triangulation>;

            friend void put(const Self &self, key_type &face_handle) {
                face_handle->info().group_number() = 0;
            }

        private:
            const Triangulation &m_triangulation;
            const FT m_tolerance;
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_WITH_SEGMENT_CONSTRAINTS_PROPERTY_MAP_H