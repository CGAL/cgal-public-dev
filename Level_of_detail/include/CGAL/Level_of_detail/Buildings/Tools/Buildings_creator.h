#ifndef CGAL_LEVEL_OF_DETAIL_BUILDINGS_CREATOR_H
#define CGAL_LEVEL_OF_DETAIL_BUILDINGS_CREATOR_H

// STL includes.
#include <map>

namespace CGAL {

	namespace Level_of_detail {

        template<class InputKernel, class InputBuilding>
        class Buildings_creator {
        
        public:
            using Kernel   = InputKernel;
            using Building = InputBuilding;

            using Point_2 = typename Kernel::Point_2;
            
            using Buildings                = std::map<int, Building>;
            using Const_buildings_iterator = typename Buildings::const_iterator;

            template<class Triangulation, class Output>
            void create(const Triangulation &triangulation, Output &output) const {
                
                Buildings buildings;
                create_buildings(triangulation, buildings);
                create_output(buildings, output);
            }

        private:
            template<class Triangulation>
            void create_buildings(const Triangulation &triangulation, Buildings &buildings) const {
                
                buildings.clear();
                using Triangulation_faces_iterator = typename Triangulation::Finite_faces_iterator;

                Building new_building;
                for (Triangulation_faces_iterator tf_it = triangulation.finite_faces_begin(); tf_it != triangulation.finite_faces_end(); ++tf_it) {
                    
                    const int building_number = tf_it->info().group_number();
                    if (building_number < 0) continue;

                    if (is_new_building(building_number, buildings)) buildings[building_number] = new_building;
                    add_floor_face_to_building(tf_it, buildings.at(building_number));
                }
            }

            bool is_new_building(const int building_number, const Buildings &buildings) const {

                const Const_buildings_iterator cb_it = buildings.find(building_number);
                return cb_it == buildings.end();
            }

            template<class Triangulation_face_handle>
            void add_floor_face_to_building(const Triangulation_face_handle &face_handle, Building &building) const {

                const Point_2 &p1 = face_handle->vertex(0)->point();
                const Point_2 &p2 = face_handle->vertex(1)->point();
                const Point_2 &p3 = face_handle->vertex(2)->point();
			    
                building.add_floor_face(p1, p2, p3);
                building.add_floor_face_handle(face_handle);
            }

            template<class Output>
            void create_output(const Buildings &buildings, Output &output) const {
                
                output.clear();
                for (Const_buildings_iterator cb_it = buildings.begin(); cb_it != buildings.end(); ++cb_it) {
                    
                    const Building &building = (*cb_it).second;
                    output.push_back(building);
                }
            }
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDINGS_CREATOR_H