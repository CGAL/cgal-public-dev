#ifndef CGAL_LEVEL_OF_DETAIL_BUILDINGS_HEIGHT_SETTER_H
#define CGAL_LEVEL_OF_DETAIL_BUILDINGS_HEIGHT_SETTER_H

namespace CGAL {

	namespace Level_of_detail {

        class Buildings_height_setter {
        
        public:
            template<class Building_height_map, class Buildings>
            void set_heights(const Building_height_map &building_height_map, Buildings &buildings) const {

                using Buildings_iterator = typename Buildings::iterator;
				for (Buildings_iterator bu_it = buildings.begin(); bu_it != buildings.end(); ++bu_it) {
                    
                    auto &building = *bu_it;
                    put(building_height_map, building);
                }
            }
        };

    } // Level_of_detail

} // CGAL

#endif // CGAL_LEVEL_OF_DETAIL_BUILDINGS_HEIGHT_SETTER_H