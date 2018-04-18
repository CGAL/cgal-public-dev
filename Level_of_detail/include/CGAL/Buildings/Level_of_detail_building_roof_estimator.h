#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_ROOF_ESTIMATOR_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_ROOF_ESTIMATOR_H

// STL includes.
#include <map>
#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/intersections.h>

// New CGAL includes.
#include <CGAL/Level_of_detail_enum.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputBuilding, class InputBuildings>
		class Level_of_detail_building_roof_estimator {
            
        public:
            typedef KernelTraits   Kernel;
            typedef InputBuilding  Building;
            typedef InputBuildings Buildings;

            typedef typename Kernel::Intersect_3 Intersect;

            using FT      = typename Kernel::FT;
            using Line_3  = typename Kernel::Line_3;
            using Point_3 = typename Kernel::Point_3;
            using Plane_3 = typename Kernel::Plane_3;

            using Building_iterator = typename Buildings::iterator;

            using Roof              = typename Building::Roof;
            using Roofs             = typename Building::Roofs;
            using Envelope_input    = typename Building::Data_triangles;

            using Associated_planes = typename Roof::Associated_planes; 
            using Boundary          = typename Roof::Roof_boundary;

            Level_of_detail_building_roof_estimator(const FT ground_height) :
            m_ground_height(ground_height),
            m_big_value(FT(100000000000000)) { }

            void estimate(Buildings &buildings) const {
                
                if (buildings.size() == 0) return;
				for (Building_iterator bit = buildings.begin(); bit != buildings.end(); ++bit) {

                    auto &building = bit->second;
					if (building.is_valid) process_building(building);
                }
            }

        private:
            const FT m_ground_height;
            const FT m_big_value;

            void process_building(Building &building) const {
                
                Roofs &roofs = building.roofs;
                if (roofs.size() == 0) {

                    building.is_valid = false;
                    return;
                }

                for (size_t i = 0; i < roofs.size(); ++i) {
                    
                    Roof &roof = roofs[i];
                    process_roof(building, roof);
                }
            }

            void process_roof(const Building &building, Roof &roof) const {

                Boundary &points = roof.boundary;
                assert(points.size() > 2);

                const Associated_planes &associated_planes = roof.associated_planes;
                if (associated_planes.size() == 0) roof.is_valid = false;
                
                for (size_t i = 0; i < points.size(); ++i) {
                    Point_3 &p = points[i];

                    const FT x = p.x();
                    const FT y = p.y();
                          FT z = estimate_z_coordinate(building.envelope_input, associated_planes, p);

                    if (!roof.is_valid) z += building.height;
                    p = Point_3(x, y, z);
                }
            }
            
            FT estimate_z_coordinate(const Envelope_input &triangles, const Associated_planes &associated_planes, const Point_3 &p) const {
                
                if (associated_planes.size() == 0) 
                    return m_ground_height;

                assert(associated_planes.size() == 1);

                if (associated_planes.size() > 1) 
                    std::cerr << std::endl << "WARNING: Too many associated planes!" << std::endl << std::endl;

                const size_t index = associated_planes[0];
                if (triangles[index].second.is_vertical) return m_ground_height;

                const Point_3 &p1 = triangles[index].first.vertex(0);
			    const Point_3 &p2 = triangles[index].first.vertex(1);
				const Point_3 &p3 = triangles[index].first.vertex(2);

                const Plane_3 plane = Plane_3(p1, p2, p3);
                return intersect_line_with_plane(p, plane);
            }

            FT intersect_line_with_plane(const Point_3 &p, const Plane_3 &plane) const {
                
                const Point_3   q = Point_3(p.x(), p.y(), m_ground_height);
                const Line_3 line = Line_3(p, q);

				typename CGAL::cpp11::result_of<Intersect(Line_3, Plane_3)>::type result = intersection(line, plane);
                const Point_3 r = boost::get<Point_3>(*result);

				return r.z();
            }
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_ROOF_ESTIMATOR_H