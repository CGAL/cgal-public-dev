#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_ROOF_FACE_VALIDATOR_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_ROOF_FACE_VALIDATOR_H

// STL includes.
#include <map>
#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/partition_2.h>
#include <CGAL/Partition_traits_2.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputBuilding>
		class Level_of_detail_building_roof_face_validator {
            
        public:
            typedef KernelTraits  Kernel;
            typedef InputBuilding Building;

            using FT         = typename Kernel::FT;
            using Point_2    = typename Kernel::Point_2;
            using Point_3    = typename Kernel::Point_3;
            using Triangle_2 = typename Kernel::Triangle_2;
            
            using Boundary = std::vector<Point_3>;

            using Partition_traits = CGAL::Partition_traits_2<Kernel>;
            using Polygon          = typename Partition_traits::Polygon_2;
            using Polygons         = std::vector<Polygon>;

            Level_of_detail_building_roof_face_validator() { }

            bool is_valid_roof_face(const Building &building, const Boundary &boundary) const {

                Point_2 query;
                find_query_point_using_barycentre(boundary, query);

                const auto &faces = building.faces;
                for (size_t i = 0; i < faces.size(); ++i) {
                        
                    const Point_2 &p1 = faces[i]->vertex(0)->point();
                    const Point_2 &p2 = faces[i]->vertex(1)->point();
                    const Point_2 &p3 = faces[i]->vertex(2)->point();

                    const Triangle_2 triangle = Triangle_2(p1, p2, p3);
                    if (triangle.has_on_bounded_side(query) || triangle.has_on_boundary(query)) return true;
                }
                return false;
            }

        private:
            
            void find_query_point_using_barycentre(const Boundary &boundary, Point_2 &query) const {

                FT x = FT(0), y = FT(0);
                for (size_t i = 0; i < boundary.size(); ++i) {
                    const Point_3 &p = boundary[i];
                    
                    x += p.x();
                    y += p.y();
                }

                x /= static_cast<FT>(boundary.size());
                y /= static_cast<FT>(boundary.size());

                query = Point_2(x, y);
            }

            void find_query_point_using_partition(const Boundary &boundary, Point_2 &query) const {
                
                const auto &points = boundary;
                assert(points.size() > 0);

                Polygon polygon;
                for (size_t i = 0; i < points.size(); ++i)
                    polygon.push_back(Point_2(points[i].x(), points[i].y()));

                if (orientation_2(polygon.vertices_begin(), polygon.vertices_end(), Partition_traits()) != CGAL::COUNTERCLOCKWISE) 
                    polygon.reverse_orientation();

                Polygons result;
                CGAL::approx_convex_partition_2(polygon.vertices_begin(), polygon.vertices_end(), std::back_inserter(result));

                FT x = FT(0), y = FT(0);
                for (auto vit = result[0].vertices_begin(); vit != result[0].vertices_end(); ++vit) {

                    x += (*vit).x();
                    y += (*vit).y();
                }

                const FT n = static_cast<FT>(std::distance(result[0].vertices_begin(), result[0].vertices_end()));

                x /= n;
                y /= n;

                query = Point_2(x, y);
            }
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_ROOF_FACE_VALIDATOR_H