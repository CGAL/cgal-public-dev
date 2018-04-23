#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_PARTITION_NAIVE_PLANE_ASSOCIATER_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_PARTITION_NAIVE_PLANE_ASSOCIATER_H

// STL includes.
#include <map>
#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>

// CGAL includes.
#include <CGAL/utils.h>

// New CGAL includes.
#include <CGAL/intersections.h>
#include <CGAL/Level_of_detail_enum.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputBuilding>
		class Level_of_detail_building_partition_naive_plane_associater {
            
        public:
            typedef KernelTraits  Kernel;
            typedef InputBuilding Building;

            using FT      = typename Kernel::FT;
            using Point_3 = typename Kernel::Point_3;
            using Plane_3 = typename Kernel::Plane_3;
            using Line_3  = typename Kernel::Line_3;

            using Roof              = typename Building::Roof;
            using Intersect         = typename Kernel::Intersect_3;
            using Envelope_input    = typename Building::Data_triangles;
            using Associated_planes = typename Roof::Associated_planes;
            
            using Boundary = std::vector<Point_3>;
            using Log      = CGAL::LOD::Mylog;

            Level_of_detail_building_partition_naive_plane_associater() :
            m_big_value(FT(100000000000000)) 
            { }

            void find_associated_planes(const Boundary &boundary, const Envelope_input &triangles, const FT reference_height, Associated_planes &associated_planes) const {
                associated_planes.clear();

                Point_3 barycentre;
                compute_barycentre(boundary, barycentre);
                
                Point_3 result;
                std::vector<size_t> indices;
                std::vector<FT> heights;

                for (size_t i = 0; i < triangles.size(); i += 2) {
                    if (triangles[i].second.is_vertical) continue;

                    const Point_3 &p1 = triangles[i].first.vertex(0);
			        const Point_3 &p2 = triangles[i].first.vertex(1);
				    const Point_3 &p3 = triangles[i].first.vertex(2);

                    const Plane_3 plane = Plane_3(p1, p2, p3);
                    intersect_line_with_plane(barycentre, plane, result);

                    const FT z = result.z();
                    if (z < reference_height) continue;

                    indices.push_back(i);
                    heights.push_back(z);
                }

                int final_index = -1;
                FT minz = m_big_value;

                for (size_t i = 0; i < indices.size(); ++i) {
                    if (heights[i] < minz) {

                        minz        = heights[i];
                        final_index = indices[i];
                    }
                }

                assert(final_index != -1);
                if (final_index != -1) associated_planes.push_back(static_cast<size_t>(final_index));
            }

        private:
            const FT m_big_value;

            void compute_barycentre(const Boundary &boundary, Point_3 &barycentre) const {

                FT x = FT(0), y = FT(0);
                for (size_t i = 0; i < boundary.size(); ++i) {
                    const Point_3 &p = boundary[i];

                    x += p.x();
                    y += p.y();
                }

                x /= static_cast<FT>(boundary.size());
                y /= static_cast<FT>(boundary.size());

                barycentre = Point_3(x, y, boundary[0].z());
            }

            void intersect_line_with_plane(const Point_3 &p, const Plane_3 &plane, Point_3 &r) const {
                
                const Point_3   q = Point_3(p.x(), p.y(), p.z() + FT(10));
                const Line_3 line = Line_3(p, q);

				typename CGAL::cpp11::result_of<Intersect(Line_3, Plane_3)>::type result = intersection(line, plane);
                r = boost::get<Point_3>(*result);
            }
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_PARTITION_NAIVE_PLANE_ASSOCIATER_H