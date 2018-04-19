#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_PARTITION_INPUT_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_PARTITION_INPUT_H

#if defined(WIN32) || defined(_WIN32) 
#define PSR "\\" 
#else 
#define PSR "/" 
#endif 

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
#include <CGAL/IO/Color.h>
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Level_of_detail_enum.h>
#include <CGAL/Buildings/Level_of_detail_building_envelope_input.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputContainer, class InputBuildings, class EstimationStrategy>
		class Level_of_detail_building_partition_input {
            
        public:
            typedef KernelTraits       Kernel;
            typedef InputContainer     Input;
            typedef InputBuildings     Buildings;
            typedef EstimationStrategy Strategy;

            using FT         = typename Kernel::FT;
            using Point_2    = typename Kernel::Point_2;
            using Point_3    = typename Kernel::Point_3;
            using Triangle_3 = typename Kernel::Triangle_3;
            using Segment_3  = typename Kernel::Segment_3;

            using Building          = typename Strategy::Building;
            using Building_iterator = typename Buildings::iterator;

            using Data           = typename Building::Data;
            using Data_triangle  = typename Building::Data_triangle;
            using Data_triangles = typename Building::Data_triangles;

            using Partition_input   = typename Building::Partition_input;
            using Partition_element = typename Building::Partition_element;

            using Color = CGAL::Color;
            using Roofs = typename Building::Roofs;

            using Boundary = std::vector<Point_3>;
            using Points_3 = std::vector<Point_3>;

            using Log         = CGAL::LOD::Mylog;
            using Intersect_3 = typename Kernel::Intersect_3;

            using Envelope_input = CGAL::LOD::Level_of_detail_building_envelope_input<Kernel, Input, Buildings, Strategy>;
            
            using Index   = int;
			using Indices = std::vector<Index>;

            Level_of_detail_building_partition_input(const Input &input) :
            m_input(input),
            m_big_value(FT(100000000000000)),
            m_alpha(-FT(1)),
            m_envelope_input(input)
            { }

            void create(Buildings &buildings) {
                create_envelope_input(buildings);

                if (buildings.size() == 0) return;
				for (Building_iterator bit = buildings.begin(); bit != buildings.end(); ++bit) {

                    Building &building = bit->second;
					if (building.is_valid) process_building(building);
                }
            }

            void set_alpha(const FT new_value) {
                
                assert(new_value > FT(0));
                m_alpha = new_value;
            }

        private:
            const Input &m_input;
            const FT m_big_value;

            FT m_alpha;
            Envelope_input m_envelope_input;

            void create_envelope_input(Buildings &buildings) {

                m_envelope_input.set_alpha(m_alpha);
                m_envelope_input.create(buildings);
            }

            void process_building(Building &building) const {
                set_roofs_min_height(building);

                building.clear_partition_input();
                Partition_input &input = building.partition_input;

                add_walls_segments(building, input);
                add_roofs_segments(building, input);
            }

            void set_roofs_min_height(Building &building) const {
                
                FT minz = m_big_value;
                const auto &shapes = building.shapes;

				for (size_t i = 0; i < shapes.size(); ++i) {
                    
                    const Indices &indices = shapes[i];
                    for (size_t j = 0; j < indices.size(); ++j) {
                        
                        const Point_3 &p = m_input.point(indices[j]);
                        minz = CGAL::min(minz, p.z());
                    }
                }
                building.roofs_min_height = minz;
            }

            void add_walls_segments(const Building &building, Partition_input &input) const {

                assert(!building.is_oriented);
                const auto &boundary = building.boundaries[0];

                for (size_t i = 0; i < boundary.size(); i += 2) {
                    const size_t ip = (i + 1) % boundary.size();

                    const Point_2 &v1 = boundary[i]->point();
                    const Point_2 &v2 = boundary[ip]->point();

                    const Point_3 p1 = Point_3(v1.x(), v1.y(), FT(0));
                    const Point_3 p2 = Point_3(v2.x(), v2.y(), FT(0));

                    add_segment(p1, p2, input);
                }
            }

            void add_segment(const Point_3 &p1, const Point_3 &p2, Partition_input &input) const {

                Partition_element element;
                element.segment = Segment_3(p1, p2);

                input.push_back(element);
            }

            void add_roofs_segments(const Building &building, Partition_input &input) const {

                const Data_triangles &triangles = building.envelope_input;
                const size_t num_triangles = triangles.size();
                
                if (num_triangles < 4) return;
                assert(num_triangles % 2 == 0);

                for (size_t i = 0; i < num_triangles; ++i) {
                    for (size_t j = 0; j < num_triangles; ++j) {

                        if (!triangles[i].second.is_vertical && !triangles[j].second.is_vertical) {
                            
                            const Triangle_3 &tri1 = triangles[i].first;
                            const Triangle_3 &tri2 = triangles[j].first;

                            intersect_two_triangles(tri1, tri2, input);
                        }
                    }
                }
            }

            void intersect_two_triangles(const Triangle_3 &tri1, const Triangle_3 &tri2, Partition_input &input) const {

                typename CGAL::cpp11::result_of<Intersect_3(Triangle_3, Triangle_3)>::type result = intersection(tri1, tri2);
                if (result) {

                    if (const Segment_3* segment = boost::get<Segment_3>(&*result)) {
                     
                        const Point_3 &source = (*segment).source();
                        const Point_3 &target = (*segment).target();

                        if (is_valid_point(source) && is_valid_point(target))
                            add_segment(source, target, input);
                    }
                }
            }

            bool is_valid_point(const Point_3 &p) const {

                const double x = CGAL::to_double(p.x());
                const double y = CGAL::to_double(p.y());
                const double z = CGAL::to_double(p.z());

                bool is_valid = true;
                if (!std::isfinite(x)) is_valid = false;
                if (!std::isfinite(y)) is_valid = false;
                if (!std::isfinite(z)) is_valid = false;

                return is_valid;
            }
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_PARTITION_INPUT_H