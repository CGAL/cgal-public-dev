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
#include <CGAL/Region_growing/Level_of_detail_linear_region_growing.h>

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
            using Segment_2  = typename Kernel::Segment_2;
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

            using Segment_region_growing = CGAL::LOD::Level_of_detail_linear_region_growing<Kernel>;
            
            using States   = typename Segment_region_growing::States;
            using Segments = typename Segment_region_growing::Segments;
            using Regions  = typename Segment_region_growing::Output;

            Level_of_detail_building_partition_input(const Input &input) :
            m_input(input),
            m_big_value(FT(100000000000000)),
            m_alpha(-FT(1)),
            m_tolerance(FT(1) / FT(100000)),
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

            FT       m_alpha;
            const FT m_tolerance;

            Envelope_input         m_envelope_input;
            Segment_region_growing m_region_growing;

            void create_envelope_input(Buildings &buildings) {

                m_envelope_input.set_alpha(m_alpha);
                m_envelope_input.create(buildings);
            }

            void process_building(Building &building) const {
                
                set_roofs_min_and_max_height(building);

                building.clear_partition_input();
                building.clear_partition_segments();

                Partition_input &input = building.partition_input;
                Segments &segments     = building.partition_segments;

                add_walls_segments(building, input);
                add_roofs_segments(building, input);

                States states;
                set_segments(input, segments, states);
                
                merge_segments(states, segments);
                std::cout << std::endl;
            }

            void set_roofs_min_and_max_height(Building &building) const {
                
                FT minz = m_big_value, maxz = -m_big_value;
                const auto &shapes = building.shapes;

				for (size_t i = 0; i < shapes.size(); ++i) {
                    const Indices &indices = shapes[i];

                    for (size_t j = 0; j < indices.size(); ++j) {    
                        const Point_3 &p = m_input.point(indices[j]);
                        
                        minz = CGAL::min(minz, p.z());
                        maxz = CGAL::max(maxz, p.z());
                    }
                }

                building.roofs_min_height = minz;
                building.roofs_max_height = maxz;
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

                    add_segment(p1, p2, true, input);
                }
            }

            void add_segment(const Point_3 &p1, const Point_3 &p2, const bool state, Partition_input &input) const {

                Partition_element element;
                
                element.to_be_used = state;
                element.segment    = Segment_3(p1, p2);

                input.push_back(element);
            }

            void add_roofs_segments(const Building &building, Partition_input &input) const {

                const Data_triangles &triangles = building.envelope_input;
                const size_t num_triangles = triangles.size();
                
                if (num_triangles < 4) return;
                assert(num_triangles % 2 == 0);

                for (size_t i = 0; i < num_triangles; ++i) {
                    for (size_t j = i; j < num_triangles; ++j) {

                        if (!triangles[i].second.is_vertical && !triangles[j].second.is_vertical && i != j) {
                            
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
                            add_segment(source, target, true, input);
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

            void set_segments(const Partition_input &partition_input, Segments &segments, States &states) const {

                segments.clear();
                segments.resize(partition_input.size());

                states.clear();
                states.resize(partition_input.size(), false);

                for (size_t i = 0; i < partition_input.size(); ++i) {
                    
                    const Partition_element &element = partition_input[i];
                    const Segment_3 &segment = element.segment;

                    const Point_3 &source = segment.source();
                    const Point_3 &target = segment.target();

                    const Point_2 p1 = Point_2(source.x(), source.y());
                    const Point_2 p2 = Point_2(target.x(), target.y());

                    segments[i] = Segment_2(p1, p2);
                    if (element.to_be_used) states[i] = true;
                }
            }

            void merge_segments(const States &states, Segments &segments) const {
                
                // Find groups of segments.
                Regions result;
                m_region_growing.find_connected_segments(segments, states, result);

                std::cout << std::endl << "num regions: " << result.size() << std::endl;

                // Merge segments.
                Segments new_segments;
                for (size_t i = 0; i < segments.size(); ++i)
                    if (!states[i]) add_separate_segment(segments[i], new_segments);

                std::cout << "num segments before: " << new_segments.size() << std::endl;

                for (size_t i = 0; i < result.size(); ++i)
                    add_merged_segment(segments, result[i], new_segments);

                segments = new_segments;

                std::cout << "num segments after: " << new_segments.size() << std::endl;
            }

            void add_separate_segment(const Segment_2 &segment, Segments &new_segments) const {
                
                Segment_2 new_segment;
                get_shortened_segment(segment, new_segment);

                if (new_segment.squared_length() > m_tolerance * m_tolerance)
                    new_segments.push_back(new_segment);
            }

            void get_shortened_segment(const Segment_2 &segment, Segment_2 &new_segment) const {
                
                const Point_2 &source = segment.source();
                const Point_2 &target = segment.target();

                const FT b1 = FT(9) / FT(10);
                const FT b2 = FT(1) - b1;

                const FT x1 = b1 * source.x() + b2 * target.x();
                const FT y1 = b1 * source.y() + b2 * target.y();

                const FT x2 = b2 * source.x() + b1 * target.x();
                const FT y2 = b2 * source.y() + b1 * target.y();

                new_segment = Segment_2(Point_2(x1, y1), Point_2(x2, y2));
            }

            void add_merged_segment(const Segments &segments, const Indices &indices, Segments &new_segments) const {

                Segment_2 merged_segment;
                if (get_merged_segment(segments, indices, merged_segment))
                    add_separate_segment(merged_segment, new_segments);
            }

            bool get_merged_segment(const Segments &segments, const Indices &indices, Segment_2 &merged_segment) const {

                size_t count = 0;
                std::vector<Point_2> input_points(indices.size() * 2);

                for (size_t i = 0; i < indices.size(); ++i) {
                    
                    input_points[count++] = segments[indices[i]].source();
                    input_points[count++] = segments[indices[i]].target();
                }

                std::vector<Point_2> output_points;
                for (size_t i = 0; i < input_points.size(); ++i) {
                    
                    if (!exists(i, input_points[i], input_points))
                        output_points.push_back(input_points[i]);
                }

                if (output_points.size() == 2) {
                    
                    merged_segment = Segment_2(output_points[0], output_points[1]);
                    return true;
                }
                return false;
            }

            bool exists(const size_t index, const Point_2 &point, const std::vector<Point_2> &input_points) const {

                for (size_t i = 0; i < input_points.size(); ++i)
                    if (i != index && are_equal(point, input_points[i])) 
                        return true;
                return false;
            }

            bool are_equal(const Point_2 &p, const Point_2 &q) const {

                const FT eps = m_tolerance;
                if (CGAL::abs(p.x() - q.x()) < eps && CGAL::abs(p.y() - q.y()) < eps) return true;
                return false;
            }
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_PARTITION_INPUT_H