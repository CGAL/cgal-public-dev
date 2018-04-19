#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_PARTITION_CREATOR_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_PARTITION_CREATOR_H

// STL includes.
#include <map>
#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>

// CGAL includes.
#include <CGAL/utils.h>

// New CGAL includes.
#include <CGAL/Regularizer/Level_of_detail_polygonizer_jean_philippe.h>
#include <CGAL/Container/Level_of_detail_container.h>
#include <CGAL/Level_of_detail_enum.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/intersections.h>
#include <CGAL/partition_2.h>
#include <CGAL/Mylog/Mylog.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputBuilding, class InputBuildings>
		class Level_of_detail_building_partition_creator {
            
        public:
            typedef KernelTraits   Kernel;
            typedef InputBuilding  Building;
            typedef InputBuildings Buildings;

            using FT         = typename Kernel::FT;
            using Point_2    = typename Kernel::Point_2;
            using Point_3    = typename Kernel::Point_3;
            using Triangle_2 = typename Kernel::Triangle_2;
            using Triangle_3 = typename Kernel::Triangle_3;
            using Segment_2  = typename Kernel::Segment_2;
            using Segment_3  = typename Kernel::Segment_3;
            using Plane_3    = typename Kernel::Plane_3;
            using Line_3     = typename Kernel::Line_3;

            using Intersect         = typename Kernel::Intersect_3;
            using Building_iterator = typename Buildings::iterator;

            using Roof              = typename Building::Roof;
            using Data              = typename Building::Data;
            using Envelope_input    = typename Building::Data_triangles;
            using Envelope_element  = typename Building::Data_triangle;
            using Partition_input   = typename Building::Partition_input;
            using Partition_element = typename Building::Partition_element;
            using Associated_planes = typename Roof::Associated_planes;
            
            using Boundary = std::vector<Point_3>;

            using Partition_traits = CGAL::Partition_traits_2<Kernel>;
            using Polygon          = typename Partition_traits::Polygon_2;
            using Polygons         = std::vector<Polygon>;
            
            using Log      = CGAL::LOD::Mylog;
            using Segments = std::vector<Segment_2>;

            using Data_structure = CGAL::LOD::Level_of_detail_container<Kernel>;
			using Polygonizer    = CGAL::LOD::Level_of_detail_polygonizer_jean_philippe<Kernel, Data_structure>;

            using DS_container  = typename Data_structure::Container;
            using DS_containers = typename Data_structure::Containers;
            
            using DS_polygon = typename DS_container::Polygon;
            using DS_polygon_vertices_iterator = typename DS_polygon::Vertex_const_iterator;

            Level_of_detail_building_partition_creator(const FT ground_height) :
            m_ground_height(ground_height), 
            m_num_intersections(1),
            m_min_face_width(-FT(1)),
            m_debug(false),
            m_big_value(FT(100000000000000))
            { }

            void create(Buildings &buildings) const {
                
                if (buildings.size() == 0) return;
				for (Building_iterator bit = buildings.begin(); bit != buildings.end(); ++bit) {

                    Building &building = bit->second;
					if (building.is_valid) process_building(building);
                }
            }

            void set_min_face_width(const FT new_value) {
                
                assert(new_value > FT(0));
                m_min_face_width = new_value;
            }

        private:
            const FT     m_ground_height;
            const size_t m_num_intersections;

            FT m_min_face_width;
            
            const bool m_debug;
            const FT m_big_value;
            
            void process_building(Building &building) const {
                
                const Partition_input &partition_input = building.partition_input;
                if (partition_input.size() < 3) {
                    
                    building.is_valid = false;
                    return;
                }

                Segments segments;
                set_segments(partition_input, segments);

                Data_structure data_structure;
                apply_polygonizer(segments, data_structure);

                update_roofs(data_structure, building);
            }

            void set_segments(const Partition_input &partition_input, Segments &segments) const {

                segments.clear();
                segments.resize(partition_input.size());

                for (size_t i = 0; i < partition_input.size(); ++i) {
                    
                    const Partition_element &element = partition_input[i];
                    const Segment_3 &segment = element.segment;

                    const Point_3 &source = segment.source();
                    const Point_3 &target = segment.target();

                    const Point_2 p1 = Point_2(source.x(), source.y());
                    const Point_2 p2 = Point_2(target.x(), target.y());

                    segments[i] = Segment_2(p1, p2);
                }
            }

            void apply_polygonizer(Segments &segments, Data_structure &data_structure) const {

                assert(m_num_intersections > 0);
                assert(m_min_face_width > FT(0));

				Polygonizer polygonizer;
				polygonizer.make_silent(true);

				polygonizer.set_number_of_intersections(m_num_intersections);
				polygonizer.set_min_face_width(m_min_face_width);

				polygonizer.polygonize(segments, data_structure);
            }

            void update_roofs(const Data_structure &data_structure, Building &building) const {

                Roof roof;
                building.clear_roofs();
                
                Boundary &boundary                   = roof.boundary;
                Associated_planes &associated_planes = roof.associated_planes;

                const DS_containers &containers = data_structure.containers();
                for (size_t i = 0; i < containers.size(); ++i) {
					
                    boundary.clear();
                    associated_planes.clear();

                    const DS_polygon &polygon = containers[i].polygon;
                    for (DS_polygon_vertices_iterator vit = polygon.vertices_begin(); vit != polygon.vertices_end(); ++vit) {
						const Point_2 &p = *vit;

                        const FT x = p.x();
                        const FT y = p.y();
                        const FT z = building.height + m_ground_height;

                        boundary.push_back(Point_3(x, y, z));
                    }

                    if (is_valid_roof_face(building, boundary)) {
                        const FT reference_height = building.roofs_min_height + FT(1) / FT(2);
                        
                        find_associated_planes(boundary, building.envelope_input, reference_height, associated_planes);
                        building.roofs.push_back(roof);
                    }
				}

                if (m_debug) save_polygons(data_structure);
            }

            bool is_valid_roof_face(const Building &building, const Boundary &boundary) const {

                Point_2 query;
                find_query_point(boundary, query);

                const auto &faces = building.faces;
                for (size_t i = 0; i < faces.size(); ++i) {
                        
                    const Point_2 &p1 = faces[i]->vertex(0)->point();
                    const Point_2 &p2 = faces[i]->vertex(1)->point();
                    const Point_2 &p3 = faces[i]->vertex(2)->point();

                    const Triangle_2 triangle = Triangle_2(p1, p2, p3);
                    if (triangle.has_on_bounded_side(query)) return true;
                }
                return false;
            }

            void find_query_point(const Boundary &boundary, Point_2 &query) const {
                
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

            void save_polygons(const Data_structure &data_structure) const {

                const DS_containers &containers = data_structure.containers();
                assert(containers.size() > 0);

                Log exporter;
                exporter.save_polygons<DS_containers, DS_polygon, Kernel>(containers, "tmp/lod_2/polygonizer_debug");
            }
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_PARTITION_CREATOR_H