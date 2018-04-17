#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_ENVELOPE_INPUT_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_ENVELOPE_INPUT_H

// STL includes.
#include <map>
#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>

// CGAL includes.
#include <CGAL/utils.h>
#include <CGAL/number_utils.h>
#include <CGAL/Simple_cartesian.h>

// New CGAL includes.
#include <CGAL/IO/Color.h>
#include <CGAL/Level_of_detail_enum.h>
#include <CGAL/Buildings/Level_of_detail_diagonalize_traits.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class InputContainer, class InputBuildings, class EstimationStrategy>
		class Level_of_detail_building_envelope_input {
            
        public:
            typedef KernelTraits       Kernel;
            typedef InputContainer     Input;
            typedef InputBuildings     Buildings;
            typedef EstimationStrategy Strategy;

            typename Kernel::Compute_squared_distance_3 squared_distance_3;

            using FT         = typename Kernel::FT;
            using Point_2    = typename Kernel::Point_2;
            using Point_3    = typename Kernel::Point_3;
            using Plane_3    = typename Kernel::Plane_3;
            using Vector_3   = typename Kernel::Vector_3;
            using Triangle_3 = typename Kernel::Triangle_3;

            using Local_kernel       = CGAL::Simple_cartesian<double>;
            using Diagonalize_traits = CGAL::LOD::Eigen_diagonalize_traits_lod<double, 3>;

			using Point_3ft = Local_kernel::Point_3;
			using Plane_3ft = Local_kernel::Plane_3;

            using Building          = typename Strategy::Building;
            using Building_iterator = typename Buildings::iterator;

            using Index   = int;
			using Indices = std::vector<Index>;

            using Data           = typename Building::Data;
            using Data_triangle  = typename Building::Data_triangle;
            using Data_triangles = typename Building::Data_triangles;

            using Color = CGAL::Color;
            using Roofs = typename Building::Roofs;

            using Boundary = std::vector<Point_3>;
            using Points_3 = std::vector<Point_3>;

            Level_of_detail_building_envelope_input(const Input &input) :
            m_input(input), 
            m_strategy(input),
            m_alpha(-FT(1)),
            m_use_min_scale(false),
            m_big_value(FT(100000000000000)) { 

                assert(m_strategy.name() == "box");
            }

            void create(Buildings &buildings) {
                
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
            Strategy     m_strategy;
            
            FT         m_alpha;
            const bool m_use_min_scale;
            const FT   m_big_value;

            void process_building(Building &building) {
                
                const auto &shapes = building.shapes;
                if (shapes.size() == 0) {
                 
                    building.is_valid = false;
                    return;
                }

                building.clear_roofs();
				for (size_t i = 0; i < shapes.size(); ++i) {
                    
                    const Indices &indices = shapes[i];
                    process_roof(indices, building);
                }
            }

            void process_roof(const Indices &indices, Building &building) {
                
                assert(m_alpha > FT(0));
                assert(indices.size() > 2);

                Plane_3 plane;
                fit_plane_to_roof_points(indices, plane);

                Points_3 points;   
                project_points_onto_plane(indices, plane, points);

                m_strategy.set_alpha(m_alpha);
                m_strategy.estimate_roof(points, plane, building);

                create_envelope_input(building);
            }

            void fit_plane_to_roof_points(const Indices &indices, Plane_3 &plane) const {
                
                Point_3ft centroid;
                std::vector<Point_3ft> points;
                set_points_and_centroid(indices, points, centroid);

                Plane_3ft tmp_plane;
				CGAL::linear_least_squares_fitting_3(points.begin(), points.end(), tmp_plane, centroid, CGAL::Dimension_tag<0>(), Local_kernel(), Diagonalize_traits());
				plane = Plane_3(static_cast<FT>(tmp_plane.a()), static_cast<FT>(tmp_plane.b()), static_cast<FT>(tmp_plane.c()), static_cast<FT>(tmp_plane.d()));
            }

            void set_points_and_centroid(const Indices &indices, std::vector<Point_3ft> &points, Point_3ft &centroid) const {
                assert(indices.size() > 2);

                points.clear();
                points.resize(indices.size());

                double bx = 0.0, by = 0.0, bz = 0.0;
				for (size_t i = 0; i < indices.size(); ++i) {

					const Point_3 &p = m_input.point(indices[i]);

					const double x = CGAL::to_double(p.x());
					const double y = CGAL::to_double(p.y());
					const double z = CGAL::to_double(p.z());

					points[i] = Point_3ft(x, y, z);

                    bx += x;
                    by += y;
                    bz += z;
				}

                bx /= static_cast<double>(indices.size());
                by /= static_cast<double>(indices.size());
                bz /= static_cast<double>(indices.size());

                centroid = Point_3ft(bx, by, bz);
            }

            void project_points_onto_plane(const Indices &indices, const Plane_3 &plane, std::vector<Point_3> &points) const {
                assert(indices.size() > 2);

                points.clear();
                points.resize(indices.size());

                for (size_t i = 0; i < indices.size(); ++i) {			
					
                    const Point_3 &p = m_input.point(indices[i]);
					points[i] = plane.projection(p);
                }
            }

            void create_envelope_input(Building &building) const {
                
                size_t index = 0;
                
                building.clear_envelope_input();
                Data_triangles &triangles = building.envelope_input;
                
                add_roofs_triangles(building, triangles, index);
                add_walls_triangles(building, triangles, index);
            }

            void add_roofs_triangles(const Building &building, Data_triangles &triangles, size_t &index) const {

                const Roofs &roofs = building.roofs;
                for (size_t i = 0; i < roofs.size(); ++i) {
                    
                    Boundary roof = roofs[i].boundary;
                    assert(roof.size() == 4);

                    scale_roof(roof);
                    add_two_triangles(roof, building.color, false, triangles, index);
                }
            }

            void add_walls_triangles(const Building &building, Data_triangles &triangles, size_t &index) const {

                assert(!building.is_oriented);
                const auto &boundary = building.boundaries[0];

                Boundary wall(4);
                const FT height = building.height;

                assert(height > FT(0));

                for (size_t i = 0; i < boundary.size(); i += 2) {
                    const size_t ip = (i + 1) % boundary.size();

                    const Point_2 &v1 = boundary[i]->point();
                    const Point_2 &v2 = boundary[ip]->point();

                    wall[0] = Point_3(v1.x(), v1.y(), FT(0));
                    wall[1] = Point_3(v2.x(), v2.y(), FT(0));
                    wall[2] = Point_3(v2.x(), v2.y(), height);
                    wall[3] = Point_3(v1.x(), v1.y(), height);

                    add_two_triangles(wall, building.color, true, triangles, index);
                }
            }

            void scale_roof(Boundary &boundary) const {
                assert(boundary.size() == 4);

                Point_3 p1 = boundary[0];
                Point_3 p2 = boundary[1];
                Point_3 p3 = boundary[2];
                Point_3 p4 = boundary[3];

                const FT scale = compute_scale(boundary);

                Vector_3 v1 = Vector_3(p2, p1), v2 = Vector_3(p1, p2), v3 = Vector_3(p4, p3), v4 = Vector_3(p3, p4);
                Vector_3 v5 = Vector_3(p4, p1), v6 = Vector_3(p3, p2), v7 = Vector_3(p2, p3), v8 = Vector_3(p1, p4);
                
                normalize(v1); normalize(v2); normalize(v3); normalize(v4);
                normalize(v5); normalize(v6); normalize(v7); normalize(v8);

                move_point(p1, v1, scale); move_point(p2, v2, scale); move_point(p3, v3, scale); move_point(p4, v4, scale);
                move_point(p1, v5, scale); move_point(p2, v6, scale); move_point(p3, v7, scale); move_point(p4, v8, scale);

                boundary[0] = p1;
                boundary[1] = p2;
                boundary[2] = p3;
                boundary[3] = p4;
            }

            void normalize(Vector_3 &v) const {
                v /= static_cast<FT>(CGAL::sqrt(CGAL::to_double(v.squared_length())));
            }

            void move_point(Point_3 &p, const Vector_3 &v, const FT scale) const {
                assert(scale >= FT(0));

                const FT x = p.x() + scale * v.x();
                const FT y = p.y() + scale * v.y();
                const FT z = p.z() + scale * v.z();

                p = Point_3(x, y, z);
            }

            FT compute_scale(const Boundary &points) const {
                
                if (m_use_min_scale) return compute_min_scale(points);
                return compute_max_scale(points);
            }

            FT compute_min_scale(const Boundary &points) const {
                
                Point_3 query;
                compute_barycentre(points, query);

                return compute_min_distance(points, query) / FT(2);
            }

            void compute_barycentre(const Boundary &points, Point_3 &b) const {
                
                FT bx = FT(0), by = FT(0), bz = FT(0);
                for (size_t i = 0; i < points.size(); ++i) {

                    bx += points[i].x();
                    by += points[i].y();
                    bz += points[i].z();
                }

                bx /= static_cast<FT>(points.size());
                by /= static_cast<FT>(points.size());
                bz /= static_cast<FT>(points.size());

                b = Point_3(bx, by, bz);
            }

            FT compute_min_distance(const Boundary &points, const Point_3 &query) const {
                
                FT mindist = m_big_value;
                for (size_t i = 0; i < points.size(); ++i) {
                    const size_t ip = (i + 1) % points.size();
                    
                    const Point_3 &p1 = points[i];
                    const Point_3 &p2 = points[ip];

                    const FT x = (p1.x() + p2.x()) / FT(2);
                    const FT y = (p1.y() + p2.y()) / FT(2);
                    const FT z = (p1.z() + p2.z()) / FT(2);

                    const Point_3 point = Point_3(x, y, z);

                    const FT dist = static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_distance_3(point, query))));
                    mindist = CGAL::min(mindist, dist);
                }

                return mindist;
            }

            FT compute_max_scale(const Boundary &points) const {

                FT max_scale = -m_big_value;
                for (size_t i = 0; i < points.size(); ++i) {
                    const size_t ip = (i + 1) % points.size();
                 
                    const Point_3 &p1 = points[i];
                    const Point_3 &p2 = points[ip];

                    const FT dist = static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_distance_3(p1, p2))));
                    max_scale = CGAL::max(max_scale, dist);
                }

                return max_scale;
            }

            void add_two_triangles(const Boundary &boundary, const Color &color, const bool is_vertical, Data_triangles &triangles, size_t &index) const {
                assert(boundary.size() == 4);
                
                const Point_3 &p1 = boundary[0];
                const Point_3 &p2 = boundary[1];
                const Point_3 &p3 = boundary[2];
                const Point_3 &p4 = boundary[3];

                Data data1, data2;
                data1.color = color; data1.index = index; data1.is_vertical = is_vertical; ++index;
                data2.color = color; data2.index = index; data2.is_vertical = is_vertical; ++index;

                const Data_triangle tri1 = std::make_pair(Triangle_3(p1, p2, p3), data1);
                const Data_triangle tri2 = std::make_pair(Triangle_3(p3, p4, p1), data2);

                triangles.push_back(tri1);
                triangles.push_back(tri2);
            }
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_ENVELOPE_INPUT_H