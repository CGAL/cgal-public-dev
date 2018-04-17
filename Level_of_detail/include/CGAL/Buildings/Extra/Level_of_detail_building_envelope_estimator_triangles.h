#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_ENVELOPE_ESTIMATOR_TRIANGLES_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_ENVELOPE_ESTIMATOR_TRIANGLES_H

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

#include <CGAL/envelope_3.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Env_triangle_traits_3.h>
#include <CGAL/Env_surface_data_traits_3.h>
#include <CGAL/Partition_traits_2.h>
#include <CGAL/intersections.h>
#include <CGAL/partition_2.h>


// New CGAL includes.
#include <CGAL/IO/Color.h>
#include <CGAL/Mylog/Mylog.h>
#include <CGAL/Level_of_detail_enum.h>
#include <CGAL/Buildings/Level_of_detail_diagonalize_traits.h>

namespace CGAL {

	namespace LOD {

		template<class KernelTraits, class ContainerInput, class BuildingsInput, class EstimationStrategy>
		class Level_of_detail_building_envelope_estimator_triangles {
            
        public:
            typedef KernelTraits       Kernel;
            typedef ContainerInput     Input;
            typedef BuildingsInput     Buildings;
            typedef EstimationStrategy Strategy;

            typename Kernel::Compute_squared_distance_2 squared_distance_2;
            typename Kernel::Compute_squared_distance_3 squared_distance_3;

            typename Kernel::Compute_squared_length_3 		  squared_length;
			typename Kernel::Compute_scalar_product_3 		  dot_product;
			typename Kernel::Construct_cross_product_vector_3 cross_product;

            typedef typename Kernel::Intersect_3 Intersect;

            using FT         = typename Kernel::FT;
            using Line_3     = typename Kernel::Line_3;
            using Point_2    = typename Kernel::Point_2;
            using Point_3    = typename Kernel::Point_3;
            using Plane_3    = typename Kernel::Plane_3;
            using Vector_3   = typename Kernel::Vector_3;
            using Triangle_2 = typename Kernel::Triangle_2;

            using Local_kernel       = CGAL::Simple_cartesian<double>;
            using Diagonalize_traits = CGAL::LOD::Eigen_diagonalize_traits_lod<double, 3>;

			using Point_3ft = Local_kernel::Point_3;
			using Plane_3ft = Local_kernel::Plane_3;
			using Line_3ft  = Local_kernel::Line_3;

            using Building          = typename Strategy::Building;
            using Building_iterator = typename Buildings::iterator;

            using Index   = int;
			using Indices = std::vector<Index>;

            struct Data {
                size_t index;
                CGAL::Color color;
            };

            using Exact_type     = CGAL::Exact_rational;
            using Exact_kernel   = CGAL::Cartesian<Exact_type>;
            using Env_traits     = CGAL::Env_triangle_traits_3<Exact_kernel>;
            using Exact_point    = typename Exact_kernel::Point_3;
            using Exact_triangle = typename Env_traits::Surface_3;
            using Data_traits    = CGAL::Env_surface_data_traits_3<Env_traits, Data>;
            using Data_triangle  = typename Data_traits::Surface_3;
            using Env_diagram    = CGAL::Envelope_diagram_2<Data_traits>;

            using Face_iterator       = typename Env_diagram::Face_const_iterator;
            using Halfedge_circulator = typename Env_diagram::Ccb_halfedge_const_circulator;
            using Surface_iterator    = typename Env_diagram::Surface_const_iterator;

            using Data_triangles = std::vector<Data_triangle>;

            using Roofs = typename Building::Roofs;
            using Roof  = typename Building::Roof;

            using Boundary = std::vector<Point_3>;
            using Points_3 = std::vector<Point_3>;
            
            using Log   = CGAL::LOD::Mylog;
            using Color = CGAL::Color;

            using Partition_traits = CGAL::Partition_traits_2<Kernel>;
            using Polygon          = typename Partition_traits::Polygon_2;
            using Polygons         = std::vector<Polygon>;
            
            using Group_indices  = std::vector< std::vector<size_t> >;
            using Groups         = std::vector<Vector_3>;
            using Size_t_indices = std::vector<size_t>; 

            Level_of_detail_building_envelope_estimator_triangles(const Input &input, const FT ground_height) :
            m_input(input), 
            m_strategy(input), 
            m_alpha(-FT(1)),
            m_ground_height(ground_height),
            m_big_value(FT(100000000000000)),
            m_min_height( m_big_value),
            m_max_height(-m_big_value),
            m_use_points_based_height(true),
            m_angle_threshold(FT(5)),
            m_use_min_scale(false) { }

            void estimate(Buildings &buildings) {
                
                if (buildings.size() == 0) return;
				for (Building_iterator bit = buildings.begin(); bit != buildings.end(); ++bit) {

                    auto &building = bit->second;
					if (building.is_valid) process_building(building);
                }
            }

            void set_alpha(const FT new_value) {
                assert(new_value > FT(0));
                m_alpha = new_value;
            }

            bool is_face_based() const {
				return false;
			}

        private:
            const Input &m_input;
            Strategy     m_strategy;
            FT           m_alpha;
            
            const FT m_ground_height;
            const FT m_big_value;

            FT m_min_height;
            FT m_max_height;

            const bool m_use_points_based_height;
            const FT   m_angle_threshold;
            const bool m_use_min_scale;

            void process_building(Building &building) {
                
                m_min_height =  m_big_value;
                m_max_height = -m_big_value;

                process_initial_roofs(building);
                find_envelope(building);
            }

            void process_initial_roofs(Building &building) {
                
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
                assert(indices.size() > 2);

                Plane_3 plane;
                fit_plane_to_roof_points(indices, plane);

                Points_3 points;   
                project_points_onto_plane(indices, plane, points);
                update_min_max_heights(points);

                m_strategy.set_alpha(m_alpha);
                m_strategy.estimate_roof(points, plane, building);
            }

            void fit_plane_to_roof_points(const Indices &indices, Plane_3 &plane) const {
                
                assert(indices.size() > 2);
                double bx = 0.0, by = 0.0, bz = 0.0;

                std::vector<Point_3ft> points(indices.size());
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

				Plane_3ft tmp_plane;
                Point_3ft centroid = Point_3ft(bx, by, bz);

				CGAL::linear_least_squares_fitting_3(points.begin(), points.end(), tmp_plane, centroid, CGAL::Dimension_tag<0>(), Local_kernel(), Diagonalize_traits());
				plane = Plane_3(static_cast<FT>(tmp_plane.a()), static_cast<FT>(tmp_plane.b()), static_cast<FT>(tmp_plane.c()), static_cast<FT>(tmp_plane.d()));
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

            void update_min_max_heights(const Points_3 &points) {
                
                FT minz = m_big_value, maxz = -m_big_value;
                for (size_t i = 0; i < points.size(); ++i) {

                    const Point_3 &p = points[i];

                    minz = CGAL::min(minz, p.z());
                    maxz = CGAL::max(maxz, p.z());
                }
                
                m_min_height = CGAL::min(m_min_height, minz);
                m_max_height = CGAL::max(m_max_height, maxz);
            }

            void find_envelope(Building &building) const {

                Data_triangles triangles;
                create_envelope_input(building, triangles);
                
                Env_diagram diagram;
                compute_envelope(triangles, diagram);

                if (m_use_points_based_height) update_roofs_using_points(diagram, triangles, building);
                else update_roofs(diagram, triangles, building);

                set_building_roofs_min_height(building);
            }

            void create_envelope_input(const Building &building, Data_triangles &triangles) const {
                
                size_t index = 0;
                triangles.clear();
                
                add_roofs_triangles(building, triangles, index);
                add_walls_triangles(building, triangles, index);
            }

            void add_roofs_triangles(const Building &building, Data_triangles &triangles, size_t &index) const {

                const Roofs &roofs = building.roofs;
                for (size_t i = 0; i < roofs.size(); ++i) {
                    
                    Boundary roof = roofs[i].boundary;
                    assert(roof.size() == 4);

                    scale_roof(roof);
                    add_two_triangles(roof, building.color, triangles, index);
                }
            }

            void add_walls_triangles(const Building &building, Data_triangles &triangles, size_t &index) const {

                assert(!building.is_oriented);
                const auto &boundary = building.boundaries[0];

                Boundary wall(4);
                const FT height = building.height;

                for (size_t i = 0; i < boundary.size(); i += 2) {
                    const size_t ip = (i + 1) % boundary.size();

                    const Point_2 &v1 = boundary[i]->point();
                    const Point_2 &v2 = boundary[ip]->point();

                    wall[0] = Point_3(v1.x(), v1.y(), FT(0));
                    wall[1] = Point_3(v2.x(), v2.y(), FT(0));
                    wall[2] = Point_3(v2.x(), v2.y(), height);
                    wall[3] = Point_3(v1.x(), v1.y(), height);

                    add_two_triangles(wall, building.color, triangles, index);
                }
            }

            void add_two_triangles(const Boundary &boundary, const Color &color, Data_triangles &triangles, size_t &index) const {
                assert(boundary.size() == 4);
                
                const Point_3 &p1 = boundary[0];
                const Point_3 &p2 = boundary[1];
                const Point_3 &p3 = boundary[2];
                const Point_3 &p4 = boundary[3];

                const Exact_point ep1 = Exact_point(CGAL::to_double(p1.x()), CGAL::to_double(p1.y()), CGAL::to_double(p1.z()));
                const Exact_point ep2 = Exact_point(CGAL::to_double(p2.x()), CGAL::to_double(p2.y()), CGAL::to_double(p2.z()));
                const Exact_point ep3 = Exact_point(CGAL::to_double(p3.x()), CGAL::to_double(p3.y()), CGAL::to_double(p3.z()));
                const Exact_point ep4 = Exact_point(CGAL::to_double(p4.x()), CGAL::to_double(p4.y()), CGAL::to_double(p4.z()));

                const Exact_triangle tri1 = Exact_triangle(ep1, ep2, ep3);
                const Exact_triangle tri2 = Exact_triangle(ep3, ep4, ep1);

                Data data1, data2;
                data1.color = color; data1.index = index; ++index;
                data2.color = color; data2.index = index; ++index;

                const Data_triangle data_tri1 = Data_triangle(tri1, data1);
                const Data_triangle data_tri2 = Data_triangle(tri2, data2);

                triangles.push_back(data_tri1);
                triangles.push_back(data_tri2);
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
                
                if (m_use_min_scale)
                    return compute_min_scale(points);
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

            void compute_envelope(const Data_triangles &triangles, Env_diagram &diagram) const {
                CGAL::lower_envelope_3(triangles.begin(), triangles.end(), diagram);
            }

            void set_building_roofs_min_height(Building &building) const {
                
                FT roofs_min_height = m_big_value;
                for (size_t i = 0; i < building.roofs.size(); ++i) {
                    
                    const auto &boundary = building.roofs[i].boundary;
                    for (size_t j = 0; j < boundary.size(); ++j) {
                        
                        const Point_3 &p = boundary[j];
                        roofs_min_height = CGAL::min(roofs_min_height, p.z());
                    }
                }

                building.roofs_min_height = roofs_min_height - m_ground_height;
            }


            //////////////////
            // TEMPORARY CODE!
            

            void update_roofs_using_points(const Env_diagram &diagram, const Data_triangles &triangles, Building &building) const {
                building.clear_roofs();
                
                Roof roof;
                Boundary &boundary = roof.boundary;

                for (Face_iterator fit = diagram.faces_begin(); fit != diagram.faces_end(); ++fit) {
                    
                    boundary.clear();
                    if (fit->is_unbounded()) continue;

                    Halfedge_circulator ccb = fit->outer_ccb();
                    do {

                        const double x = CGAL::to_double(ccb->target()->point().x());
                        const double y = CGAL::to_double(ccb->target()->point().y());

                        const Point_2 p = Point_2(static_cast<FT>(x), static_cast<FT>(y));
                        const FT z      = get_height_tmp(p, fit, triangles);               
                        
                        // const FT z = get_height_using_points(p, building);

                        boundary.push_back(Point_3(p.x(), p.y(), z));
                        ++ccb;

                    } while (ccb != fit->outer_ccb());

                    if (is_valid_roof_face(building, boundary))
                        building.roofs.push_back(roof);
                }
            }

            FT get_height_tmp(const Point_2 &p, const Face_iterator &fit, const Data_triangles &triangles) const {

                for (Surface_iterator sit = fit->surfaces_begin(); sit != fit->surfaces_end(); ++sit) {
                    const size_t index = sit->data().index;

                    const Exact_point &a = triangles[index].vertex(0);
					const Exact_point &b = triangles[index].vertex(1);
					const Exact_point &c = triangles[index].vertex(2);

                    Point_3 p1, p2, p3;
                    get_corrected_points(a, b, c, p1, p2, p3);

                    const Plane_3 plane = Plane_3(p1, p2, p3);
                    return intersect_with_plane(p, plane);
                }

                return FT(0); // be careful here, it may cause problems with negative heights, see the half-tile for example
            }

            FT intersect_with_plane(const Point_2 &query, const Plane_3 &plane) const {
                
                const Point_3 p = Point_3(query.x(), query.y(), FT(0));
                const Point_3 q = Point_3(query.x(), query.y(), FT(1));

                const Line_3 line = Line_3(p, q);

				typename CGAL::cpp11::result_of<Intersect(Line_3, Plane_3)>::type result = intersection(line, plane);
                const Point_3 r = boost::get<Point_3>(*result);

				return r.z();
            }

            FT get_height_using_points(const Point_2 &p, const Building &building) const {
                
                const FT big_value      = FT(100000000000000);
                const FT default_height = m_ground_height + building.height;
                
                FT mindist = big_value;
                FT height  = default_height;

                update_height_using_points(p, building, mindist, height);
                return height;
            }

            void update_height_using_points(const Point_2 &p, const Building &building, FT &mindist, FT &height) const {

                for (size_t i = 0; i < building.shapes.size(); ++i) {
                    const auto &indices = building.shapes[i];

                    for (size_t j = 0; j < indices.size(); ++j) {
                        const Point_3 &tmp = m_input.point(indices[j]);

                        const Exact_point q = Exact_point(CGAL::to_double(tmp.x()), CGAL::to_double(tmp.y()), CGAL::to_double(tmp.z()));
                        update_height(p, q, mindist, height);
                    }
                }
            }

            void update_height(const Point_2 &p, const Exact_point &q, FT &mindist, FT &height) const {
                
                const FT x = static_cast<FT>(CGAL::to_double(q.x()));
                const FT y = static_cast<FT>(CGAL::to_double(q.y()));

                const FT dist = squared_distance_2(p, Point_2(x, y));
                if (dist < mindist) {
                        
                    mindist = dist;
                    height  = static_cast<FT>(CGAL::to_double(q.z()));
                }
            }


            ////////////
            // NEW CODE!


            void update_roofs(const Env_diagram &diagram, const Data_triangles &triangles, Building &building) const {
                
                Roof roof;
                building.clear_roofs();
                
                for (Face_iterator fit = diagram.faces_begin(); fit != diagram.faces_end(); ++fit) {
                    if (fit->is_unbounded()) continue;

                    if (!process_roof_face(fit, triangles, building, roof)) continue;
                    building.roofs.push_back(roof);
                }
                
                postprocess_roofs(building);
                if (building.roofs.size() == 0) building.is_valid = false;
            }

            bool process_roof_face(const Face_iterator &fit, const Data_triangles &triangles, const Building &building, Roof &roof) const {

                Boundary &boundary = roof.boundary;
                boundary.clear();

                Vector_3 target_normal;
                const bool success = compute_face_target_normal(fit, triangles, target_normal);
                if (!success) return false;

                get_face_boundary(fit, building, boundary);
                if (boundary.size() < 3) return false;

                if (!is_valid_roof_face(building, boundary)) return false;
                roof.tmp = boundary;

                const Vector_3 source_normal = Vector_3(FT(0), FT(0), FT(1));
                rotate_face(source_normal, target_normal, boundary);

                return true;
            }

            bool compute_face_target_normal(const Face_iterator &fit, const Data_triangles &triangles, Vector_3 &target_normal) const {

                FT wrong = FT(0), count = FT(0);
                FT x = FT(0), y = FT(0), z = FT(0);

                for (Surface_iterator sit = fit->surfaces_begin(); sit != fit->surfaces_end(); ++sit, count += FT(1)) {
                    const size_t index = sit->data().index;

                    const Exact_point &a = triangles[index].vertex(0);
					const Exact_point &b = triangles[index].vertex(1);
					const Exact_point &c = triangles[index].vertex(2);

                    Point_3 p1, p2, p3;
                    get_corrected_points(a, b, c, p1, p2, p3);

                    Vector_3 normal;
                    compute_triangle_normal(p1,p2, p3, normal);

                    if (!is_valid_normal(normal)) {
                     
                        wrong += FT(1);
                        continue;
                    }

                    x += normal.x();
                    y += normal.y();
                    z += normal.z();
                }

                if (wrong == count) return false;
                if (count == FT(0)) return false;
                
                x /= count;
                y /= count;
                z /= count;

                target_normal = Vector_3(x, y, z);
                return true;
            }

            void get_corrected_points(const Exact_point &a, const Exact_point &b, const Exact_point &c, Point_3 &p1, Point_3 &p2, Point_3 &p3) const {

                p1 = Point_3(static_cast<FT>(CGAL::to_double(a.x())), static_cast<FT>(CGAL::to_double(a.y())), static_cast<FT>(CGAL::to_double(a.z())));
                p2 = Point_3(static_cast<FT>(CGAL::to_double(b.x())), static_cast<FT>(CGAL::to_double(b.y())), static_cast<FT>(CGAL::to_double(b.z())));
                p3 = Point_3(static_cast<FT>(CGAL::to_double(c.x())), static_cast<FT>(CGAL::to_double(c.y())), static_cast<FT>(CGAL::to_double(c.z())));
            }

            void compute_triangle_normal(const Point_3 &p1, const Point_3 &p2, const Point_3 &p3, Vector_3 &normal) const {

                const Vector_3 v1 = Vector_3(p1, p2);
                const Vector_3 v2 = Vector_3(p1, p3);

                normal = cross_product(v1, v2);
                normalize(normal);
            }

            bool is_valid_normal(const Vector_3 &normal) const {
                
                const Vector_3 base_normal = Vector_3(FT(0), FT(0), FT(1));

                FT angle;
                Vector_3 stub;

                compute_angle_and_axis(base_normal, normal, angle, stub);
                const FT angle_diff = FT(90) - CGAL::abs(angle);

                if (angle_diff < m_angle_threshold) return false;
                return true;
            }

            void get_face_boundary(const Face_iterator &fit, const Building &building, Boundary &boundary) const {
                
                Halfedge_circulator ccb = fit->outer_ccb();
                do {

                    const FT x = static_cast<FT>(CGAL::to_double(ccb->target()->point().x()));
                    const FT y = static_cast<FT>(CGAL::to_double(ccb->target()->point().y()));
                    const FT z = compute_min_roof_height(building);

                    boundary.push_back(Point_3(x, y, z));
                    ++ccb;

                } while (ccb != fit->outer_ccb());
            }

            FT compute_min_roof_height(const Building &building) const {

                FT min_height = m_big_value;
                for (size_t i = 0; i < building.shapes.size(); ++i) {
                    
                    const auto &indices = building.shapes[i];
                    for (size_t j = 0; j < indices.size(); ++j) {
                        
                        const Point_3 &p = m_input.point(indices[j]);
                        min_height = CGAL::min(min_height, p.z());
                    }
                }

                return min_height;
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

            void rotate_face(const Vector_3 &source_normal, const Vector_3 &target_normal, Boundary &boundary) const {
                
                Point_3 b;
                compute_barycentre(boundary, b);

                FT angle;
                Vector_3 axis;
                compute_angle_and_axis(source_normal, target_normal, angle, axis);

                for (size_t i = 0; i < boundary.size(); ++i) {   
                    Point_3 &p = boundary[i];

                    Point_3 tmp = Point_3(p.x() - b.x(), p.y() - b.y(), p.z() - b.z());
                    rotate_point(angle, axis, tmp);
                    p = Point_3(tmp.x() + b.x(), tmp.y() + b.y(), tmp.z() + b.z());
                }
            }

            void compute_angle_and_axis(const Vector_3 &m, const Vector_3 &n, FT &angle, Vector_3 &axis) const {
				
				const auto cross = cross_product(m, n);
				const FT length  = static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_length(cross))));
				const FT dot     = dot_product(m, n);

				angle = static_cast<FT>(std::atan2(CGAL::to_double(length), CGAL::to_double(dot)));
				axis  = cross / length;

                const FT half_pi = static_cast<FT>(CGAL_PI) / FT(2);
                if (angle > half_pi) {
                    
                    angle = static_cast<FT>(CGAL_PI) - angle;
                    axis = -axis;
                }
			}

            void rotate_point(const FT angle, const Vector_3 &axis, Point_3 &p) const {

				const double tmp_angle = CGAL::to_double(angle);

				const FT c = static_cast<FT>(std::cos(tmp_angle));
				const FT s = static_cast<FT>(std::sin(tmp_angle));

				const FT C = FT(1) - c;

				const FT x = axis.x();
				const FT y = axis.y();
				const FT z = axis.z();

				p = Point_3((x * x * C + c)     * p.x() + (x * y * C - z * s) * p.y() + (x * z * C + y * s) * p.z(),
					  		(y * x * C + z * s) * p.x() + (y * y * C + c)     * p.y() + (y * z * C - x * s) * p.z(),
					  		(z * x * C - y * s) * p.x() + (z * y * C + x * s) * p.y() + (z * z * C + c)     * p.z());
			}

            void postprocess_roofs(Building &building) const {

                Groups        groups; 
                Group_indices group_indices;
                
                create_groups(building, groups, group_indices);
                rotate_groups(groups, group_indices, building);
                translate_groups(groups, group_indices, building);
            }

            void create_groups(const Building &building, Groups &groups, Group_indices &group_indices) const {

                const Roofs &roofs = building.roofs;
                const size_t num_roofs = roofs.size();

                for (size_t i = 0; i < num_roofs; ++i) {        
                    const Boundary &boundary = roofs[i].boundary;

                    assert(boundary.size() > 2);

                    const Point_3 &p1 = boundary[0];
                    const Point_3 &p2 = boundary[1];
                    const Point_3 &p3 = boundary[2];

                    Vector_3 roof_normal;
                    compute_triangle_normal(p1, p2, p3, roof_normal);

                    const int group_index = find_group(groups, roof_normal);
                    if (group_index < 0) {
                        
                        std::vector<size_t> tmp(1);
                        tmp[0] = i;

                        groups.push_back(roof_normal);
                        group_indices.push_back(tmp);

                        continue;
                    } 
                    group_indices[group_index].push_back(i);
                }
            }

            int find_group(const Groups &groups, const Vector_3 &n) const {

                const FT eps = FT(1) / FT(1000000);
                for (size_t i = 0; i < groups.size(); ++i) {
                    const Vector_3 &m = groups[i];
                    
                    const FT dx = CGAL::abs(m.x() - n.x());
                    const FT dy = CGAL::abs(m.y() - n.y());
                    const FT dz = CGAL::abs(m.z() - n.z());

                    if (dx < eps && dy < eps && dz < eps) return static_cast<int>(i);
                }
                return -1;
            }

            void rotate_groups(const Groups &groups, const Group_indices &group_indices, Building &building) const {

                for (size_t i = 0; i < group_indices.size(); ++i)
                    rotate_group(groups[i], group_indices[i], building);
            }

            void rotate_group(const Vector_3 &target_normal, const Size_t_indices &indices, Building &building) const {

                Point_3 b;
                compute_group_barycentre(building, indices, b);
                const Vector_3 source_normal = Vector_3(FT(0), FT(0), FT(1));

                FT angle;
                Vector_3 axis;
                compute_angle_and_axis(source_normal, target_normal, angle, axis);

                Roofs &roofs = building.roofs;
                for (size_t i = 0; i < indices.size(); ++i) {
                    
                    const size_t index = indices[i];
                    Boundary &boundary = roofs[index].tmp;

                    for (size_t j = 0; j < boundary.size(); ++j) {
                        Point_3 p = boundary[j];

                        Point_3 tmp = Point_3(p.x() - b.x(), p.y() - b.y(), p.z() - b.z());
                        rotate_point(angle, axis, tmp);
                        p = Point_3(tmp.x() + b.x(), tmp.y() + b.y(), tmp.z() + b.z());

                        roofs[index].boundary[j] = p;
                    }
                }
            }

            void compute_group_barycentre(const Building &building, const Size_t_indices &indices, Point_3 &barycentre) const {
                
                FT x = FT(0), y = FT(0), z = FT(0);

                const Roofs &roofs = building.roofs;
                for (size_t i = 0; i < indices.size(); ++i) {
                    
                    const size_t index = indices[i];
                    const Boundary &boundary = roofs[index].tmp;

                    Point_3 b;
                    compute_barycentre(boundary, b);

                    x += b.x();
                    y += b.y();
                    z += b.z();
                }

                x /= static_cast<FT>(indices.size());
                y /= static_cast<FT>(indices.size());
                z /= static_cast<FT>(indices.size());

                barycentre = Point_3(x, y, z);
            }

            void translate_groups(const Groups &groups, const Group_indices &group_indices, Building &building) const {

                const FT min_roofs_height = compute_min_roof_height(building);

                for (size_t i = 0; i < group_indices.size(); ++i)
                    translate_group(groups[i], group_indices[i], min_roofs_height, building);
            }

            void translate_group(const Vector_3 &, const Size_t_indices &indices, const FT min_roofs_height, Building &building) const {

                const FT diff = compute_group_translation(indices, min_roofs_height, building);

                Roofs &roofs = building.roofs;
                for (size_t i = 0; i < indices.size(); ++i) {
                    
                    const size_t index = indices[i];
                    Boundary &boundary = roofs[index].boundary;

                    for (size_t j = 0; j < boundary.size(); ++j)
                        boundary[j] = Point_3(boundary[j].x(), boundary[j].y(), boundary[j].z() - diff);
                }
            }

            FT compute_group_translation(const Size_t_indices &indices, const FT min_roofs_height, const Building &building) const {

                FT minz = m_big_value;
                const Roofs &roofs = building.roofs;

                for (size_t i = 0; i < indices.size(); ++i) {    
                    const size_t index = indices[i];

                    const Boundary &boundary = roofs[index].boundary;
                    for (size_t j = 0; j < boundary.size(); ++j) minz = CGAL::min(minz, boundary[j].z());
                }

                return minz - min_roofs_height;
            }
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_ENVELOPE_ESTIMATOR_TRIANGLES_H