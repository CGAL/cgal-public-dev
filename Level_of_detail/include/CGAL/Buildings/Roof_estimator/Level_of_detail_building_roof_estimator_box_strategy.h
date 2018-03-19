#ifndef CGAL_LEVEL_OF_DETAIL_BUILDING_ROOF_ESTIMATOR_BOX_STRATEGY_H
#define CGAL_LEVEL_OF_DETAIL_BUILDING_ROOF_ESTIMATOR_BOX_STRATEGY_H

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
#include <CGAL/Level_of_detail_enum.h>

namespace CGAL {

	namespace LOD {

		// Main class.
		template<class KernelTraits, class ContainerInput, class BuildingInput>
		class Level_of_detail_building_roof_estimator_box_strategy {
            
        public:
            typedef KernelTraits   Kernel;
            typedef ContainerInput Input;
            typedef BuildingInput  Building;

            typename Kernel::Compute_squared_distance_2 	  squared_distance;
			typename Kernel::Compute_squared_length_3 		  squared_length;
			typename Kernel::Compute_scalar_product_3 		  dot_product;
			typename Kernel::Construct_cross_product_vector_3 cross_product;

            using FT       = typename Kernel::FT;
            using Point_3  = typename Kernel::Point_3;
            using Plane_3  = typename Kernel::Plane_3;
            using Line_3   = typename Kernel::Line_3;
            using Vector_3 = typename Kernel::Vector_3;
            using Points   = std::vector<Point_3>;

            using Local_kernel = CGAL::Simple_cartesian<double>;
			using Point_3ft    = Local_kernel::Point_3;
			using Line_3ft     = Local_kernel::Line_3;

            using Roof = typename Building::Roof;

            Level_of_detail_building_roof_estimator_box_strategy(const Input &input) : m_input(input) { }

            void estimate_roof(Points &roof_points, const Plane_3 &plane, Building &building) const {
                if (roof_points.size() < 2) return;

                Vector_3 roof_direction;
                estimate_roof_direction(roof_points, roof_direction);

                Vector_3 y_direction;
                set_y_direction(y_direction);
                
				Vector_3 roof_normal;
				set_plane_normal(plane, roof_normal);

				Vector_3 ground_normal;
				set_ground_normal(ground_normal);

				FT angle1; Vector_3 axis;
				compute_angle_and_axis(roof_normal, ground_normal, angle1, axis);
							
				FT angle2; Vector_3 stub;
				compute_angle_and_axis(roof_direction, y_direction, angle2, stub);

                Point_3 barycentre;
                set_roof_barycentre(roof_points, barycentre);

                estimate_building_roof(angle1, angle2, axis, barycentre, roof_points, building);
            }

			void set_alpha(const FT) { }

			bool is_face_based() const {
				return false;
			}

        private:
            const Input &m_input;

            void estimate_roof_direction(const Points &points, Vector_3 &direction) const {
                assert(points.size() > 1);

				std::vector<Point_3ft> tmp_points(points.size());
				for (size_t i = 0; i < points.size(); ++i) {
					
					const Point_3 &p = points[i];

					const double x = CGAL::to_double(p.x());
					const double y = CGAL::to_double(p.y());
					const double z = CGAL::to_double(p.z());

					tmp_points[i] = Point_3ft(x, y, z);
				}

                Line_3ft tmp_line;
				CGAL::linear_least_squares_fitting_3(tmp_points.begin(), tmp_points.end(), tmp_line, CGAL::Dimension_tag<0>());

				const Line_3 line = Line_3(
                Point_3(static_cast<FT>(tmp_line.point(0).x()), static_cast<FT>(tmp_line.point(0).y()), static_cast<FT>(tmp_line.point(0).z())), 
				Point_3(static_cast<FT>(tmp_line.point(1).x()), static_cast<FT>(tmp_line.point(1).y()), static_cast<FT>(tmp_line.point(1).z())));

                direction = line.to_vector();
            }

            void set_y_direction(Vector_3 &direction) const {
                direction = Vector_3(FT(0), FT(1), FT(0));
            }

            void set_plane_normal(const Plane_3 &plane, Vector_3 &m) const {
				m = plane.orthogonal_vector();
			}

			void set_ground_normal(Vector_3 &n) const {
				const Plane_3 ground = Plane_3(FT(0), FT(0), FT(1), FT(0));				
				n = ground.orthogonal_vector();
			}

            void compute_angle_and_axis(const Vector_3 &m, const Vector_3 &n, FT &angle, Vector_3 &axis) const {
				
				const auto cross = cross_product(m, n);
				const FT length  = static_cast<FT>(CGAL::sqrt(CGAL::to_double(squared_length(cross))));
				const FT dot     = dot_product(m, n);

				angle = static_cast<FT>(std::atan2(CGAL::to_double(length), CGAL::to_double(dot)));
				axis = cross / length;
			}

            void set_roof_barycentre(const Points &points, Point_3 &barycentre) const {

				FT bx = 0, by = 0;
				for (size_t i = 0; i < points.size(); ++i) {
							
					const Point_3 &p = points[i];

					bx += p.x();
					by += p.y();
				}

				bx /= static_cast<FT>(points.size());
                by /= static_cast<FT>(points.size());

                barycentre = Point_3(bx, by, FT(0));
            }

            void estimate_building_roof(const FT angle1, const FT angle2, const Vector_3 &axis, const Point_3 &bary,
            Points &roof_points, Building &building) const {

                const FT big_value = FT(100000000000000);

				FT minx =  big_value, miny =  big_value;
				FT maxx = -big_value, maxy = -big_value;

				FT z;
				for (size_t i = 0; i < roof_points.size(); ++i) {
					Point_3 &p = roof_points[i];

					rotate_point_1(angle1, axis, p);	
					rotate_point_2(angle2, bary, p);

					minx = CGAL::min(minx, p.x());
					miny = CGAL::min(miny, p.y());

				    maxx = CGAL::max(maxx, p.x());
					maxy = CGAL::max(maxy, p.y());

					z = p.z();
				}

                Points boundary(4);

                boundary[0] = Point_3(minx, miny, z);
				boundary[1] = Point_3(maxx, miny, z);
				boundary[2] = Point_3(maxx, maxy, z);
				boundary[3] = Point_3(minx, maxy, z);

				rotate_point_2(-angle2, bary, boundary[0]);
				rotate_point_1(-angle1, axis, boundary[0]);

                rotate_point_2(-angle2, bary, boundary[1]);
				rotate_point_1(-angle1, axis, boundary[1]);

                rotate_point_2(-angle2, bary, boundary[2]);
				rotate_point_1(-angle1, axis, boundary[2]);

                rotate_point_2(-angle2, bary, boundary[3]);
				rotate_point_1(-angle1, axis, boundary[3]);

				if (std::isnan(CGAL::to_double(boundary[0].x())) ||
					std::isnan(CGAL::to_double(boundary[0].y())) ||
					std::isnan(CGAL::to_double(boundary[0].z()))  ) return;

                if (std::isnan(CGAL::to_double(boundary[1].x())) ||
					std::isnan(CGAL::to_double(boundary[1].y())) ||
					std::isnan(CGAL::to_double(boundary[1].z()))  ) return;

                if (std::isnan(CGAL::to_double(boundary[2].x())) ||
					std::isnan(CGAL::to_double(boundary[2].y())) ||
					std::isnan(CGAL::to_double(boundary[2].z()))  ) return;

                if (std::isnan(CGAL::to_double(boundary[3].x())) ||
					std::isnan(CGAL::to_double(boundary[3].y())) ||
					std::isnan(CGAL::to_double(boundary[3].z()))  ) return;

                Roof roof;

                roof.boundary = boundary;
                building.roofs.push_back(roof);
            }

			void rotate_point_1(const FT angle, const Vector_3 &axis, Point_3 &p) const {

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

            void rotate_point_2(const FT angle, const Point_3 &bary, Point_3 &p) const {

				FT x = p.x();
				FT y = p.y();

				x -= bary.x();
				y -= bary.y();

				p = Point_3(x, y, p.z());

                const double tmp_angle = CGAL::to_double(angle);

                const FT c = static_cast<FT>(std::cos(tmp_angle));
				const FT s = static_cast<FT>(std::sin(tmp_angle));

				x = p.x() * c - p.y() * s;
				y = p.y() * c + p.x() * s;

				x += bary.x();
				y += bary.y();

				p = Point_3(x, y, p.z());
			} 
        };
    }
}

#endif // CGAL_LEVEL_OF_DETAIL_BUILDING_ROOF_ESTIMATOR_BOX_STRATEGY_H